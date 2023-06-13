//
// Created by ixiaohu on 2023/5/6.
//

#include "comp_seeding.h"

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>

#include "../Matching/SimplePgMatcher.h"
#include "../Utils/helper.h"
#include "../FM_index/bntseq.h"
#include "../bwalib/bwa.h"
#include "../bwalib/ksw.h"
#include "../cstl/kvec.h"
#include "../cstl/kthread.h"
#include "../cstl/kbtree.h"
#include "../cstl/ksort.h"
#include "../cstl/kstring.h"
#include "../bwalib/utils.h"

/********************************
 * Bookmark 1: Seeding with SST *
 ********************************/
CompAligner::CompAligner() {
	uint64_t stamp = __rdtsc(); sleep(1); cpu_frequency = __rdtsc() - stamp;
	opt = mem_opt_init();
}

void CompAligner::load_index(const char *fn) {
	bwa_idx = bwa_idx_load_from_shm(fn);
	if (bwa_idx == nullptr) {
		bwa_idx = bwa_idx_load(fn, BWA_IDX_ALL);
		if (bwa_idx == nullptr) {
			fprintf(stderr, "Load the FM-index `%s` failed\n", fn);
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "Load the FM-index from disk\n");
		}
	} else {
		fprintf(stderr, "Load the FM-index from shared memory\n");
	}
	bns = bwa_idx->bns;
	bwt = bwa_idx->bwt;
	pac = bwa_idx->pac;
}

int CompAligner::collect_mem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux &aux) {
	aux.super_mem.clear();
	if (seq[pivot] > 3) return pivot + 1;
	bwtintv_t ik, next[4];
	int node_id = aux.forward_sst->query_forward_child(0, seq[pivot]);
	ik = aux.forward_sst->get_intv(node_id);
	ik.info = pivot + 1;
	aux.prev_intv.clear();
	int ret_pivot = len;
	for (int i = pivot + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = aux.forward_sst->query_forward_child(node_id, c);
			next[c] = aux.forward_sst->get_intv(node_id);
			if (next[c].x[2] != ik.x[2]) {
				aux.prev_intv.push_back(ik);
				if (next[c].x[2] < min_hits) {
					ret_pivot = i;
					break;
				}
			}
			ik = next[c];
			ik.info = i + 1;
		} else {
			aux.prev_intv.push_back(ik);
			ret_pivot = i + 1;
			break;
		}
	}
	if (ret_pivot == len) aux.prev_intv.push_back(ik);
	if (pivot == 0) { // Quickly return full-length matched reads
		aux.super_mem.push_back(aux.prev_intv.back());
		return ret_pivot;
	}
	std::reverse(aux.prev_intv.begin(), aux.prev_intv.end());

	// Collect node indexes for LEPs in backward SST
	for (auto &p : aux.prev_intv) {
		node_id = 0;
		for (int j = (int)p.info - 1; j >= pivot + 1; j--) {
			node_id = aux.backward_sst->add_empty_child(node_id, seq[j]);
		}
		node_id = aux.backward_sst->add_lep_child(node_id, seq[pivot], p.x);
		p.info |= (1UL * node_id << 32);
	}

	for (int i = pivot - 1; i >= -1; i--) {
		int c = (i == -1) ?4 :seq[i];
		aux.curr_intv.clear();
		for (int j = 0; j < aux.prev_intv.size(); j++) {
			ik = aux.prev_intv[j];
			node_id = ik.info >> 32;
			if (c < 4) {
				node_id = aux.backward_sst->query_backward_child(node_id, c);
				next[c] = aux.backward_sst->get_intv(node_id);
			}
			if (c > 3 or next[c].x[2] < min_hits) {
				if (aux.super_mem.empty() or i + 1 < aux.super_mem.back().info >> 32) {
					ik.info = (1UL * (i + 1) << 32) | (int32_t)ik.info;
					aux.super_mem.push_back(ik);
				}
			} else if (aux.curr_intv.empty() or next[c].x[2] != aux.curr_intv.back().x[2]) {
				next[c].info = (1UL * node_id << 32) | (int32_t)ik.info;
				aux.curr_intv.push_back(next[c]);
			}
		}
		if (aux.curr_intv.empty()) break;
		std::swap(aux.prev_intv, aux.curr_intv);
	}
	return ret_pivot;
}

int CompAligner::tem_forward_sst(const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux &aux) {
	if (seq[start] > 3) return start + 1;
	memset(mem, 0, sizeof(bwtintv_t));
	int node_id = aux.forward_sst->query_forward_child(0, seq[start]);
	bwtintv_t ik = aux.forward_sst->get_intv(node_id);
	for (int i = start + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = aux.forward_sst->query_forward_child(node_id, c);
			ik = aux.forward_sst->get_intv(node_id);
			if (ik.x[2] < opt->max_mem_intv && i - start >= opt->min_seed_len ) {
				*mem = ik;
				mem->info = (1UL * start << 32) | (i + 1);
				return i + 1;
			}
		} else return i + 1;
	}
	return len;
}

static inline bool mem_cmp(const bwtintv_t &a, const bwtintv_t &b) { return a.info < b.info; }

static inline int mem_beg(const bwtintv_t &a) { return a.info >> 32; }

static inline int mem_end(const bwtintv_t &a) { return (int)a.info; }

static inline int mem_len(const bwtintv_t &a) { return mem_end(a) - mem_beg(a); }

static std::string mem_str(const bwtintv_t &a) {
	char buf[256];
	sprintf(buf, "[%d, %d) [%ld, %ld, %ld]", mem_beg(a), mem_end(a), a.x[0], a.x[1], a.x[2]);
	return std::string(buf);
}

/***********************
 * Bookmark 2 Chaining *
 ***********************/

#define chain_cmp(a, b) (((b).anchor < (a).anchor) - ((a).anchor < (b).anchor))
KBTREE_INIT(chn, seed_chain, chain_cmp)

bool CompAligner::add_seed_to_chain(seed_chain *c, const seed_hit &s) {
	if (s.rid != c->rid) return false; // Seed is not on same chromosome with chain
	if (s.qbeg >= c->que_beg() and s.qbeg + s.len <= c->que_end()  and
		s.rbeg >= c->ref_beg() and s.rbeg + s.len <= c->ref_end()) {
		return true; // Seed in contained in chain
	}
	const auto &last = c->seeds[c->n - 1];
	if ((last.rbeg < bns->l_pac or c->ref_beg() < bns->l_pac) and s.rbeg >= bns->l_pac) {
		return false; // Seed is on different strand from chain
	}
	int64_t x = s.qbeg - last.qbeg; // Non-negative; guaranteed by sorting matches
	int64_t y = s.rbeg - last.rbeg;
	if (// Chained seeds should be co-linear on both read and reference
		y >= 0 and
		// Minimum insertion/deletions should be within DP matrix bandwidth
		std::abs(x - y) <= opt->w and
		// Gap between chain and seed on both read and reference is limited
		x - last.len < opt->max_chain_gap and y - last.len < opt->max_chain_gap
	) {
		c->push_back(s); // Grow the chain
		return true;
	} else {
		return false;
	}
}

void CompAligner::print_chains_to(const std::vector<seed_chain> &chains, kstring_t *s) {
	for (int i = 0; i < chains.size(); i++) {
		const auto &c = chains[i];
		ksprintf(s, "* Found CHAIN(%d): n=%d; weight=%d", i, c.n, c.calc_weight());
		for (int j = 0; j < c.n; j++) {
			bwtint_t pos; int is_rev;
			pos = bns_depos(bns, c.seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= c.seeds[j].len - 1;
			ksprintf(s, "\t%d;%d;%d,%ld(%s:%c%ld)",
				c.seeds[j].score, c.seeds[j].len,
				c.seeds[j].qbeg, (long)c.seeds[j].rbeg,
				bns->anns[c.rid].name, "+-"[is_rev],
				(long)(pos - bns->anns[c.rid].offset) + 1);
		}
		ksprintf(s, "\n");
	}
}

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(flt, seed_chain, flt_lt)

std::vector<seed_chain> CompAligner::chaining(const std::vector<seed_hit> &seed) {
	kbtree_chn_t *tree = kb_init_chn(KB_DEFAULT_SIZE) ;
	for (const auto &s : seed) {
		if (s.rid < 0) continue;
		bool need_new_chain = true;
		seed_chain tmp(s.rbeg), *lower, *upper;
		if (kb_size(tree) > 0) {
			kb_intervalp_chn(tree, &tmp, &lower, &upper);
			if (lower) {
				// Check if the seed could merged into the closet chain
				need_new_chain = not add_seed_to_chain(lower, s);
			}
		}
		if (need_new_chain) {
			tmp.add_first_seed(s);
			tmp.is_alt = bns->anns[s.rid].is_alt;
			kb_putp_chn(tree, &tmp);
		}
	}

	std::vector<seed_chain> chain;
	chain.reserve(kb_size(tree));
	#define traverse_func(p) (chain.push_back(*(p)))
		__kb_traverse(seed_chain, tree, traverse_func);
	#undef traverse_func
	kb_destroy(chn, tree);

	if (chain.empty()) return chain;

	// Drop chains with the weight smaller than the specified minimum
	const int DOMINATE = 3; // Dominating chains with large weight
	const int SHADOWED = 2; // Shadowed by another chain
	const int DISCARD = 0; // Dropped chain
	int kept_chain_n = 0;
	for (auto &c : chain) {
		c.kept = DISCARD;
		c.first_cover = -1;
		c.w = c.calc_weight();
		if (c.w < opt->min_chain_weight) c.destroy();
		else chain[kept_chain_n++] = c;
	}
	chain.resize(kept_chain_n);
	// Sort chains by decreasing weight
	auto *chain_array = chain.data();
	ks_introsort_flt(kept_chain_n, chain_array);
	// Pairwise comparing chains
	std::vector<int> kept_index; // Storing the indexes of kept chains
	chain[0].kept = DOMINATE;
	kept_index.push_back(0);
	for (int i = 1; i < chain.size(); i++) {
		bool masked = false, dropped = false;
		for (auto j : kept_index) { // Check whether the current chain is overlapped by any kept chain
			int beg_max = std::max(chain[i].que_beg(), chain[j].que_beg());
			int end_min = std::min(chain[i].que_end(), chain[j].que_end());
			if (end_min > beg_max and not (chain[j].is_alt and not chain[i].is_alt)) { // Overlap found
				// When the kept chain is alternative and the current chain is primary, don't consider overlap
				int min_len = std::min(
					chain[i].que_end() - chain[i].que_beg(),
					chain[j].que_end() - chain[j].que_beg()
				);
				if ((float)(end_min - beg_max) >= (float)min_len * opt->mask_level
					and min_len < opt->max_chain_gap) {
					masked = true;
					// Record the first or heaviest chain that a kept chain shadows for more accurate MAPQ
					if (chain[j].first_cover < 0) chain[j].first_cover = i;
					if ((float)chain[i].w < (float) chain[j].w * opt->drop_ratio
						and chain[j].w - chain[i].w >= opt->min_seed_len * 2) {
						dropped = true;
						break;
					}
				}
			}
		}
		if (not dropped) {
			chain[i].kept = masked ?SHADOWED :DOMINATE;
			kept_index.push_back(i);
		}
	}
	for (auto i : kept_index) {
		const auto &c = chain[i];
		if (c.first_cover >= 0) {
			chain[c.first_cover].kept = 1; // the function of value 1?
		}
	}
	// Do not extend too many insignificant chains
	int to_extend = 0;
	for (auto &i : kept_index) {
		auto &c = chain[i];
		if (c.kept == DOMINATE) continue;
		to_extend++;
		if (to_extend >= opt->max_chain_extend) c.kept = DISCARD;
	}
	kept_chain_n = 0;
	for (const auto &c : chain) {
		if (c.kept == DISCARD) c.destroy();
		else chain[kept_chain_n++] = c;
	}
	chain.resize(kept_chain_n);
	return chain;
}

void CompAligner::display_profile() {
	thread_aux total;
	for (int i = 0; i < threads_n; i++) total += thr_aux[i];
	auto full_read_match = total.full_read_match;
	auto shortcut = total.shortcut;
	fprintf(stderr, "Input %d reads in total\n", total_reads_count);
	fprintf(stderr, "Perfect Matched Reads: %d (%.2f %%)\n", full_read_match, 100.0 * full_read_match / total_reads_count);
	fprintf(stderr, "Reads go shortcut:     %d (%.2f %%)\n", shortcut, 100.0 * shortcut / total_reads_count);
}

void CompAligner::seed_and_extend(int _start, int _end, int tid) {
	auto &aux = thr_aux[tid];
	aux.forward_sst->clear(); aux.backward_sst->clear();
	_end = std::min(_end, (int) reads.size()); // Out of right boundary happens
	int n = _end - _start;
	long ref_position = -1; // Starting from this position might avoid many BWT-extension for fully matched reads
	for (int i = n - 1; i >= 0; i--) {
		auto &read = reads[_start + i];
		auto *bases = (uint8_t*) read.bases;
		for (int j = 0; j < read.len; j++) { // Convert ACGTN to 01234 if hasn't done so far
			if (bases[j] > 4) bases[j] = nst_nt4_table[bases[j]];
		}

		std::vector<bwtintv_t> &match = aux.match[i]; match.clear();
		bool full_match = false; // Whether this read could be full-length matched
		if (ref_position != -1 and ref_position - read.offset <= read.len / 3) {
			// Lookup forward SST to test if this read can be fully matched
			int node_id = 0;
			for (int j = ref_position - read.offset; j < read.len; j++) {
				if (bases[j] < 4) {
					if (j == ref_position - read.offset) node_id = aux.forward_sst->get_child(node_id, bases[j]);
					else node_id = aux.forward_sst->get_child(node_id, 3 - bases[j]);
				} else node_id = -1;
				if (node_id == -1) break;
			}
			if (node_id != -1) { // matched in forward SST till the end of read
				bwtintv_t ik = aux.forward_sst->get_intv(node_id), ok[4];
				// It is a potential full-match read
				for (int j = ref_position - read.offset - 1; j >= 0; j--) {
					if (bases[j] < 4) {
						bwt_extend(bwt, &ik, ok, 1);
						ik = ok[bases[j]];
					} else ik.x[2] = 0;
					if (ik.x[2] == 0) break;
				}
				if (ik.x[2] > 0) {
					full_match = true;
					aux.shortcut++;
					ik.info = read.len;
					if (mem_len(ik) >= opt->min_seed_len) match.push_back(ik);
				}
			}
		}

		if (not full_match) { // Go to the regular SMEM searching pass
			for (int j = 0; j < read.len; ) {
				j = collect_mem_with_sst(bases, read.len, j, 1, aux);
				for (const auto &m : aux.super_mem) {
					if (mem_len(m) >= opt->min_seed_len)
						match.push_back(m);
				}
			}
			ref_position = read.offset;
		}
		aux.full_read_match += (match.size() == 1 and mem_len(match[0]) == read.len);

		int old_n = (int)match.size();
		for (int j = 0; j < old_n; j++) {
			const auto &p = match[j] ;
			int beg = mem_beg(p), end = mem_end(p);
			if (end - beg < (int)(1.0 * opt->min_seed_len * opt->split_factor + .499) or p.x[2] > opt->split_width) continue;
			collect_mem_with_sst(bases, read.len, (beg + end) / 2, p.x[2] + 1, aux);
			for (const auto &m : aux.super_mem) {
				if (mem_len(m) >= opt->min_seed_len)
					match.push_back(m);
			}
		}

		if (opt->max_mem_intv > 0) {
			for (int j = 0; j < read.len; ) {
				if (bases[j] < 4) {
					bwtintv_t m;
					j = tem_forward_sst(bases, read.len, j, &m, aux);
					if (m.x[2] > 0) match.push_back(m);
				} else {
					j++;
				}
			}
		}
		std::sort(match.begin(), match.end(), mem_cmp);
	}

	// Find unique hit locations by merging and sorting
	auto unique_sal = aux.unique_sal; unique_sal.clear();
	for (int read_id = 0; read_id < n; read_id++) {
		const auto &mem = aux.match[read_id];
		auto &seed = aux.seed[read_id]; seed.clear();
		int array_id = 0;
		for (const auto &m : mem) {
			uint64_t step = m.x[2] > opt->max_occ ? m.x[2] / opt->max_occ : 1;
			for (uint64_t k = 0, count = 0; k < m.x[2] && count < opt->max_occ; k += step, count++) {
				seed_hit s;
				s.qbeg = mem_beg(m);
				s.score = s.len = mem_len(m);
				seed.push_back(s);
				unique_sal.emplace_back(sal_request(m.x[0] + k, read_id, array_id));
				array_id++;
			}
		}
	}
	std::sort(unique_sal.begin(), unique_sal.end());
	// Lookup suffix array for hit locations
	uint64_t coordinate = 0;
	for (int i = 0; i < unique_sal.size(); i++) {
		const auto &p = unique_sal[i];
		if (i == 0 or unique_sal[i-1].que_location != p.que_location) {
			coordinate = bwt_sa(bwa_idx->bwt, p.que_location);
		}
		auto &s = aux.seed[p.read_id][p.array_id];
		s.rbeg = coordinate;
		s.rid = bns_intv2rid(bwa_idx->bns, s.rbeg, s.rbeg + s.len);
	}

	// Print seeds
	for (int i = 0; i < n; i++) {
		const auto &seed = aux.seed[i];
		auto chains = chaining(seed);

		// Calculate repetition fraction of chains
		const auto &mem = aux.match[i];
		int beg = 0, end = 0, l_repetition = 0;
		for (const auto &m : mem) { // Matches are sorted by interval already
			if (m.x[2] <= opt->max_occ) continue; // Not a repetitive match
			int b = mem_beg(m), e = mem_end(m);
			if (b > end) { l_repetition += end - beg; beg = b; end = e; }
			else end = std::max(end, e);
		}
		l_repetition += end - beg;
		for (auto &c : chains) {
			c.frac_rep = (float)l_repetition / (float)reads[_start + i].len;
		}

		kstring_t *d = &debug_out[_start + i];
		print_chains_to(chains, d);
		for (auto &c : chains) c.destroy();
	}
}

void CompAligner::run(const char *fn) {
	// Prepare for thread auxiliary
	thr_aux = new thread_aux[threads_n];
	for (int i = 0; i < threads_n; i++) {
		auto &a = thr_aux[i];
		a.forward_sst = new SST(bwa_idx->bwt);
		a.backward_sst = new SST(bwa_idx->bwt);
	}

	fprintf(stderr, "Input test data with strand-corrected reads and overlapping information\n");
	std::ifstream in(fn); assert(in.is_open());
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	long off, global_offset = 0;
	long bytes = 0, round_n = 0;
	while (in >> buffer >> off) {
		global_offset += off;
		ngs_read read1;
		read1.bases = strdup(buffer);
		read1.len = strlen(buffer);
		read1.offset = global_offset;
		reads.push_back(read1);
		bytes += read1.len;
		total_reads_count++;
		if (bytes >= chunk_size * threads_n) {
			debug_out = (kstring_t*) calloc(reads.size(), sizeof(kstring_t));
			kt_for(
				threads_n,
				[](void *d, long i, int t) -> void {
					((CompAligner*)d)->seed_and_extend(i*BATCH_SIZE, (i+1)*BATCH_SIZE, t);
				},
				this,
				((int)reads.size() + BATCH_SIZE - 1) / BATCH_SIZE
			);
			for (auto &r : reads) { free(r.bases); free(r.sam); }
			fprintf(stderr, "Seed and Extend: %d reads processed\n", total_reads_count);
			for (int i = 0; i < reads.size(); i++) {
				fprintf(stdout, "%s\n", debug_out[i].s);
				free(debug_out[i].s);
			}
			free(debug_out);
			reads.clear(); bytes = 0; round_n++;
			if (round_n >= input_round_limit) break;
		}
	}
	display_profile();

	free(opt);
	bwa_idx_destroy(bwa_idx);
	for (int i = 0; i < threads_n; i++) {
		auto &a = thr_aux[i];
		delete a.forward_sst;
		delete a.backward_sst;
	}
	delete [] thr_aux;
	in.close();
	delete [] buffer;
}

void CompAligner::bwamem(const char *fn) {
	opt->n_threads = threads_n;

	fprintf(stderr, "Input test data with strand-corrected reads and overlapping information\n");
	std::ifstream in(fn); assert(in.is_open());
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	long off, bytes = 0, round_n = 0, n_processed = 0;
	std::vector<bseq1_t> bwa_seqs;
	while (in >> buffer >> off) {
		bseq1_t read1{0};
		read1.id = total_reads_count;
		read1.name = strdup(std::to_string(total_reads_count).c_str());
		read1.seq = strdup(buffer);
		read1.l_seq = (int) strlen(buffer);
		bwa_seqs.push_back(read1);
		total_reads_count++;
		bytes += read1.l_seq;
		if (bytes >= chunk_size * threads_n) {
			int n = bwa_seqs.size();
			mem_process_seqs(opt, bwt, bns, pac, n_processed, n, bwa_seqs.data(), nullptr);
			n_processed += n;
			for (int i = 0; i < n; i++) {
				const auto &s = bwa_seqs[i];
//				fprintf(stdout, "%s", s.sam);
				free(s.name); free(s.seq); free(s.sam);
			}

			fprintf(stderr, "BWA-MEM: %d reads processed\n", total_reads_count);
			bwa_seqs.clear(); bytes = 0; round_n++;
			if (round_n >= input_round_limit) break;
		}
	}

	free(opt);
	bwa_idx_destroy(bwa_idx);
	in.close();
	delete [] buffer;
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <reads>\n");
		return 1;
	}
	CompAligner worker;
	int c; bool input_test = false;
	std::string mode = "test";
	while ((c = getopt(argc, argv, "t:m:r:")) >= 0) {
		if (c == 't') {
			worker.threads_n = strtol(optarg, nullptr, 10);
		} else if (c == 'r') {
			worker.input_round_limit = strtol(optarg, nullptr, 10);
		} else if (c == 'm') {
			mode = std::string(optarg);
		} else {
			fprintf(stderr, "Unrecognized option\n");
			exit(EXIT_FAILURE);
		}
	}
	worker.load_index(argv[optind]);
	if (mode == "comp") {
		worker.run(argv[optind + 1]);
	} else if (mode == "bwa") {
		worker.bwamem(argv[optind + 1]);
	}
	return 0;
}
