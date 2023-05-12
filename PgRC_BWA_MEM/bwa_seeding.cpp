//
// Created by ixiaohu on 2023/5/6.
//

#include "bwa_seeding.h"

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>

#include "../Matching/SimplePgMatcher.h"
#include "../Utils/helper.h"
#include "../FM_index/bntseq.h"
#include "../bwalib/bwa.h"
#include "../cstl/kvec.h"
#include "../cstl/kthread.h"

void BWA_seeding::load_all_PGs(std::ifstream &in) {
	fprintf(stderr, "1. Decompressing HQ reads...\n");
	PseudoGenomeHeader hq_pgh(in);
	auto *hq_prop = new ReadsSetProperties(in);
	assert(!confirm_text_read_mode(in));
	hq_reads_list = load_extend_reads_list(
			in,
			hq_pgh.get_max_read_length(),
			hq_pgh.get_reads_count(),
			preserve_order_mode,
			false,
			false);
	hq_reads_list->reads_count = hq_prop->reads_count;
	hq_reads_list->read_length = hq_prop->max_read_length;

	fprintf(stderr, "2. Decompressing LQ reads...\n");
	PseudoGenomeHeader lq_pgh(in);
	auto *lq_prop = new ReadsSetProperties(in);
	assert(!confirm_text_read_mode(in));
	lq_reads_list = load_extend_reads_list(
			in,
			lq_pgh.get_max_read_length(),
			lq_pgh.get_reads_count(),
			preserve_order_mode,
			true,
			true);
	lq_reads_list->reads_count = lq_prop->reads_count;
	lq_reads_list->read_length = lq_prop->max_read_length;

	fprintf(stderr, "3. Decompressing N reads...\n");
	PseudoGenomeHeader n_pgh;
	ReadsSetProperties *n_prop = nullptr;
	if (separate_N) {
		n_pgh = PseudoGenomeHeader(in);
		n_prop = new ReadsSetProperties(in);
		assert(!confirm_text_read_mode(in));
		n_reads_list = load_extend_reads_list(
				in,
				n_pgh.get_max_read_length(),
				n_pgh.get_reads_count(),
				preserve_order_mode,
				true,
				true);
		n_reads_list->reads_count = n_prop->reads_count;
		n_reads_list->read_length = n_prop->max_read_length;
	}

	read_length = hq_prop->max_read_length;
	hq_reads_count = hq_prop->reads_count;
	lq_reads_count = lq_prop->reads_count;
	non_reads_count = hq_reads_count + lq_reads_count;
	n_reads_count = n_prop ?n_prop->reads_count :0;
	total_reads_count = non_reads_count + n_reads_count;

	hq_pg_length = hq_pgh.get_pg_length();
	non_pg_length = hq_pg_length + lq_pgh.get_pg_length();
	if (preserve_order_mode) {
		joined_pg_len_std = non_pg_length + n_pgh.get_pg_length() <= UINT32_MAX;
		if (joined_pg_len_std) {
			decompress_pg_position(in, pos_pg_32, total_reads_count, single_end_mode);
		} else {
			decompress_pg_position(in, pos_pg_64, total_reads_count, single_end_mode);
		}
	} else if (not single_end_mode) {
		restore_paired_idx(in, paired_idx);
	}

	restore_all_matched_pg(in, hq_pg_length, hq_pg, lq_pg, n_pg);
}

template<typename uint_pg_len>
void BWA_seeding::apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos) {
	if (preserve_order_mode) {
		int hq_idx = 0;
		const int pairs_count = total_reads_count / 2;
		for (int i = 0; i < pairs_count; i++) {
			uint_pg_len pos = pg_pos[i];
			if (pos < hq_pg_length) hq_idx++;
		}
		for (int i = pairs_count; i < total_reads_count; i++) {
			uint_pg_len pos = pg_pos[i];
			if (pos < hq_pg_length) {
				hq_reads_list->rev_comp[hq_idx] = !hq_reads_list->rev_comp[hq_idx];
				hq_idx++;
			}
		}
	} else {
		// Reverse complement each read from file2
		for (int i = 1; i < total_reads_count; i += 2) {
			int idx = paired_idx[i];
			if (idx < hq_reads_count) {
				hq_reads_list->rev_comp[idx] = !hq_reads_list->rev_comp[idx];
			}
		}
	}
}

const int SEED_MEM = 1;
const int SEED_REM = 2;
const int SEED_TEM = 3;

struct re_seed_item {
	int position;
	long min_intv;
	re_seed_item(int p, long m): position(p), min_intv(m) {}
};
std::vector<re_seed_item> rs_items;
int global_offset;

void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a) {
	int i, k, x = 0, old_n;
	int start_width = 1;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	a->mem.n = 0;
	// first pass: find all SMEMs
	while (x < len) {
		if (seq[x] < 4) {
			x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i) {
				bwtintv_t *p = &a->mem1.a[i];
				p->type = SEED_MEM;
				int slen = (uint32_t)p->info - (p->info>>32); // seed length
				if (slen >= opt->min_seed_len)
					kv_push(bwtintv_t, a->mem, *p);
			}
		} else ++x;
	}
	// second pass: find MEMs inside a long SMEM
	old_n = a->mem.n;
	for (k = 0; k < old_n; ++k) {
		bwtintv_t *p = &a->mem.a[k];
		int start = p->info>>32, end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
		rs_items.emplace_back(re_seed_item((start + end) / 2 + global_offset, p->x[2] + 1));
		for (i = 0; i < a->mem1.n; ++i)
			if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len) {
				a->mem1.a[i].type = SEED_REM;
				kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
			}
	}
	// third pass: LAST-like
	if (opt->max_mem_intv > 0) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				if (1) {
					bwtintv_t m;
					m.type = SEED_TEM;
					x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
					if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
				} else { // for now, we never come to this block which is slower
					x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
					for (i = 0; i < a->mem1.n; ++i)
						kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
				}
			} else ++x;
		}
	}
}

struct exact_match_t {
	int b, e;
	int64_t k, l, s;
	exact_match_t(int b, int e, int64_t k, int64_t l, int64_t s):
		b(b), e(e), k(k), l(l), s(s) {}
};

void BWA_seeding::traditional_seeding(const std::string &cs, std::vector<long> &offset, std::vector<std::string> &batch) {
	// At the very beginning, convert bases into number
	for (auto &read: batch) {
		for (int i = 0; i < read.length(); i++) {
			read[i] = nst_nt4_table[(int) read[i]];
		}
	}

	// Construct graph
	long os_begin = offset.front(), os_end = offset.back() + read_length;
	long cs_length = os_end - os_begin;
	for (auto &o : offset) o -= os_begin; // For convenient access this batch
	int *coverage = (int*) calloc(cs_length, sizeof(int));
	uint64_t bit_sentry[4][cs_length]; // Graph representation
	memset(bit_sentry, 0, sizeof(bit_sentry));
	std::vector<int8_t> id_of_read(batch.size(), -1); // bit setter of each read

	uint n_id = 0;
	for (int i = 0; i < batch.size(); i++) {
		const auto &read = batch[i];
		bool overflow = false, with_n = false;
		for (int j = 0; j < read.length(); j++) {
			with_n |= (read[j] == 4);
			if (coverage[offset[i] + j] + 1 > 64) {
				overflow = true;
				break;
			}
		}
		read_with_n += with_n;
		overflow_read += overflow;
		if (with_n or overflow) continue; // Skip this read temporarily
		for (int j = 0; j < read.length(); j++) coverage[offset[i] + j]++;
		n_id = (n_id + 1) % 64;
		id_of_read[i] = n_id;
		for (int j = 0; j < read.length(); j++) {
			bit_sentry[read[j]][offset[i] + j] |= (1UL << n_id);
		}
	}
	free(coverage);

	total_cs_length += cs_length;
	for (int i = 1; i < cs_length; i++) {
		int cnt = 0;
		for (int c = 0; c < 4; c++) cnt += (bit_sentry[c][i] > 0);
		if (cnt > 1) divergent_column++;
	}


	// The original BWA-MEM seeding procedure
	rs_items.clear();
	std::vector<exact_match_t> matches;
	std::vector<exact_match_t> tems;
	for (int j = 0; j < batch.size(); j++) {
		const auto &read = batch[j];
		global_offset = offset[j];
		mem_collect_intv(mem_opt, bwa_idx->bwt, read.length(), (const uint8_t*) read.c_str(), mem_aux);
		for (int i = 0; i < mem_aux->mem.n; i++) {
			const auto &p = mem_aux->mem.a[i];
			int begin = p.info >> 32, end = (int)p.info;
			matches.emplace_back(exact_match_t(begin, end, p.x[0], p.x[1], p.x[2]));
			if (p.type == SEED_TEM) {
				tems.emplace_back(exact_match_t(begin + offset[j], end + offset[j], p.x[0], p.x[1], p.x[2]));
			}
		}
	}


	// Look insight into re-seeding
	reseed_count += rs_items.size();
	std::sort(rs_items.begin(), rs_items.end(),
		   [](const re_seed_item &l, const re_seed_item &r) -> bool
		   { return l.position < r.position; });
	for (int i = 1; i < rs_items.size(); i++) {
		if ((rs_items[i].position - rs_items[i-1].position) < 3 and
		rs_items[i].min_intv == rs_items[i-1].min_intv) {
			identical_reseed++;
		}
	}
	for (const auto &item : rs_items) {
		reseed_bin[item.min_intv]++;
	}

	// Dig the 3rd round seeding redundancy
	std::sort(tems.begin(), tems.end(),
		   [](const exact_match_t &l, const exact_match_t &r) -> bool
		   { return (l.b != r.b) ?(l.b < r.b) :(l.e < r.e); });
	tem_count += tems.size();
	for (int i = 1; i < tems.size(); i++) {
		if (tems[i].b == tems[i-1].b and tems[i].e == tems[i-1].e) {
			same_tem++;
		}
		if ((tems[i].l == tems[i-1].l and tems[i].s == tems[i-1].s) or
			(tems[i].k == tems[i-1].k and tems[i].s == tems[i-1].s)) {
			sal_tem++;
		}
	}


	// Exploit SAL redundancy
	std::sort(matches.begin(), matches.end(),
		   [] (const exact_match_t &l, const exact_match_t &r) -> bool
		   { return l.k < r.k; });
	for (const auto &m : matches) {
		int occ = m.s > mem_opt->max_occ ?mem_opt->max_occ :m.s;
		original_sal += occ;
	}
	for (int i = 0; i < matches.size(); i++) {
		const auto &curr = matches[i];
		int occ = curr.s > mem_opt->max_occ ?mem_opt->max_occ :curr.s;
		if (i == 0) {
			compressed_sal += occ;
			continue;
		}
		// If current MEM has same SA interval with previous one, skip it.
		const auto &prev = matches[i-1];
		if (prev.k == curr.k and prev.s == curr.s) {
			continue;
		}
		compressed_sal += occ;
	}
}

void BWA_seeding::seeding_SE() {
	char *buffer = new char[read_length * 2]; memset(buffer, 0, read_length * 2);
	long curr_position = 0, curr_mis_cnt = 0;
	std::vector<std::string> read_batch;
	std::vector<long> offset;
	for (int i = 0; i < hq_reads_list->reads_count; i++) {
		curr_position += hq_reads_list->off[i];
		memcpy(buffer, (hq_pg.data() + curr_position), hq_reads_list->read_length);
		if (hq_reads_list->rev_comp[i]) {
			rev_comp_in_place(buffer, hq_reads_list->read_length);
		}
		for (uint8_t j = 0; j < hq_reads_list->mis_cnt[i]; j++) {
			uint8_t mis_pos = hq_reads_list->mis_off[curr_mis_cnt];
			uint8_t mis_code = hq_reads_list->mis_sym_code[curr_mis_cnt];
			buffer[mis_pos] = code_to_mismatch(buffer[mis_pos], mis_code);
			curr_mis_cnt++;
		}

		// Reverse-complement it temporarily
		if (hq_reads_list->rev_comp[i]) {
			rev_comp_in_place(buffer, hq_reads_list->read_length);
		}
		read_batch.emplace_back(std::string(buffer));
		offset.push_back(curr_position);
		if (read_batch.size() >= COMP_BATCH_SIZE or i == hq_reads_list->reads_count - 1) {
			traditional_seeding(hq_pg, offset, read_batch);
			read_batch.clear();
			offset.clear();
		}
	}

	fprintf(stderr, "Graph construction: columns %ld, divergent %ld, space %.2f\n",
		 total_cs_length, divergent_column, 1.0 * total_cs_length / divergent_column);

	fprintf(stderr, "Reads with N %.2f, overflow read %.2f\n",
		 100.0 * read_with_n / hq_reads_count,
		 100.0 * overflow_read / hq_reads_count);

	fprintf(stderr, "Re-seeding redundancy: original %ld, identical %ld, percentage %.2f\n",
		 reseed_count, identical_reseed, 100.0 * identical_reseed / reseed_count);
	for (int i = 0; i < 500; i++) {
		if (reseed_bin[i] == 0) continue;
		fprintf(stderr, "  %.2f\n", 100.0 * reseed_bin[i] / reseed_count);
	}

	fprintf(stderr, "Third round seeding redundancy: %.2f, with same SAL: %.2f\n",
		 100.0 * same_tem / tem_count, 100.0 * sal_tem / tem_count);


	fprintf(stderr, "SAL redundancy: original %ld, compressed SAL %ld, percentage %.2f\n",
		 original_sal, compressed_sal, 100.0 * (original_sal - compressed_sal) / original_sal);


	curr_position = 0;
	auto *coverage = new int[hq_pg.length()];
	memset(coverage, 0, hq_pg.length() * sizeof(int));
	assert(hq_reads_list->off[0] == 0);
	for (int i = 0; i < hq_reads_list->reads_count; i++) {
		curr_position += hq_reads_list->off[i];
		for (int j = 0; j < hq_reads_list->read_length; j++) {
			coverage[curr_position + j]++;
		}
	}
	int fraction = 0;
	for (int i = 0; i < hq_pg.length(); i++) {
		if (coverage[i] > 64) fraction++;
	}
	fprintf(stderr, "Coverage > 64: %.2f\n", 100.0 * fraction / hq_pg.length());
	delete [] coverage;


	curr_position = 0;
	for (int i = 0; i < lq_reads_list->reads_count; i++) {
		curr_position += lq_reads_list->off[i];
		memcpy(buffer, (lq_pg.data() + curr_position), lq_reads_list->read_length);
	}

	curr_position = 0;
	for (int i = 0; i < n_reads_list->reads_count; i++) {
		curr_position += n_reads_list->off[i];
		memcpy(buffer, (n_pg.data() + curr_position), n_reads_list->read_length);
	}

	delete [] buffer;
}

mem_opt_t *mem_opt_init() {
	auto *o = (mem_opt_t*) calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

smem_aux_t *smem_aux_init() {
	smem_aux_t *a;
	a = (smem_aux_t*) calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	return a;
}

void smem_aux_destroy(smem_aux_t *a) {
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

void BWA_seeding::compressive_seeding() {
	// Loading FM-index
	bwa_idx = bwa_idx_load_from_shm(index_name.c_str());
	if (bwa_idx == nullptr) {
		bwa_idx = bwa_idx_load(index_name.c_str(), BWA_IDX_ALL);
		if (bwa_idx == nullptr) {
			fprintf(stderr, "Load the FM-index failed\n");
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "Load the FM-index from disk\n");
		}
	} else {
		fprintf(stderr, "Load the FM-index from shared memory\n");
	}
	mem_opt = mem_opt_init();
	mem_aux = smem_aux_init();


	std::ifstream in(archive_name); assert(in.is_open());
	for (int i = 0; i < strlen(PGRC_HEADER); i++) assert(in.get() == PGRC_HEADER[i]);
	char version_mode = in.get(); assert(version_mode == '#');
	char version_major = in.get();
	char version_minor = in.get();
	char version_revision = in.get();
	if (version_major != PGRC_VERSION_MAJOR or version_minor != PGRC_VERSION_MINOR) {
		fprintf(stderr, "Archive is packed with a different version PgRC %d.%d\n",
		        version_major, version_minor);
		exit(EXIT_FAILURE);
	}
	compression_level = in.get();
	char pgrc_mode = in.get();
	if (pgrc_mode != PGRC_SE_MODE) {
		fprintf(stderr, "PgRC-BWA-MEM now supports SE mode only\n");
		exit(EXIT_FAILURE);
	}
	const std::string five_modes[] = {"SE", "PE", "SE_ORD", "PE_ORD"};
	fprintf(stderr, "Compression mode: %s\n", five_modes[pgrc_mode].c_str());
	separate_N = (bool) in.get();
	if (pgrc_mode == PGRC_PE_MODE or pgrc_mode == PGRC_ORD_PE_MODE) {
		rev_comp_pair = (bool) in.get();
	}
	std::string tmp_directory; in >> tmp_directory; tmp_directory += "/";
	in.get();

	preserve_order_mode = (pgrc_mode == PGRC_ORD_SE_MODE or pgrc_mode == PGRC_ORD_PE_MODE);
	single_end_mode = (pgrc_mode == PGRC_SE_MODE or pgrc_mode == PGRC_ORD_SE_MODE);

	load_all_PGs(in);

	if (rev_comp_pair) {
		fprintf(stderr, "Reverse complement second read from a pair\n");
		if (joined_pg_len_std) apply_rc_pair_to_pg(pos_pg_32);
		else apply_rc_pair_to_pg(pos_pg_64);
	}

	if (pgrc_mode == PGRC_SE_MODE) seeding_SE();

	fprintf(stderr, "Decompressed %d reads in total\n", total_reads_count);
	in.close();
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <pgrc_archieve>\n");
		return 1;
	}
	BWA_seeding worker;
	worker.set_index_name(argv[1]);
	worker.set_archive_name(argv[2]);
	worker.compressive_seeding();
}
