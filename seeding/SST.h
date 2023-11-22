//
// Created by ixiaohu on 2023/6/8.
//

#ifndef COMP_SEED_SST_H
#define COMP_SEED_SST_H

#include "../FM_index/bwt.h"
#include <vector>

/** SMEM Search Tree (SST) Node */
struct SST_Node_t {
	bwtintv_t match; // SA interval in FMD-index
	int children[4]; // Indexes of four children
	SST_Node_t() { children[0] = children[1] = children[2] = children[3] = -1; }
};

class SST {
private:
	std::vector<SST_Node_t> nodes;
	const bwt_t *bwt;
	bwtintv_t next[4];

public:
	explicit SST(const bwt_t *b);

	int64_t bwt_call = 0; /** Real call times of BWT-extend */

	/** At parent node with prefix p, query child node for p + base.
	 * If the child node does not exist, query it from FMD-index. */
	int query_child(int parent, uint8_t base, bool is_back);

	inline int query_forward_child(int parent, uint8_t base);

	inline int query_backward_child(int parent, uint8_t base);

	/** Add [pivot, LEP] to backward SST. */
	inline int add_lep_child(int parent, uint8_t base, const uint64_t *x);

	/** Add suffixes of [pivot + 1, LEP] to backward SST. But SA interval of
	 * suffixes is unknown so {0,0,0} is used to mark such nodes in SST. */
	inline int add_empty_child(int parent, uint8_t base);

	inline bwtintv_t get_intv(int id) { return nodes[id].match; }

	int get_child(int parent, uint8_t base) { return nodes[parent].children[base]; }

	/** Only keep root and its four children */
	inline void clear() {
		bwt_call = 0;
		nodes.resize(5);
		for (int i = 1; i <= 4; i++) {
			for (int &c : nodes[i].children) {
				c = -1;
			}
		}
	}
};

inline int SST::query_forward_child(int parent, uint8_t base) {
	if (nodes[parent].children[base] == -1) {
		bwt_extend(bwt, &nodes[parent].match, next, 0);
		bwt_call++;

		SST_Node_t child;
		child.match = next[base];
		nodes[parent].children[base] = nodes.size();
		nodes.push_back(child);
	}
	return nodes[parent].children[base];
}

inline int SST::query_backward_child(int parent, uint8_t base) {
	if (nodes[parent].children[base] == -1) {
		bwt_extend(bwt, &nodes[parent].match, next, 1);
		bwt_call++;
		SST_Node_t child;
		child.match = next[base];
		nodes[parent].children[base] = nodes.size();
		nodes.push_back(child);
	}
	auto &c = nodes[nodes[parent].children[base]];
	// If find an empty node, BWT query is required
	if (c.match.x[0] or c.match.x[1] or c.match.x[2]) {
		return nodes[parent].children[base];
	} else {
		bwt_extend(bwt, &nodes[parent].match, next, 1);
		bwt_call++;
		c.match = next[base];
		return nodes[parent].children[base];
	}
}

inline int SST::add_lep_child(int parent, uint8_t base, const uint64_t *x) {
	if (nodes[parent].children[base] == -1) {
		nodes[parent].children[base] = nodes.size();
		SST_Node_t child;
		child.match.x[0] = x[0];
		child.match.x[1] = x[1];
		child.match.x[2] = x[2];
		nodes.push_back(child);
	} else {
		auto &c = nodes[nodes[parent].children[base]];
		c.match.x[0] = x[0];
		c.match.x[1] = x[1];
		c.match.x[2] = x[2];
	}
	return nodes[parent].children[base];
}

inline int SST::add_empty_child(int parent, uint8_t base) {
	if (nodes[parent].children[base] == -1) {
		nodes[parent].children[base] = nodes.size();
		SST_Node_t child;
		child.match.x[0] = child.match.x[1] = child.match.x[2] = 0;
		nodes.push_back(child);
	} // else do nothing whatever the child node is empty or not
	return nodes[parent].children[base];
}

#endif //COMP_SEED_SST_H
