//
// Created by ixiaohu on 2023/6/22.
//

#include "SST.h"

SST::SST(const bwt_t *b) {
	bwt = b;
	SST_Node_t root;
	for (int i = 0; i < 4; i++) root.children[i] = i + 1;
	nodes.push_back(root);

	SST_Node_t child;
	for (uint8_t c = 0; c < 4; c++) {
		bwt_set_intv(bwt, c, child.match);
		nodes.push_back(child);
	}
}

int SST::query_child(int parent, uint8_t base, bool is_back) {
	if (nodes[parent].children[base] == -1) {
		bwt_extend(bwt, &nodes[parent].match, next, is_back);
		bwt_calls++;

		for (uint8_t c = 0; c < 4; c++) { // Is it necessary?
			SST_Node_t child;
			child.match = next[c];
			nodes[parent].children[c] = nodes.size();
			nodes.push_back(child);
		}
	}
	auto &c = nodes[nodes[parent].children[base]];
	// If find an empty node, BWT query is required
	if (c.match.x[0] + c.match.x[1] + c.match.x[2] == 0) {
		bwt_extend(bwt, &nodes[parent].match, next, is_back);
		bwt_calls++;
		c.match = next[base];
	}
	return nodes[parent].children[base];
}