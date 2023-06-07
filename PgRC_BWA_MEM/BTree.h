//
// Created by ixiaohu on 2023/6/2.
//

#ifndef PGRC_LEARN_BTREE_H
#define PGRC_LEARN_BTREE_H

#include <stdint.h>
#include <stdlib.h>
#include <cassert>
#include <stack>
#include <map>

template <class T>
class BTree {
private:
	/** A node in B-Tree can be internal or external. An internal node stores n keys and n + 1
	 * pointers to children with the rule that keys in child[i] < key[i] < keys in child[i+1].
	 * External nodes have only n keys. */
	struct node_t { int32_t is_internal: 1, n: 31; };

	node_t *root = nullptr;
	int off_ptr = 0; /** offset for pointers to children: sizeof(node) + (2t-1) * sizeof(T) */
	int t = 0; /** Order of B-tree t; number of keys in a node is in range [t-1, 2t-1] */
	int m = 0; /** m = 2t-1, the maximum keys in node */
	int in_mem = 0; /** Required memory of internal nodes for storing keys and pointers */
	int ex_mem = 0; /** Required memory of external nodes for storing keys */
	int keys_n = 0; /** Total number of keys or size of the tree */

	const int DEFAULT_BLOCK_SIZE = 512;

	/** Initialize a tree with the given block size for storing one node */
	inline void init(int size) {
		t = ((size - sizeof(node_t) - sizeof(node_t*)) / (sizeof(node_t*) + sizeof(T)) + 1) / 2;
		assert(t >= 2); // At least be order of 2
		m = 2 * t - 1;
		in_mem = (sizeof(node_t) + m * sizeof(T) + (m + 1) * sizeof(node_t*) + 3) / 4 * 4;
		ex_mem = (sizeof(node_t) + m * sizeof(T) + 3) / 4 * 4;
		off_ptr = sizeof(node_t) + m * sizeof(T);
		root = (node_t*) calloc(1, in_mem);
	}

	int(*cmp)(const T&, const T&);

	/** Return pointer to children[i] */
	inline node_t** PTR(const node_t* s) { return ((node_t**)((char*)s + off_ptr)); }

	// Return pointer for elements
	inline T* KEY(const node_t *s) { return ((T*)((char*)s + 4)); }

	/** For debug and print */
	std::map<const node_t*, int> pointer_dict;

	inline int node_id(const node_t *x) {
		if (pointer_dict.find(x) == pointer_dict.end()) {
			pointer_dict[x] = pointer_dict.size();
		}
		return pointer_dict[x];
	}

	inline std::string node_str(const node_t *x) {
		char buf[1024]; char *p = buf;
		p += sprintf(p, "n=%d, p=%d {", x->n, node_id(x));
		if (x->is_internal) {
			for (int i = 0; i < x->n; i++) {
				p += sprintf(p, "%d < %d <", node_id(PTR(x)[i]), KEY(x)[i]);
			}
			p += sprintf(p, " %d", node_id(PTR(x)[x->n]));
		} else {
			for (int i = 0; i < x->n; i++) {
				p += sprintf(p, "%s%d", i ? " " : "", KEY(x)[i]);
			}
		}
		sprintf(p, "}");
		return std::string(buf);
	}

	/** Split the node y which has the maximum 2t-1 keys.
	 * Move the second half y[t, 2t-1) to a new node z, and y keeps the first half [0, t-1) itself.
	 * Put y[t-1] to the vacated position x[i+1], the both sides of which point to y and z, respectively.
	 * Here the B-tree rule is hold: y < x[i+1] < z.
	 * Note: x is the parent of y, and it is guaranteed to be not full. */
	inline void split(node_t *x, int i, node_t *y) {
		node_t *z = (node_t*) calloc(1, y->is_internal ? in_mem : ex_mem);
		z->is_internal = y->is_internal; z->n = t - 1;
		memcpy(KEY(z), KEY(y) + t, sizeof(T) * (t - 1));
		if (y->is_internal) memcpy(PTR(z), PTR(y) + t, sizeof(T*) * t);

		y->n = t - 1;
		PTR(x)[i] = y; // Pointer to the small half

		memmove(PTR(x) + i + 2, PTR(x) + i + 1, sizeof(T*) * (x->n - i));
		PTR(x)[i + 1] = z; // Pointer to the big half
		memmove(KEY(x) + i + 1, KEY(x) + i, sizeof(T) * (x->n - i));
		KEY(x)[i] = KEY(y)[t-1];
		x->n++;
	}

	/** Binary search k in node x. Returned index is the last key < k or the first key = k in node x. */
	inline int binary_search_node(const node_t *x, const T &k) {
		int begin = 0, end = x->n;
		while (begin < end) {
			int mid = (begin + end) / 2;
			if (cmp(KEY(x)[mid], k) < 0) begin = mid + 1;
			else end = mid;
		}
		if (begin == x->n) return x->n-1; // No value >= k
		if (cmp(KEY(x)[begin], k) > 0) return begin-1; // No value = k
		return begin; // Value = k (maybe it is okay to return begin-1, too)
	}

	/** In practice, DFS is faster than stack-implemented */
	inline void traverse_dfs(const node_t *x, std::vector<T> &ans) {
		if (x->is_internal) {
			for (int i = 0; i < x->n; i++) {
				traverse_dfs(PTR(x)[i], ans);
				ans.push_back(KEY(x)[i]);
			}
			traverse_dfs(PTR(x)[x->n], ans);
		} else {
			for (int i = 0; i < x->n; i++) {
				ans.push_back(KEY(x)[i]);
			}
		}
	}

public:
	BTree(int size, int(*cmp)(const T&, const T&)): cmp(cmp) { init(size); }
	explicit BTree(int(*cmp)(const T&, const T&)): cmp(cmp) { init(DEFAULT_BLOCK_SIZE); }

	inline void add(const T &k) {
		keys_n++;
		if (root->n == m) { // Special judge for overflow in root node
			node_t *s = (node_t*) calloc(1, in_mem); // New root
			s->is_internal = 1; s->n = 0;
			split(s, 0, root);
			root = s;
		}

		node_t *x = root;
		while (x->is_internal != 0) { // Recursively locate the external node to put k
			// The binary search result r is the last value < k or the first value = k
			// i = r + 1 means the equal value is always put to the right child of key[r]
			int i = binary_search_node(x, k) + 1;
			if (PTR(x)[i]->n == m) { // If the child has been full, split it
				split(x, i, PTR(x)[i]);
				if (cmp(k, KEY(x)[i]) > 0) i++; // k is in the right child of x[i]
			}
			x = PTR(x)[i];
		}
		int i = binary_search_node(x, k);
		memmove(KEY(x) + i + 2, KEY(x) + i + 1, (x->n - i - 1) * sizeof(T));
		KEY(x)[i+1] = k; x->n++; // Extern node is guaranteed to be not full
	}

	inline void print() {
		std::stack<node_t*> st; st.push(root);
		while (!st.empty()) {
			node_t *x = st.top(); st.pop();
			if (x->is_internal !=  0) {
				fprintf(stdout, "internal %d %d", node_id(x), x->n);
				for (int i = 0; i < x->n; i++) {
					fprintf(stdout, " %d", KEY(x)[i]);
				}
				for (int i = 0; i <= x->n; i++) {
					fprintf(stdout, " %d", node_id(PTR(x)[i]));
				}
				fprintf(stdout, "\n");

				for (int i = 0; i <= x->n; i++) {
					st.push(PTR(x)[i]);
				}
			} else {
				fprintf(stdout, "external %d %d", node_id(x), x->n);
				for (int i = 0; i < x->n; i++) {
					fprintf(stdout, " %d", KEY(x)[i]);
				}
				fprintf(stdout, "\n");
			}
		}
	}

	inline void traverse_stack(std::vector<T> &ans) {
		ans.clear();
		struct Item {
			const node_t *x; // The visited node
			int child; // which child of x is being visited
			Item(const node_t *x, int child): x(x), child(child) {}
		};
		std::stack<Item> st; st.push(Item(root, 0));
		while (!st.empty()) {
			auto item = st.top(); // Do not remove this node, which might have children
			// Depth First Search along this child
			while (item.x->is_internal != 0 and item.child <= item.x->n) {
				const node_t *next = PTR(item.x)[item.child];
				st.push(Item(next, 0));
				item = st.top();
			}
			if (item.x->is_internal == 0) { // reach a leaf node
				for (int i = 0; i < item.x->n; i++) {
					ans.push_back(KEY(item.x)[i]);
				}
			} // else all children of the node has been visited

			st.pop(); // Remove the processed node
			if (!st.empty()) { // Move up to the parent node
				auto parent = st.top(); st.pop();
				if (parent.child < parent.x->n) {
					ans.push_back(KEY(parent.x)[parent.child]);
				}
				parent.child++; // Move to next child
				st.push(parent);
			}
		}
	}

	inline void traverse(std::vector<T> &ans) {
		ans.clear();
		traverse_dfs(root, ans);
	}
};

#endif //PGRC_LEARN_BTREE_H
