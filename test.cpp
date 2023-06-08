//
// Created by ixiaohu on 2023/6/2.
//

#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include <sys/resource.h>
#include <sys/time.h>
#include <algorithm>
#include "cstl/kbtree.h"
#include "PgRC_BWA_MEM/BTree.h"

std::map<const kbnode_t*, int> node_pointer;

int node_id(const kbnode_t *x) {
	if (node_pointer.find(x) == node_pointer.end()) {
		node_pointer[x] = node_pointer.size();
	}
	return node_pointer[x];
}

std::string node_str(const kbnode_t *x) {
	char buf[1024] = {0};
	char *p = buf;
	p += sprintf(p, "[%d]%c(", node_id(x), x->is_internal ?'I' :'E');
	for (int i = 0; i < x->n; i++) {
		if (i > 0) p += sprintf(p, " ");
		p += sprintf(p, "%d", __KB_KEY(int, x)[i]);
	}
	sprintf(p, ")");
	return std::string(buf);
}

static double localtime() {
	timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

struct Key {
	int value, index;
	Key() = default;
	Key(int v , int i): value(v), index(i) {}
};

#define cmp(a, b) ((a).value - (b).value)

KBTREE_INIT(qqq, Key, cmp)

int main() {
	std::vector<Key> array;
	int N = 500000000, MAX_VALUE = 5000;
	for (int i = 0; i < N; i++) {
		array.emplace_back(Key(rand() % MAX_VALUE, i));
	}
	double t_start, total_start;
	std::vector<Key> sorted;
	std::vector<Key> truth;

	total_start = t_start = localtime();
	kbtree_t(qqq) *t = kb_init(qqq, KB_DEFAULT_SIZE);
	for (auto a : array) kb_put(qqq, t, a);
	fprintf(stderr, "Add: %.2f\n", localtime() - t_start);

	t_start = localtime();
	#define traversal_func(a) (truth.push_back(*(a)))
	__kb_traverse(Key, t, traversal_func);
	#undef traversal_func
	fprintf(stderr, "traverse: %.2f\n", localtime() - t_start);

	t_start = localtime();
	__kb_destroy(t);
	fprintf(stderr, "destroy: %.2f\n", localtime() - t_start);
	fprintf(stderr, "C B-Tree: %.2f\n", localtime() - total_start);

	// B-Tree (stack)
	total_start = t_start = localtime();
	BTree<Key> tree(512, [](const Key &a, const Key &b) -> int { return a.value - b.value; });
	for (auto a : array) tree.add(a);
	fprintf(stderr, "Add: %.2f\n", localtime() - t_start);

	t_start = localtime();
	tree.traverse(sorted);
	fprintf(stderr, "traverse: %.2f\n", localtime() - t_start);

	t_start = localtime();
	tree.destroy();
	fprintf(stderr, "destroy: %.2f\n", localtime() - t_start);
	fprintf(stderr, "C++ B-Tree: %.2f\n", localtime() - total_start);

	assert(sorted.size() == N); assert(truth.size() == N);
	for (int i = 0; i < truth.size(); i++) {
		assert(truth[i].value == sorted[i].value);
		assert(truth[i].index == sorted[i].index);
	}
}

