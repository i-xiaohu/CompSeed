//
// Created by ixiaohu on 2023/6/2.
//

#include <cstdio>
#include <string>
#include <map>
#include "cstl/kbtree.h"
#include "PgRC_BWA_MEM/BTree.h"

#define cmp(a, b) ((a) - (b))

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

typedef struct {
	kbnode_t *root;
	int off_key, off_ptr, ilen, elen;
	int n, t;
	int n_keys, n_nodes;
} kbtree_test_t;
kbtree_test_t *kb_init_test(int size) {
	kbtree_test_t *b;
	b = (kbtree_test_t *) calloc(1, sizeof(kbtree_test_t));
	b->t = ((size - 4 - sizeof(void *)) / (sizeof(void *) + sizeof(int)) + 1) >> 1;
	if (b->t < 2) {
		free(b);
		return 0;
	}
	b->n = 2 * b->t - 1;
	b->off_ptr = 4 + b->n * sizeof(int);
	b->ilen = (4 + sizeof(void *) + b->n * (sizeof(void *) + sizeof(int)) + 3) >> 2 << 2;
	b->elen = (b->off_ptr + 3) >> 2 << 2;
	fprintf(stderr, "ilen = %d, elen = %d\n", b->ilen, b->elen);
	b->root = (kbnode_t *) calloc(1, b->ilen);
	++b->n_nodes;
	return b;
}
static inline int __kb_getp_aux_test(const kbnode_t *x, const int *k, int *r) {
	int tr, *rr, begin = 0, end = x->n;
	if (x->n == 0)return -1;
	rr = r ? r : &tr;
	while (begin < end) {
		int mid = (begin + end) >> 1;
		if (((((int *) ((char *) x + 4))[mid]) - (*k)) < 0)begin = mid + 1; else end = mid;
	}
	if (begin == x->n) {
		*rr = 1;
		return x->n - 1;
	}
	if ((*rr = ((*k) - (((int *) ((char *) x + 4))[begin]))) < 0)--begin;
	return begin;
}
static int *kb_getp_test(kbtree_test_t *b, const int *k) {
	int i, r = 0;
	kbnode_t *x = b->root;
	while (x) {
		i = __kb_getp_aux_test(x, k, &r);
		if (i >= 0 && r == 0)return &((int *) ((char *) x + 4))[i];
		if (x->is_internal == 0)return 0;
		x = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1];
	}
	return 0;
}
static inline int *kb_get_test(kbtree_test_t *b, const int k) { return kb_getp_test(b, &k); }
static void kb_intervalp_test(kbtree_test_t *b, const int *k, int **lower, int **upper) {
	int i, r = 0;
	kbnode_t *x = b->root;
	*lower = *upper = 0;
	while (x) {
		i = __kb_getp_aux_test(x, k, &r);
		if (i >= 0 && r == 0) {
			*lower = *upper = &((int *) ((char *) x + 4))[i];
			return;
		}
		if (i >= 0)*lower = &((int *) ((char *) x + 4))[i];
		if (i < x->n - 1)*upper = &((int *) ((char *) x + 4))[i + 1];
		if (x->is_internal == 0)return;
		x = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1];
	}
}
static inline void kb_interval_test(kbtree_test_t *b, const int k, int **lower, int **upper) {
	kb_intervalp_test(b, &k, lower, upper);
}
static void __kb_split_test(kbtree_test_t *b, kbnode_t *x, int i, kbnode_t *y) {
	kbnode_t *z;
	z = (kbnode_t *) calloc(1, y->is_internal ? b->ilen : b->elen);
	++b->n_nodes;
	z->is_internal = y->is_internal;
	z->n = b->t - 1;
	memcpy(((int *) ((char *) z + 4)), ((int *) ((char *) y + 4)) + b->t, sizeof(int) * (b->t - 1));
	if (y->is_internal)
		memcpy(((kbnode_t **) ((char *) z + b->off_ptr)), ((kbnode_t **) ((char *) y + b->off_ptr)) + b->t,
		       sizeof(void *) * b->t);
	y->n = b->t - 1;
	memmove(((kbnode_t **) ((char *) x + b->off_ptr)) + i + 2, ((kbnode_t **) ((char *) x + b->off_ptr)) + i + 1,
	        sizeof(void *) * (x->n - i));
	((kbnode_t **) ((char *) x + b->off_ptr))[i + 1] = z;
	memmove(((int *) ((char *) x + 4)) + i + 1, ((int *) ((char *) x + 4)) + i, sizeof(int) * (x->n - i));
	((int *) ((char *) x + 4))[i] = ((int *) ((char *) y + 4))[b->t - 1];
	++x->n;
	fprintf(stderr, "x = %s, y = %s, z = %s\n", node_str(x).c_str(), node_str(y).c_str(), node_str(z).c_str());
}
static void __kb_putp_aux_test(kbtree_test_t *b, kbnode_t *x, const int *k) {
	int i = x->n - 1;
	if (x->is_internal == 0) {
		i = __kb_getp_aux_test(x, k, 0);
		if (i != x->n - 1)
			memmove(((int *) ((char *) x + 4)) + i + 2, ((int *) ((char *) x + 4)) + i + 1,
			        (x->n - i - 1) * sizeof(int));
		((int *) ((char *) x + 4))[i + 1] = *k;
		++x->n; // Extern node is guaranteed to be not full
		fprintf(stderr, "Put key %d into external node %s\n", *k, node_str(x).c_str());
	}
	else {
		i = __kb_getp_aux_test(x, k, 0) + 1;
		fprintf(stderr, "Position of key %d in node %s is %d, which points to %d\n",
		  *k, node_str(x).c_str(), i, node_id(__KB_PTR(b, x)[i]));
		if (((kbnode_t **) ((char *) x + b->off_ptr))[i]->n == 2 * b->t - 1) {
			fprintf(stderr, "But node %s has been full\n", node_str(__KB_PTR(b, x)[i]).c_str());
			__kb_split_test(b, x, i, ((kbnode_t **) ((char *) x + b->off_ptr))[i]);
			if (((*k) - (((int *) ((char *) x + 4))[i])) > 0) ++i;
		}
		__kb_putp_aux_test(b, ((kbnode_t **) ((char *) x + b->off_ptr))[i], k);
	}
}
static void kb_putp_test(kbtree_test_t *b, const int *k) {
	kbnode_t *r, *s;
	++b->n_keys;
	r = b->root;
	if (r->n == 2 * b->t - 1) {
		fprintf(stderr, "Root node %s has been full, split it into y and z\n", node_str(r).c_str());
		++b->n_nodes;
		s = (kbnode_t *) calloc(1, b->ilen);
		b->root = s;
		s->is_internal = 1;
		s->n = 0;
		((kbnode_t **) ((char *) s + b->off_ptr))[0] = r;
		__kb_split_test(b, s, 0, r);
		r = s;
	}
	__kb_putp_aux_test(b, r, k);
}

static inline void kb_put_test(kbtree_test_t *b, const int k) { kb_putp_test(b, &k); }

static int __kb_delp_aux_test(kbtree_test_t *b, kbnode_t *x, const int *k, int s) {
	int yn, zn, i, r = 0;
	kbnode_t *xp, *y, *z;
	int kp;
	if (x == 0)return *k;
	if (s) {
		r = x->is_internal == 0 ? 0 : s == 1 ? 1 : -1;
		i = s == 1 ? x->n - 1 : -1;
	}
	else i = __kb_getp_aux_test(x, k, &r);
	if (x->is_internal == 0) {
		if (s == 2)++i;
		kp = ((int *) ((char *) x + 4))[i];
		memmove(((int *) ((char *) x + 4)) + i, ((int *) ((char *) x + 4)) + i + 1, (x->n - i - 1) * sizeof(int));
		--x->n;
		return kp;
	}
	if (r == 0) {
		if ((yn = ((kbnode_t **) ((char *) x + b->off_ptr))[i]->n) >= b->t) {
			xp = ((kbnode_t **) ((char *) x + b->off_ptr))[i];
			kp = ((int *) ((char *) x + 4))[i];
			((int *) ((char *) x + 4))[i] = __kb_delp_aux_test(b, xp, 0, 1);
			return kp;
		}
		else if ((zn = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1]->n) >= b->t) {
			xp = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1];
			kp = ((int *) ((char *) x + 4))[i];
			((int *) ((char *) x + 4))[i] = __kb_delp_aux_test(b, xp, 0, 2);
			return kp;
		}
		else if (yn == b->t - 1 && zn == b->t - 1) {
			y = ((kbnode_t **) ((char *) x + b->off_ptr))[i];
			z = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1];
			((int *) ((char *) y + 4))[y->n++] = *k;
			memmove(((int *) ((char *) y + 4)) + y->n, ((int *) ((char *) z + 4)), z->n * sizeof(int));
			if (y->is_internal)
				memmove(((kbnode_t **) ((char *) y + b->off_ptr)) + y->n, ((kbnode_t **) ((char *) z + b->off_ptr)),
				        (z->n + 1) * sizeof(void *));
			y->n += z->n;
			memmove(((int *) ((char *) x + 4)) + i, ((int *) ((char *) x + 4)) + i + 1, (x->n - i - 1) * sizeof(int));
			memmove(((kbnode_t **) ((char *) x + b->off_ptr)) + i + 1,
			        ((kbnode_t **) ((char *) x + b->off_ptr)) + i + 2, (x->n - i - 1) * sizeof(void *));
			--x->n;
			free(z);
			return __kb_delp_aux_test(b, y, k, s);
		}
	}
	++i;
	if ((xp = ((kbnode_t **) ((char *) x + b->off_ptr))[i])->n == b->t - 1) {
		if (i > 0 && (y = ((kbnode_t **) ((char *) x + b->off_ptr))[i - 1])->n >= b->t) {
			memmove(((int *) ((char *) xp + 4)) + 1, ((int *) ((char *) xp + 4)), xp->n * sizeof(int));
			if (xp->is_internal)
				memmove(((kbnode_t **) ((char *) xp + b->off_ptr)) + 1, ((kbnode_t **) ((char *) xp + b->off_ptr)),
				        (xp->n + 1) * sizeof(void *));
			((int *) ((char *) xp + 4))[0] = ((int *) ((char *) x + 4))[i - 1];
			((int *) ((char *) x + 4))[i - 1] = ((int *) ((char *) y + 4))[y->n - 1];
			if (xp->is_internal)
				((kbnode_t **) ((char *) xp + b->off_ptr))[0] = ((kbnode_t **) ((char *) y + b->off_ptr))[y->n];
			--y->n;
			++xp->n;
		}
		else if (i < x->n && (y = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1])->n >= b->t) {
			((int *) ((char *) xp + 4))[xp->n++] = ((int *) ((char *) x + 4))[i];
			((int *) ((char *) x + 4))[i] = ((int *) ((char *) y + 4))[0];
			if (xp->is_internal)
				((kbnode_t **) ((char *) xp + b->off_ptr))[xp->n] = ((kbnode_t **) ((char *) y + b->off_ptr))[0];
			--y->n;
			memmove(((int *) ((char *) y + 4)), ((int *) ((char *) y + 4)) + 1, y->n * sizeof(int));
			if (y->is_internal)
				memmove(((kbnode_t **) ((char *) y + b->off_ptr)), ((kbnode_t **) ((char *) y + b->off_ptr)) + 1,
				        (y->n + 1) * sizeof(void *));
		}
		else if (i > 0 && (y = ((kbnode_t **) ((char *) x + b->off_ptr))[i - 1])->n == b->t - 1) {
			((int *) ((char *) y + 4))[y->n++] = ((int *) ((char *) x + 4))[i - 1];
			memmove(((int *) ((char *) y + 4)) + y->n, ((int *) ((char *) xp + 4)), xp->n * sizeof(int));
			if (y->is_internal)
				memmove(((kbnode_t **) ((char *) y + b->off_ptr)) + y->n, ((kbnode_t **) ((char *) xp + b->off_ptr)),
				        (xp->n + 1) * sizeof(void *));
			y->n += xp->n;
			memmove(((int *) ((char *) x + 4)) + i - 1, ((int *) ((char *) x + 4)) + i, (x->n - i) * sizeof(int));
			memmove(((kbnode_t **) ((char *) x + b->off_ptr)) + i, ((kbnode_t **) ((char *) x + b->off_ptr)) + i + 1,
			        (x->n - i) * sizeof(void *));
			--x->n;
			free(xp);
			xp = y;
		}
		else if (i < x->n && (y = ((kbnode_t **) ((char *) x + b->off_ptr))[i + 1])->n == b->t - 1) {
			((int *) ((char *) xp + 4))[xp->n++] = ((int *) ((char *) x + 4))[i];
			memmove(((int *) ((char *) xp + 4)) + xp->n, ((int *) ((char *) y + 4)), y->n * sizeof(int));
			if (xp->is_internal)
				memmove(((kbnode_t **) ((char *) xp + b->off_ptr)) + xp->n, ((kbnode_t **) ((char *) y + b->off_ptr)),
				        (y->n + 1) * sizeof(void *));
			xp->n += y->n;
			memmove(((int *) ((char *) x + 4)) + i, ((int *) ((char *) x + 4)) + i + 1, (x->n - i - 1) * sizeof(int));
			memmove(((kbnode_t **) ((char *) x + b->off_ptr)) + i + 1,
			        ((kbnode_t **) ((char *) x + b->off_ptr)) + i + 2, (x->n - i - 1) * sizeof(void *));
			--x->n;
			free(y);
		}
	}
	return __kb_delp_aux_test(b, xp, k, s);
}
static int kb_delp_test(kbtree_test_t *b, const int *k) {
	kbnode_t *x;
	int ret;
	ret = __kb_delp_aux_test(b, b->root, k, 0);
	--b->n_keys;
	if (b->root->n == 0 && b->root->is_internal) {
		--b->n_nodes;
		x = b->root;
		b->root = ((kbnode_t **) ((char *) x + b->off_ptr))[0];
		free(x);
	}
	return ret;
}
static inline int kb_del_test(kbtree_test_t *b, const int k) { return kb_delp_test(b, &k); }

static void print_btree(const kbtree_test_t *b, const kbnode_t *r) {
	if (r->is_internal == 0) { // External node; print keys only
		fprintf(stdout, "external %d %d", node_id(r), r->n);
		for (int i = 0; i < r->n; i++) {
			fprintf(stdout, " %d", __KB_KEY(int, r)[i]);
		}
		fprintf(stdout, "\n");
	} else {
		fprintf(stdout, "internal %d %d", node_id(r), r->n);
		for (int i = 0; i < r->n; i++) {
			fprintf(stdout, " %d", __KB_KEY(int, r)[i]);
		}
		for (int i = 0; i <= r->n; i++) {
			fprintf(stdout, " %d", node_id(__KB_PTR(b, r)[i]));
		}
		fprintf(stdout, "\n");

		for (int i = 0; i <= r->n; i++) {
			print_btree(b, __KB_PTR(b, r)[i]);
		}
	}
}

const int BLOCK = 1024;

int main() {
	kbtree_test_t *t = kb_init_test(64);
	fprintf(stderr, "Internal nodes range: [%d,%d]\n", t->t, t->n);
	// An normal random array
//	int array[] = {1, 3, 7, 2, 18, 2, 65, 4, 9, 11, 13};
	// An extreme case: all same value
	int array[] = {9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9};
	for (auto key : array) fprintf(stderr, "%d ", key); fprintf(stderr, "\n");
	for (auto key : array) {
		kb_put_test(t, key);
	}
	for (int i = 0; i < 50; i++) kb_put_test(t, 9);
	print_btree(t, t->root);

}

