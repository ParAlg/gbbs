#include <cstdlib>

template<class T>
struct SplayNode {
	SplayNode *child[2], *father;
	int size;
	T val;
	SplayNode(T val = 0): val(val), size(1), father(NULL) {
		child[0] = child[1] = NULL;
	}
};

template<class T>
struct SplayTree {
public:
	typedef SplayNode<T> Splay;
	Splay *root;
	Splay *begindummy;
	Splay *enddummy;
	SplayTree(): root(new Splay), begindummy(root), enddummy(new Splay) {
		insertAfter(enddummy, begindummy);
	}

	void insertAfter(Splay *x, Splay *y) { // Insert x after y
		// You should guarantee y != NULL
		splay(y, NULL); // Used to push down the labels along the path from the root to y!
		Splay *z = y->child[1];
		if (z == NULL) {
			y->child[1] = x;
			x->father = y;
			refresh(y);
		}
		else {
			pushDown(z);
			while (z->child[0] != NULL)
				z = z->child[0], pushDown(z);
			z->child[0] = x;
			x->father = z;
			while (z != NULL)
				refresh(z), z = z->father;
		}
		x->child[0] = x->child[1] = NULL;
		splay(x, NULL);
	}

	Splay *selectKth(int k) { // Return the k-th element (indexing from 0)
		Splay *tree = root, *last;
		while (tree != NULL) {
			pushDown(tree);
			int leftSize = (tree->child[0] != NULL) ? tree->child[0]->size : 0;
			last = tree;
			if (leftSize == k) {
				splay(tree, NULL);
				return tree;
			}
			else if (leftSize > k) tree = tree->child[0];
			else k -= leftSize + 1, tree = tree->child[1];
		}
		splay(last, NULL);
		return NULL; // K-th element does not exist (the tree has no greater than k elements)
	}

	Splay *neighbor(Splay *x, bool dir) {
		splay(x, NULL); // Used to push down the labels along the path from the root to x!
		if (x->child[dir] == NULL) return NULL;
		x = x->child[dir];
		pushDown(x);
		while (x->child[!dir] != NULL) x = x->child[!dir], pushDown(x);
		return x;
	}

	Splay *prev(Splay *x) {
		return neighbor(x, 0);
	}

	Splay *succ(Splay *x) {
		return neighbor(x, 1);
	}

	void del(Splay *x) { // Delete x from the tree
		splay(x, NULL);
		if (x->child[0] == NULL) {
			root = x->child[1];
			if (x->child[1] != NULL) x->child[1]->father = NULL;
		}
		else {
			Splay *y = prev(x);
			splay(y, x);
			y->child[1] = x->child[1];
			y->father = NULL;
			if (x->child[1] != NULL) x->child[1]->father = y;
			root = y;
			refresh(y);
		}
	}

	int rank(Splay *x) { // Return the ranking of x (indexing from 0)
		splay(x, NULL);
		if (x->child[0] == NULL) return 0;
		return x->child[0]->size;
	}

private:
	void refresh(Splay *x) {
		x->size = 1;
		if (x->child[0] != NULL) x->size += x->child[0]->size;
		if (x->child[1] != NULL) x->size += x->child[1]->size;
		// Refresh other information here
	}

	void pushDown(Splay *x) {
	}

	void rotate(Splay *x, bool dir) {
		// x != NULL, and x->father != NULL
	/*	y												 x
		 / \	 	rotate(x, 0)		 	/ \
		o   x   ------------->   y   o
			 / \  <-------------  / \
			o   o  rotate(y, 1)  o   o   */
		Splay *y = x->father;
		pushDown(y);
		pushDown(x);
		y->child[!dir]=x->child[dir];
		if (x->child[dir] != NULL) x->child[dir]->father = y;
		x->father = y->father;
		if (y->father != NULL)
			if (y->father->child[0] == y) y->father->child[0] = x;
			else y->father->child[1] = x;
		x->child[dir] = y;
		y->father = x;
		if (y == root) root = x;
		refresh(y);
		refresh(x);
	}

	void splay(Splay *x, Splay *f) {
		if (x != NULL) pushDown(x);
		if (x == f || x == NULL) return;
		while (x->father != f) {
			if (x->father->father == f) {
				pushDown(x->father);
				pushDown(x);
				rotate(x, x->father->child[0] == x);
			}
			else {
				Splay *y = x->father;
				Splay *z = y->father;
				pushDown(z);
				pushDown(y);
				pushDown(x);
				if (z->child[0] == y)
					if (y->child[0] == x)
						rotate(y, 1), rotate(x, 1);
					else
						rotate(x, 0), rotate(x, 1);
				else
					if (y->child[0] == x)
						rotate(x, 1), rotate(x, 0);
					else
						rotate(y, 0), rotate(x, 0);
			}
		}
		if (f == NULL) root = x;
	//	if (f != NULL) refresh(f); // Is it useful?
	}
};
