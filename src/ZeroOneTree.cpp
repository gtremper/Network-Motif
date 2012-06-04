#include "ZeroOneTree.h"

#define MAX 1000

tree::tree(int subgraphSize) {
	this->subgraphSize = subgraphSize;
	root = new Node();
	leaf_num = 0;
	allocate_node();
	allocate_leaf();
}

tree::~tree() {
/*    delete [] node;
    delete [] leaf;
    delete root;*/
}

Node * tree::add_node() {
	Node * n = &node[node_ptr];
	node_ptr++;
	if(node_ptr == MAX)
		allocate_node();
	return n;
}

Node * tree::return_root() {
	return root;
}

Leaf * tree::add_leaf() {
	Leaf * l = &leaf[leaf_ptr];
	leaf_ptr++;
	leaf_num++;
	if(leaf_ptr == MAX)
		allocate_leaf();
	return l;
}


Leaf * tree::update_leaf(Leaf * l, int val) {
	l->count += val;
	return l;
}

void tree::allocate_node() {
	node = new Node [MAX];
	node_ptr = 0;
}

void tree::allocate_leaf() {
	leaf = new Leaf [MAX];
	leaf_ptr = 0;
}

int tree::get_leafnum() {
	return leaf_num;
}

void tree::init_cur_node() {
	cur_node = return_root();
}

void tree::insert_one_main() {
	if(cur_node->right)
		cur_node = cur_node->right;
	else {
		cur_node->right = add_node();
		cur_node = cur_node->right;
	}
}

void tree::insert_one_rand() {
	if(cur_node->right)
		cur_node = cur_node->right;
	else
		return;
}

void tree::insert_zero_main() {
	if(cur_node->left) 
		cur_node = cur_node->left;
	else {
		cur_node->left = add_node();
		cur_node = cur_node->left;
	}	
}

void tree::insert_zero_rand() {
	if(cur_node->left) 
		cur_node = cur_node->left;
	else
		return;
}

void tree::update_one_main(int val) {
	Leaf * leaf;
	if(cur_node->right) {
		leaf = (Leaf*)cur_node->right;
		update_leaf(leaf, val);
	}
	else {
		cur_node->right = add_leaf();
		leaf = (Leaf*)cur_node->right;
	}
}

void tree::update_one_rand(int val) {
	Leaf * leaf;
	if(cur_node->right) {
		leaf = (Leaf*)cur_node->right;
		update_leaf(leaf, val);
	}
	else
		return;
}

void tree::update_zero_main(int val) {
	Leaf * leaf;
	if(cur_node->left) {
		leaf = (Leaf*)cur_node->left;
		update_leaf(leaf, val);
	}
	else {
		cur_node->left = add_leaf();
		leaf = (Leaf*)cur_node->left;
	}
}

void tree::update_zero_rand(int val) {
	Leaf * leaf;
	if(cur_node->left) {
		leaf = (Leaf*)cur_node->left;
		update_leaf(leaf, val);
	}
	else
		return;
}

double * tree::extract() {
	C = new double[get_leafnum()];
	head = 0;
	DFS(return_root());
	return C;
}

void tree::DFS(Node * cur) {
	if(!cur->left && !cur->right) {
		Leaf * cur_leaf = (Leaf *)cur;
		if(cur_leaf->count > 0) {
			C[head] = leaf->count;
			head++;
			cur_leaf->count = 0;
		}
		return;
	}
	if(cur->left)
		DFS(cur->left);
	if(cur->right)
		DFS(cur->right);
}

value tree::destroy() {
	register int i;
	v.C = new double[get_leafnum()];
	v.adj_str = new char * [get_leafnum()];
	for(i = 0; i < get_leafnum(); i++)
		v.adj_str[i] = new char[subgraphSize*(subgraphSize-1)];
	head = 0;
	char * str = new char[subgraphSize*(subgraphSize-1)];
	DFS_value(return_root(), str, 0);
	return v;
}

void tree::DFS_value(Node * cur, char * str, int lev) {
	if(!cur->left && !cur->right) {
		Leaf * leaf = (Leaf *)cur;
		v.C[head] = leaf->count;
		v.adj_str[head] = str;
		head++;
		leaf->count = 0;
		return;
	}
	if(cur->left) {
		str[lev] = 0;
		DFS_value(cur->left, str, lev+1);
	}
	if(cur->right) {
		str[lev] = 1;
		DFS_value(cur->right, str, lev+1);
	}
}


