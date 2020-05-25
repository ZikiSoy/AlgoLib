
// $Treap 平衡树
// Treap
// split(root, left, right, interval)
// 0<=interval<=len(root)
// root = merge(left, right)
namespace Treap {
const int N = 5e5 + 10;
using int64 = long long;
using value_type = int64;
struct Node;
inline int len(Node* x) { return x ? x->s : 0; }
// inline value_type S(Node* x) { return x ? x->sum : 0; }
struct Node {
  // left, right, father
  Node *l, *r, *fa;
  // seed, size, boolReverse
  int w, s, rev;
  // value, lazyPushFlag
  int v, sav;
  // addmore
  // int sum;
  void init(int t = 0) {
    l = r = fa = NULL;
    w = rand();
    v = t;
    s = 1;
    sav = 0;
  }
  void update() {
    s = len(l) + len(r) + 1;
    /*sum=S(l)+S(r)+sav*s+v;*/
  }
  void push_down() {
    if (sav) {
      v += sav;
      if (l) {
        l->sav += sav; /*l->sum+=sav*len(l);*/
      }
      if (r) {
        r->sav += sav; /*r->sum+=sav*len(r);*/
      }
      sav = 0;
    }
    if (rev) {
      swap(l, r);
      if (l) l->rev ^= 1;
      if (r) r->rev ^= 1;
      rev = 0;
    }
  }
} p[N];
//0<=size<=len(p)
void split(Node *p, Node *&a, Node *&b, int sz) {
  p->push_down();
  if (sz <= 0)
    a = NULL, b = p;
  else if (sz >= len(p))
    a = p, b = NULL;
  else if (sz <= len(p->l)) {
    split(p->l, a, p->l, sz), b = p;
    a->fa = NULL;
    if (p->l) p->l->fa = p;
  } else {
    split(p->r, p->r, b, sz - len(p->l) - 1), a = p;
    b->fa = NULL;
    if (p->r) p->r->fa = p;
  }
  p->update();
}
// node a, b could be NULL
Node *merge(Node *a, Node *b) {
  if (!a || !b) return a ? a : b;
  a->push_down();
  b->push_down();
  if (a->w > b->w) {
    a->r = merge(a->r, b), a->r->fa = a;
    a->update();
    return a;
  } else {
    b->l = merge(a, b->l), b->l->fa = b;
    b->update();
    return b;
  }
}

// 返回 a 在二叉平搜索中是排在第几位
// 返回值为 0~len(root)-1
int getpos(Node *a) {
  int ret = len(a->l);
  if (a->rev) ret = len(a) - 1 - ret;
  for (Node *b = a->fa; b; a = b, b = b->fa) {
    if (a == b->r) ret += len(b->l) + 1;
    if (b->rev) ret = len(b) - 1 - ret;
  }
  return ret;
}

// // Give number from 0~len(root)-1, return node
// Node* getpos(Node* root, int pos) {
//   int ret = len(a->l);
//   if (a->rev) ret = len(a) - 1 - ret;
//   for (Node *b = a->fa; b; a = b, b = b->fa) {
//     if (a == b->r) ret += len(b->l) + 1;
//     if (b->rev) ret = len(b) - 1 - ret;
//   }
//   return ret;
// }
// // Delete node from position
// void delete_pos(Node* root, int pos) {
// }
// 递归调试 可能影响结构
void trace(Node *a) {
  if (!a) return;
  a->push_down();
  trace(a->l);
  printf("%d\n", a->v);
  trace(a->r);
}
}  // namespace Treap