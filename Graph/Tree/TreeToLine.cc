// https://www.codechef.com/submit/CHEFTRIP
#include <bits/stdc++.h>
 
// headers {{{
using namespace std;
// using namespace rel_ops;
 
typedef long long int64;
typedef unsigned long long uint64;
const double pi=acos(-1.0);
const double eps=1e-11;
const int inf=0x7FFFFFFF;
template<class T> inline bool checkmin(T &a,T b){return b<a?a=b,1:0;}
template<class T> inline bool checkmax(T &a,T b){return b>a?a=b,1:0;}
template<class T> inline T sqr(T x){return x*x;}
typedef pair<int,int> PII;
typedef vector<int> VI;
typedef vector<PII> VII;
#define MP(A,B) make_pair(A,B)
#define PB(X) push_back(X)
#define mem(a,b) memset(a, b, sizeof(a))
#define clr(a) memset(a, 0, sizeof(a))
#define rep(i,n) for(int i=0; i<n; i++)
#define repit(i,v) for(typeof(v.begin()) i=v.begin(); i!=v.end(); i++)
#define iter(v) typeof(v.begin())
#define ff first
#define ss second
#ifdef LOCAL
#define dbg(args...) printf(args); //##__VA_ARGS__
#define dout cout
#define out(x) (cout<<#x<<": "<<x<<endl)
template<class T>void show(T a, int n){for(int i=0; i<n; ++i) cout<<a[i]<<' '; cout<<endl;}
template<class T>void show(T a, int r, int l){for(int i=0; i<r; ++i)show(a[i],l);cout<<endl;}
#else
#define dbg(...)
#define dout if(true);else cout
#define out(...)
#define show(...)
#endif
// }}}

// $树链剖分
// 树链剖分模版 将树映射到链
namespace TTL {
const int MAXN = 1e5 + 10;
//==================================================
int to[MAXN * 2], nxt[MAXN * 2], head[MAXN];
int son[MAXN], fa[MAXN], siz[MAXN],
    dep[MAXN], top[MAXN];
int id[MAXN], ec, cnt;  //Need initialize
void init_tree(int n = MAXN) {
  memset(head, -1, sizeof(int)*n);
  ec = 0;
  cnt = 0;
}
inline void add(int u, int v) {
  to[ec] = v; nxt[ec] = head[u]; head[u] = ec++;
  swap(u, v);
  to[ec] = v; nxt[ec] = head[u]; head[u] = ec++;
}
void dfs(int u) {
  siz[u] = 1;
  son[u] = -1;  //dep[u] = 0; fa[u] = -1;
  for (int i = head[u]; i != -1; i = nxt[i]) {
    int &v = to[i];
    if (v != fa[u]) {
      fa[v] = u;
      dep[v] = dep[u] + 1;
      dfs(v);
      siz[u] += siz[v];
      if (son[u] == -1 || siz[v] > siz[son[u]])
        son[u] = v;
    }
  }
}
//from top->son 1.....a;
// a+1..b,
//call maketree(root,root)
void maketree(int u, int tp) {
  //maketree(thisnode,topnode)
  id[u] = cnt++;
  top[u] = tp;
  if (son[u] != -1) maketree(son[u], tp);
  for (int i = head[u]; i != -1; i = nxt[i]) {
    if (to[i] != fa[u] && to[i] != son[u])
      maketree(to[i], to[i]);
  }
}
int lca(int u, int v) {
  int fu = top[u], fv = top[v];
  while (fu != fv) {
    if (dep[fu] > dep[fv])
      u = fa[fu], fu = top[u];
    else
      v = fa[fv], fv = top[v];
  }
  return dep[u] < dep[v] ? u : v;
}
// 得到LogN个链。每个链first->second连续，可在线段树上操作
PII path[MAXN], pv[MAXN];
int path_cnt;
int getpath(int u, int v) {
  vector<PII> up, vp;
  int fu = top[u], fv = top[v];
  int cu = 0, cv = 0;
  while (fu != fv) {
      // out (u<<' '<<v<<' ' <<fu<<' '<<fv);
    if (dep[fu] > dep[fv]) {
      path[cu++] = MP(id[u], id[fu]);
      u = fa[fu], fu = top[u];
    }
    else {
      pv[cv++] = MP(id[fv], id[v]);
      v = fa[fv], fv = top[v];
    }
  }
  path[cu++] = MP(id[u], id[v]);
  for (int i = cv-1; i>=0; i--) path[cu++] = pv[i];
  path_cnt=cu;
  return cu;
}
};  // namespace TreeToLine

struct Node {
    int lv, rv, mx;
    bool invalid;
    Node rev() {
        Node ret = *this;
        swap(ret.lv, ret.rv);
        return ret;
    }
    void print() {
        dout << "Value " << lv << ' ' << mx << ' ' << rv << ' ' << invalid << endl;
    }
};
Node join(const Node& lhs, const Node& rhs) {
    Node ret;
    ret.invalid = lhs.invalid || rhs.invalid;
    if (ret.invalid) return ret;
    ret.lv = lhs.lv; ret.rv = rhs.rv;
    ret.mx = max(lhs.mx, rhs.mx);
    if (lhs.mx > lhs.rv && ((rhs.mx >= lhs.rv) || rhs.mx > rhs.lv))  ret.invalid = true;
    if (lhs.rv == rhs.lv) ret.invalid = true;
    if (lhs.mx == lhs.rv && rhs.mx > rhs.lv && rhs.lv <= lhs.rv) {
        ret.invalid = true;
    }
    return std::move(ret);
}
const int N = 1e5+10;
int h[N];
struct SegNode {
    int l, r;
    SegNode *lc, *rc;
    Node v;
    // void init(int ll, int rr);
    void init(const int& ll, const int& rr);
    Node query(int ll, int rr);
};
SegNode root, pool[N*2];
int cnt = 0;
void SegNode::init(const int& ll, const int& rr) {
    l = ll; r = rr;
    if (ll == rr) {
        v.lv = v.rv = v.mx = h[ll];
        v.invalid = false;
        dout << ll << ' ' << rr; v.print();
        return;
    }
    int mid = (ll+rr)/2;
    lc = pool+cnt++;
    rc = pool+cnt++;
    lc->init(ll, mid);
    rc->init(mid+1, rr);
    v = join(lc->v, rc->v);
    dout << ll << ' ' << rr; v.print();
}
Node SegNode::query(int ll, int rr) {
    if (ll <= l && rr >= r) return v;
    int mid = (l+r)/2;
    if (ll > mid) return rc->query(ll, rr);
    if (rr <= mid) return lc->query(ll, rr);
    return join(lc->query(ll, rr), rc->query(ll, rr));
}

Node query() {
    auto &p = TTL::path[0];
    Node tmp;
    if (p.second >= p.first) tmp = root.query(p.first, p.second);
    else tmp = root.query(p.second, p.first).rev();
    if (tmp.invalid) return move(tmp);
    for (int i = 1; i < TTL::path_cnt; i++) {
        auto &p = TTL::path[i];
        if (p.second >= p.first) tmp = join(tmp, root.query(p.first, p.second));
        else tmp = join(tmp, root.query(p.second, p.first).rev());
        if (tmp.invalid) return move(tmp);
    }
    return tmp;
    // ans.print();
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    cout.tie(NULL); 
    int T; cin >> T; rep(cas, T) {
        int N, Q;
        cin >> N >> Q;
        TTL::init_tree(N+1);
        out("init_tree()");
        rep(i, N-1) {
            int u, v;
            cin >> u >> v;
            TTL::add(u, v);
        }
        out("dfs()");
        TTL::dfs(1);
        out("maketree");
        TTL::maketree(1, 1);
        for (int i = 1; i<=N; i++) {
            dout << "$ " << i << ' ' <<TTL::id[i] <<' ' << TTL::fa[i] << ' ' << TTL::dep[i] << endl;
        }
        rep(i, N) {
            cin >> h[TTL::id[i+1]];
        }
        cnt = 0;
        out("init");
        root.init(0, N-1);
        string ans;
        rep(i, Q) {
            int x, y;
            cin >> x >> y;
            out("path");
            auto v = TTL::getpath(x, y);
            auto p = query();
            if (p.invalid) ans += '0';
            else ans += '1';
            out(ans);
        }
        cout << ans << endl;
    }
    return 0;
}

