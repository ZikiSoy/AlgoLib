// $Dinic 最大流
struct MaxFlow {
  constexpr int N = 2005;
  constexpr int E = N * 100;
  typedef long long Cost;

  int nv, ec, pnt[E], nxt[E];
  int que[N], lev[N], head[N];
  int flow, cap[E];
  int src, snk;
  // cap invcap
  void add(int u, int v, int c, int invf) {
    pnt[ec] = v;
    cap[ec] = c;
    nxt[ec] = head[u];
    head[u] = ec++;
    swap(u, v);
    pnt[ec] = v;
    cap[ec] = invf;
    nxt[ec] = head[u];
    head[u] = ec++;
  }
  void init(int _nv = N) {
    nv = _nv;
    memset(head, -1, sizeof(int) * (nv));
    ec = 0;
  }
  int _bfs() {  //bfs level
    int l = 0, r = 0, u;
    memset(lev, -1, sizeof(int) * (nv));
    lev[src] = 0;
    que[r++] = src;
    while (l < r) {
      u = que[l++];
      for (int p = head[u]; ~p; p = nxt[p])
        if (cap[p] && lev[pnt[p]] == -1) {
          lev[(que[r++] = pnt[p])] = lev[u] + 1;
          if (pnt[p] == snk) return 1;
        }
    }
    return 0;
  }

  int _find(int u, int in) {
    if (u == snk) return in;
    int t, w = 0;
    for (int p = head[u]; p >= 0 && w < in; p = nxt[p]) {
      if (cap[p] > 0 && lev[pnt[p]] == lev[u] + 1) {
        if (t = _find(pnt[p], min(cap[p], in - w))) {
          cap[p] -= t;
          cap[p ^ 1] += t;
          w += t;
        }
      }
    }
    if (w < in) lev[u] = -1;
    return w;
  }
  int dinic(int S, int T) {
    src = S, snk = T;
    int t, res = 0;
    while (_bfs())
      while ((t = _find(src, inf))) res += t;
    return res;
  }
} g;