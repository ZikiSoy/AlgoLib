// $.vimrc
// $vimrc
// set sw=4 ts=4 nobk ai cin
// set mp=g++\ -g\ -DLOCAL\ -O2\ -o\ %<\ %
// map <F5> :w<CR> :make<CR>
// map <F6> :!gnome-terminal -x ./%<<CR>
// map <F7> :!gnome-terminal -x gdb %<<CR>
// set gfn=Courier_New:h14
// set fencs=utf-8,gb2312,gbk,gb18030
// $header
#include <bits/stdc++.h>
using namespace std;
typedef long long int64;
#ifdef LOCAL
#define dbg(args...) printf(args);  //##__VA_ARGS__
#define dout cout
#define out(x) (cout << #x << ": " << x << endl)
template <class T>
void show(T a, int n) {
  for (int i = 0; i < n; ++i) cout << a[i] << ' ';
  cout << endl;
}
template <class T>
void show(T a, int r, int l) {
  for (int i = 0; i < r; ++i) show(a[i], l);
  cout << endl;
}
#else
#define dbg(...)
#define dout \
  if (true)  \
    ;        \
  else       \
    cout
#define out(...)
#define show(...)
#endif
const int inf = INT32_MAX/2;
#define rep(i, n) for (int i = )

// $$数学
// $质数筛法
#define Z 1000000
// factor primes phi(x)欧拉函数 prime_count
int f[Z], prime[Z], phi[Z], prc;
void GenPrime() {
  int i, p, tmp;
  for (p = 2; p < Z; ++p) {
    if (!f[p]) f[p] = p, prime[prc++] = p, phi[p] = p - 1;
    for (i = 0; i < prc && prime[i] <= f[p] && (tmp = prime[i] * p) < Z; ++i) {
      f[tmp] = prime[i];
      if (prime[i] == f[p])
        phi[tmp] = phi[p] * prime[i];
      else
        phi[tmp] = phi[p] * (prime[i] - 1);
    }
  }
}
// $Miller_Rabin + Pollard_Rho
const int MAXT = 100;
const int MAXN = 30;
int64 len, dig, limit;
int64 random() {
  int64 a = rand() * rand();
  return a * a;
}
bool Miller_Rabin(int64 n) {
  if (n < 2) return false;
  if (n == 2) return true;
  if (!(n & 1)) return false;
  int64 i, j, k, m, a;
  for (m = n - 1, k = 0; !(m & 1); m >>= 1, k++)
    ;
  for (i = 0; i < MAXT; i++) {
    a = power(random() % (n - 1) + 1, m, n);
    if (a == 1) continue;
    for (j = 0; j < k; j++) {
      if (a == n - 1) break;
      a = mul(a, a, n);
    }
    if (j == k) return false;
  }
  return true;
}
int64 gcd(int64 a, int64 b) {
  return b ? gcd(b, a % b) : a;
}
int64 f(int64 x, int64 n) {
  return (mul(x, x, n) + 1) % n;
}
int64 Pollard_Rho(int64 n) {
  if (n <= 2) return 0;
  if (!(n & 1)) return 2;
  for (int64 i = 1; i < MAXT; i++) {
    int64 x = random() % n;
    int64 xx = f(x, n);
    int64 p = gcd((xx + n - x) % n, n);
    while (p == 1) {
      x = f(x, n);
      xx = f(f(xx, n), n);
      p = gcd((xx + n - x) % n, n) % n;
    }
    if (p) return p;
  }
  return 0;
}
int64 factor[MAXN], m;
int64 Prime(int64 a) {
  if (Miller_Rabin(a)) return 0;
  int64 t = Pollard_Rho(a);
  int64 p = Prime(t);
  if (p) return p;
  return t;
}
//所有素约数保存于 factor[]中, m 为其个数, 返回分解后的素因子个数
int factorize(int64 x, int64 factor[]) {
  int m = 0;
  while (x > 1) {
    if (Miller_Rabin(x)) break;
    int64 t = Prime(x);
    factor[m++] = t;
    x /= t;
  }
  if (x > 0) factor[m++] = x;
  //所有素约数保存于 factor[]中,m 为其个数
  return m;
}

namespace romberg {
// $龙贝格积分公式
template <typename T>
T fun(T x) {  //被积函数
  return 3 * x * x;
}

template <typename T>
T romberg(T a, T b, T eps = 1e-5) {
  vector<T> R;
  int k = -1, pow2 = 1;
  double A;
  T r = 0.5 * (b - a) * (fun(a) + fun(b));
  R.push_back(r);
  do {
    A = R[0];
    k += 1;
    r = 0.0;
    for (int i = 0; i < pow2; i++)
      r += fun(a + (b - a) * (i + 0.5) / pow2);
    r *= (b - a) / (2.0 * pow2);
    r += 0.5 * R[k];
    R.push_back(r);
    double pow4 = 4.0;
    for (int m = 0; m <= k; m++, pow4 *= 4.0)
      R[k - m] = (pow4 * R[k + 1 - m] - R[k - m]) / (pow4 - 1);
    pow2 *= 2;
  } while (fabs(R[0] - A) > eps);
  return R[0];
}
int main() {
  printf("%lf\n", romberg(0.0, 5.0, 1e-5));
}
}  // namespace romberg

// $高斯消元
// 高斯消元——一般矩阵
namespace gauss {
double mat[MAXN][MAXN + 1], answer[MAXN];
// 返回矩阵的秩
// mat[n][[[[[[[[[[[[[[[[[[m]]]]]]]]]]]]]]]]]]
int gauss(int n, int m) {
  int k = 0, t = 0;
  for (; k < m && t < n; ++k, ++t) {
    // k for column; t for row
    int dec = t;
    for (int i = t + 1; i < n; ++i)
      if (fabs(mat[i][k]) > fabs(mat[dec][k]))
        dec = i;
    if (fabs(mat[dec][k]) < 1e-8) {
      --t;
      continue;
    }
    if (dec != t) {
      for (int j = k; j <= m; ++j)
        swap(mat[t][j], mat[dec][j]);
    }
    for (int i = t + 1; i < n; ++i) {
      double r = mat[i][k] / mat[t][k];
      mat[i][k] = 0.0;
      for (int j = k + 1; j <= m; ++j)
        mat[i][j] -= mat[t][j] * r;
    }
  }
  for (int i = t; i < n; ++i)
    if (fabs(mat[i][m]) > 1e-8)
      return -1;
  if (t == m)
    for (int k = m - 1; k >= 0; --k) {
      answer[k] = mat[k][m];
      for (int i = k + 1; i < m; ++i)
        answer[k] -= mat[k][i] * answer[i];
      answer[k] /= mat[k][k];
    }

  return t;
}
}  // namespace gauss

// $组合函数
//错位排列数
//D[n]表示 n 个相异的元素排成一排 a1, a2, ..., an,
//且 ai(i = 1, 2, ..., n)不在第 i 位的排列个数
int64 D[MAXN + 1];
void Derangement(int n) {
  D[1] = 0, D[2] = 1;
  for (int i = 3; i <= n; i++) {
    D[i] = (i - 1) * (D[i - 1] + D[i - 2]);
  }
}
//第一类斯特灵数
//StirlingS1[n, m]是有正负的,其绝对值是
// 包含 n 个元素的集合分成 m 个环排列的方法数目
//又即有 n 个人分成 m 组,
// 每组内再按特定顺序围圈的分组方法的数目
//S1(n, m) = (n - 1) * S1(n - 1, m)
// + S1(n - 1, m - 1)
//S1(n, 0) = 0, S1(1, 1) = 1, S1(4, 2) = 11
int64 S1[MAXN + 1][MAXN + 1];
void StirlingS1(int n) {
  S1[0][0] = 1;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= i; j++) {
      S1[i][j] = (i - 1) * S1[i - 1][j] + S1[i - 1][j - 1];
    }
  }
}
//第二类斯特林数及贝尔数
//StirlingS2[n, m]给出把 n 个元素的集合分到 m 个
// 非空子集的分法的数目 也叫集合的分拆
//B[n]是基数为 n 的集合的划分方法的数目
//1, 1, 2, 5, 15, 52, 203, 877, 4140,
// 21147, 115975, 678570, 4213597
int64 S2[MAXN + 1][MAXN + 1];
int64 B[MAXN + 1];
void StirlingS2(int n) {
  S2[0][0] = B[0] = 1;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= i; j++) {
      S2[i][j] = j * S2[i - 1][j] + S2[i - 1][j - 1];
      B[i] += S2[i][j];
    }
  }
}
//卡特兰数
//Catalan[n]是 n 个+1 和 n 个-1 构成 2n 项,
// 其部分和满足 a1 + a2 + ... + ak >= 0 的序列个数
//Catalan[n] = 1 / (n + 1) * C(2 * n, n)
// = C(2 * n, n) - C(2 * n, n + 1)
//Catalan[n] = (4 * n - 2) /
// (n + 1) * C[n - 1], Catalan[1] = 1
//Catalan[n] = sigma(catalan[i] * catalan[n - 1 - i])
// (0 <= i < n)
//前几个 Catalan 数是:
//1, 1, 2, 5, 14, 42, 132, 429, 1430,
// 4862, 16796, 58786, 208012, 742900
//广义卡特兰数

//求由 n 个 1、m 个 0 组成,并且任意前缀中 1 的个数不少于 0 的个数的字符串的个数是
// C(n + m, n) - C(n + m, n + 1)
//从(0, 0)点走到(m, n)点的方案是
// C(n + m, n) - C(n + m, n + 1) =
// (n - m + 1) / (n + 1) * C(n + m, m)
int64 catalan[MAXN + 1];
void Catalan(int n) {
  catalan[0] = 1;
  for (int i = 1; i <= n; i++) {
    catalan[i] = 0;
    for (int j = 0; j < i; j++) {
      catalan[i] += catalan[j] * catalan[i - 1 - j];
    }
  }
}
//分拆函数 P
//PartP[n][m]表示数 n 分拆最大部分为 m 的分拆个数(n >= m)
//分拆函数 P 前几项是:
//1, 1, 2, 3, 5, 7, 11, 15, 22, 30,
// 42, 56, 77, 101, 135, 176, 231, 297
int64 PartP[MAXN + 1][MAXN + 1];
int64 PartitionP[MAXN + 1];
void Partition_P(int n) {
  PartP[0][0] = PartitionP[0] = 1;
  for (int i = 1; i <= n; i++) {
    PartitionP[i] = 0;
    for (int j = 1; j <= i; j++) {
      PartP[i][j] = PartP[i - 1][j - 1] + PartP[i -
                                                j][j];
      PartitionP[i] += PartP[i][j];
    }
  }
}
//分拆函数 Q
//PartQ[n][m]表示数 n 分拆为 m 个不同的正整数的方案数目
// G[x] = Product_{m=1..inf} (1 + x^m)的展开式 x^n 项系数;
//分拆函数 Q 前几项是:
//1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12,
// 15, 18, 22, 27, 32, 38, 46, 54, 64, 76, 89
int64 PartQ[MAXN + 1][MAXN + 1];
int64 PartitionQ[MAXN + 1];
void Partition_Q(int n) {
  Partition_P(n);
  for (int i = 0; i <= n; i++) {
    PartitionQ[i] = 0;
    for (int j = 0; j <= i; j++) {
      if (i - j * (j - 1) / 2 < 0) break;
      PartQ[i][j] = PartP[i - j * (j - 1) / 2][j];
      PartitionQ[i] += PartQ[i][j];
    }
  }
}
//自然数前 n 项 k 次方和
//1^k+2^k+...+n^k
// = sigma(P[k][i] * C[n + 1][i + 1]);//i = 1...k;
// = sigma((sigma((-1)^j * C[i][j] * (i - j)^k))
// * C[n + 1][i + 1]); //i = 1..k, j = 0..i - 1;
//其中 C[i][j]为组合数,
// C[i][j] = C[i - 1][j] + C[i - 1][j - 1];
//P[i][0] = 0; P[i][1] = 1;
//P[i][i] = 1 * 2 * ... * i = i !;
//P[i][j] = j * (P[i - 1][j - 1] + P[i - 1][j]);

// (1 <= j <= i)
int64 PS[MAXN + 1][MAXN + 1];
void PowerSum(int n) {
  // Combination(n);
  PS[0][1] = 1;
  for (int i = 0; i <= n; i++) {
    for (int j = 1; j <= i; j++) {
      PS[i][j] = (j * (PS[i - 1][j - 1] + PS[i - 1][j]));
    }
  }
  //完整的加上逆元、快速幂模块用公式求即可
}
// 根号方程
// (sqrt(a) + sqrt(b))^(2n) + (sqrt(a) - sqrt(b))^(2n)
// 递推 Z[n+1] = 2 * (a + b) * Z[n]
// - (a - b)^2 * z[n - 1];
// Z0 = 2, Z1 = 2 * (a + b)
int main() {
  Partition_Q(MAXN);
  int cas, n;
  scanf("%d", &cas);
  while (cas--) {
    scanf("%d", &n);
    printf("%lld\n", PartitionQ[n]);
  }
}
// $$图论
// $Dinic 最大流
struct MaxFlow {
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
      while (t = _find(src, inf)) res += t;
    return res;
  }
} g;
// $最小费用流最大流
// 可解:最小费用流、最小费用最大流
// 引流:某些边必须要走的时候 将这些边的费用*-1e9
const int N = 205;
const int E = N * N;
typedef long long Cost;
struct MincostFlow {
  int nv, ec, pnt[E], nxt[E];
  int vis[N], que[N], head[N], pv[N], pe[N];
  int flow, cap[E];
  Cost cost, dis[E], d[N];
  void add(int u, int v, int c, Cost w) {
    pnt[ec] = v;
    cap[ec] = c;
    dis[ec] = +w;
    nxt[ec] = head[u];
    head[u] = ec++;
    swap(u, v);
    pnt[ec] = v;
    cap[ec] = 0;
    dis[ec] = -w;
    nxt[ec] = head[u];
    head[u] = ec++;
  }
  bool spfa(int src, int snk, int n) {
    int u, v, k, l, r;
    memset(pv, -1, sizeof(int) * (n));
    memset(vis, 0, sizeof(int) * (n));
    memset(d, 63, sizeof(Cost) * (n));
    d[src] = 0;
    vis[src] = 1;
    pv[src] = src;
    for (l = 0, r = 1, que[0] = src; l != r;) {
      u = que[l++];
      vis[u] = 0;
      if (N == l) l = 0;
      for (k = head[u]; ~k; k = nxt[k]) {
        v = pnt[k];
        if (cap[k] && dis[k] + d[u] < d[v]) {
          d[v] = dis[k] + d[u];
          if (!vis[v]) {
            vis[v] = 1;
            que[r++] = v;
            if (r == N) r = 0;
          }
          pv[v] = u;
          pe[v] = k;
        }
      }
    }
    return pv[snk] != -1;
  }
  int mincost(int src, int snk, int n) {
    cost = 0, flow = 0;
    int mxf, k;
    while (spfa(src, snk, n)) {
      if (d[snk] > 0) break;  //mincost but not maxflow
      for (k = snk, mxf = inf; k != src; k = pv[k])
        if (cap[pe[k]] < mxf) mxf = cap[pe[k]];
      flow += mxf;
      cost += d[snk] * mxf;
      for (k = snk; k != src; k = pv[k]) {
        cap[pe[k]] -= mxf;
        cap[pe[k] ^ 1] += mxf;
      }
    }
    return flow;
  }
  void init() {
    memset(head, -1, sizeof(head));
    ec = 0;
  }
} g;
// $二分图最大匹配
namespace bimatch {
const int MAXN = 1005, MAXM = 1005;
int nx, ny, m, ans;  //nx,ny 分别为二分图两边节点的个数,两边的节点分别用 1..nx,1..ny 编号,m 为边数
bool g[MAXN][MAXM];  //图 G 邻接矩阵 g[x][y]
bool y[MAXM];        //Y 集合中点 i 访问标记
int link[MAXM];      //link[y]表示当前与 y 节点相邻的 x 节点
void init() {
  int x, y;
  memset(g, 0, sizeof(g));
  ans = 0;
  scanf("%d%d%d", &nx, &ny, &m);
  for (int i = 1; i <= m; i++) {
    scanf("%d%d", &x, &y);
    g[x][y] = true;
  }
}
bool find(int x) {  //是否存在 X 集合中节点 x 开始的增广路
  for (int i = 1; i <= ny; i++)
    if (g[x][i] && !y[i]) {
      y[i] = true;
      if (link[i] == -1 || find(link[i])) {
        link[i] = x;
        return true;
      }
    }
  return false;
}
int MaximumMatch() {
  int ret = 0;
  memset(link, -1, sizeof(link));
  for (int i = 1; i <= nx; i++) {
    memset(y, 0, sizeof(y));
    if (find(i))
      ret++;
  }
  return ret;
}
}  // namespace bimatch

// $部分点生成树 斯坦纳树 Steiner Tree
int steiner(int n, int m, int p[]) {
  //0~n 共 n 个点里面选出 p 中的 m 个点

  memset(dp, 63, sizeof(dp));
  for (int i = 0; i < n; i++)
    SPFA(i);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      dp[1 << i][j] = dis[p[i]][j];
  for (int i = 1; i < (1 << m); i++)
    if (((i - 1) & i) != 0) {  //二进制不止一个 1
      for (int j = 0; j < n; j++) {
        dp[i][j] = MAXV;
        for (int k = (i - 1) & i; k > 0; k = (k - 1) & i)  //枚举子集
          dp[i][j] = min(dp[i][j], dp[k][j] + dp[i ^ k][j]);
      }
      for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
          dp[i][j] = min(dp[i][j], dp[i][k] + dis[k][j]);
    }
  return dp[(1 << m) - 1][p[0]];
}

namespace treecount { 
// $生成树计数 Matrix-Tree 定理
// (Kirchhoff 矩阵-树定理)
// 给一个无向图 问有多少种生成树
// A[i][j] = bool(ij 有边)
// A[i][i] = -出度
// 求 A 的 n - 1 阶行列式

// $ Prüfer 序列
#include <iostream>
using namespace std;
typedef long long int64;
const int MAXN = 252;
// 有标号无向图
// LG[n] = 2^(C[n][2])
// 点度全偶时 LG2[N] = 2^(C[n - 1][2])

// 有标号有根树计数
// LRT[n] = n^(n - 1)

// 有标号无根树计数
// LNRT[n] = n^(n - 2)

// 无标号有根树计数
// RT[n]表示 n 个相同的顶点能构成有根树的方案数
// 前几项为 1 2 4 9 20 48 115 286 719 1842 4766 12486 32973 87811 235381
int64 RT[MAXN], RTS[MAXN][MAXN];
void RootedTree(int n) {
  RT[1] = 1;
  for (int i = 1; i < n; i++) {
    for (int j = 1; j <= i; j++) {
      RTS[i][j] = RTS[i - j][j] + RT[i + 1 - j];
      RT[i + 1] += j * RT[j] * RTS[i][j];
    }
    RT[i + 1] /= i;
  }
}

// $无标号无根树计数
// NRT[n]表示 n 个相同的顶点能构成无根树的方案数
// 前几项为 11 1 2 3 6 11 23 47 106 235 551 1301 3159 7741 19320 
int64 NRT[MAXN];
void NonRootedTree(int n) {
  RootedTree(n);

  for (int i = 1; i <= n; i++) {
    NRT[i] = RT[i];
    for (int j = 1; j <= i / 2; j++) {
      NRT[i] -= RT[j] * RT[i - j];
    }
    if (i % 2 == 0) {
      NRT[i] += (RT[i / 2] + 1) * RT[i / 2] / 2;
    }
  }
}
//无标号二叉树
//BT[n] = Catalan[n] = C[2 * n][n] / (n + 1)
//标号为 k 的点度为 vk 的无根树
//NRTK[n] = (n - 2)!/(all(vk-1)!)
//无标号毛毛虫
//Caterpiller[n] = 2^(n-4) + 2 ^ [(n-4)/2]下整
//有标号连通图
//LG[1] = LG[2]=1;
//LG[n]=2^(C[n][2]) – sigma(C[n - 1][i - 1] * LG[i] * 2 ^ (C[n - i][2]));
//即枚举最后一个点所在连通分量的大小。
int main() {
  NonRootedTree(MAXN);
  int cas, n;
  scanf("%d", &cas);
  while (cas--) {
    scanf("%d", &n);
    printf("%I64d\n", NRT[n]);
  }
}
}  // namespace treecount

// 最大权闭合图:
// 定义:一个有向图的闭合图 G=(V,E)是该有向图的一个点集,且该点集的
// 所有出边都还指向该点集。
// 即闭合图内的任意点的任意后继也一定在闭合图中。
// 给每个点 v 分配一个点权(任意实数,可正可负)。最大权闭合图,是一
// 个点权之和最大的闭合图。
// 闭合图的性质恰好反映了事件之间的必要条件的关系:一个事件发生,它
// 需要的所有前提都要发生。
// 下面通过构图,我们将最大权闭合图问题转化成最小割问题,即最大流问
// 题求解。
// 定义 W[I]代表顶点 I 的权值,新增源点 S,汇点 T。
// 1、若 W[I]>0,则 S 向 I 连一条容量为 W[I]的边。
// 2、若 W[I]<0,则 I 向 T 连一条容量为-W[I]的边。
// 3、原图中的边,容量设置为正无穷。
// 这样,最小割就对应了最大权闭合图,而总盈利-最大流就是权和。

// $无向图点双联通分量
// find_scc( int node_num ); 结束后 bccno 为每个点的双联
// 通分量标号 bcc 为每个双联通分量中的点
int pre[maxn], iscut[maxn], bccno[maxn], dfsclk, bcc_cnt;
VI G[maxn], bcc[maxn];
stack<PII> S;
int dfs(int u, int fa) {
  int lowu = pre[u] = ++dfsclk;
  int child = 0;
  for (int i = 0; i < G[u].size(); i++) {
    int v = G[u][i];

    PII e = MP(u, v);
    if (!pre[v]) {
      S.push(e);
      child++;
      int lowv = dfs(v, u);
      lowu = min(lowu, lowv);
      if (lowv >= pre[u]) {
        iscut[u] = true;
        bcc_cnt++;
        bcc[bcc_cnt].clear();
        for (;;) {
          PII x = S.top();
          S.pop();
          if (bccno[x.ff] != bcc_cnt) {
            bcc[bcc_cnt].push_back(x.ff);
            bccno[x.ff] = bcc_cnt;
          }
          if (bccno[x.ss] != bcc_cnt) {
            bcc[bcc_cnt].push_back(x.ss);
            bccno[x.ss] = bcc_cnt;
          }
          if (x.ff == u && x.ss == v) break;
        }
      } else if (pre[v] < pre[u] && v != fa) {
        S.push(e);
        lowu = min(lowu, pre[v]);  // 反向边
      }
    }
  }
  if (fa < 0 && child == 1) iscut[u] = 0;
  return lowu;
}
void find_bcc(int n) {
  mem(pre, 0);
  mem(iscut, 0);
  mem(bccno, 0);
  dfsclk = bcc_cnt = 0;
  for (int i = 0; i < n; i++)
    if (!pre[i]) dfs(i, -1);
  //memset(pre,0,sizeof(pre));
  //memset(iscut,0,sizeof(iscut));
}
// $无向图边双联通分量
// find_scc( int node_num );
// 结束后 bccno 为每个点的双联通分量标号
// bcc 为每个双联通分量中的点
int pre[maxn], iscut[maxn], bccno[maxn], dfsclk,
    bcc_cnt;
VI bcc[maxn];
int dfs(int u, int fe) {
  int lowu = pre[u] = ++dfsclk;
  for (int i = 0; ~i; i = nxt[i])
    if (i ^ fe ^ 1) {  //不是反向边
      int v = to[i];
      if (!pre[v]) {
        int lowv = dfs(v, i);
        if (lowv > pre[u]) {
          iscut[i] = iscut[i ^ 1] = true;
        }
      } else
        lowu = min(lowu, pre[v]);
    }
  return low[u] = lowu;
}
void dfs2(int u) {
  bccno[u] = bcc_cnt;
  bcc[bcc_cnt].push_back(u);
  for (int i = 0; ~i; i = nxt[i])
    if (!iscut[i]) {
      int v = to[i];
      if (!bccno[u]) dfs2(v);
    }
}
void find_bcc(int n) {
  mem(pre, 0);
  mem(iscut, 0);
  mem(bccno, 0);
  dfsclk = bcc_cnt = 0;
  for (int i = 0; i < n; i++)
    if (!pre[i]) dfs(i, -1);
  //memset(pre,0,sizeof(pre));
  //memset(iscut,0,sizeof(iscut));
}
// $Kosaraju 强连通分量
// forw 为正向边 back 为反向边
// 节点编号 0..n-1,调用 kosaraju(n)后,
// scc 为 scc 个数 num[i]为节点 scc 编号,
// sccv 存储每个 scc 所在的节点,scc 按拓扑从大到小排列。
typedef vector<int> VI;
const int MAXV = 1e4 + 10;
int dfn[MAXV], num[MAXV], scc, cnt;
VI forw[MAXV], back[MAXV];
VI sccv[MAXV];
void dfs(int u) {
  num[u] = -1;
  for (int i = 0; i < forw[u].size(); i++)
    if (num[forw[u][i]] != -1) {
      dfs(forw[u][i]);
    }
  dfn[cnt++] = u;
}
void ndfs(int u) {
  num[u] = scc;
  sccv[scc].push_back(u);
  for (int i = 0; i < back[u].size(); i++) {
    if (num[back[u][i]] == -1) {
      ndfs(back[u][i]);
    }
  }
}
void kosaraju(int n) {
  clr(num);
  scc = cnt = 0;
  for (int i = 0; i < n; i++)
    if (!num[i]) dfs(i);
  for (int i = cnt - 1; i >= 0; i--) {
    if (num[dfn[i]] == -1) {
      sccv[scc].clear();
      ndfs(dfn[i]);
      scc++;
    }
  }
}

// $图的最小环 O(n^3):
// dist[i][j]表示距离,每次用新点去更新环的大小
// 然后更新距离

// $2-SAT
// 求可行解、最小字典序可行解
// x,y 不能共存,则连边 x->y^1, y->x^1
// 最后如果 x,x^1 在同一个强连通分量内则无解
const int maxn = 1e4 + 10;
struct TwoSAT {
  int n;
  vector<int> G[maxn * 2];
  bool mark[maxn * 2];
  bool dfs(int x) {
    if (mark[x ^ 1]) return false;
    if (mark[x]) return true;
    mark[x] = true;
    for (int i = 0; i < G[x].size(); i++)
      if (!dfs(G[x][i])) return false;
    return true;
  }
  void init(int n) {
    this->n = n;
    for (int i = 0; i < n * 2; i++) G[i].clear();
    memset(mark, 0, sizeof(mark));
  }
  // x = xval or y = yval
  void addc(int x, int xval, int y, int yval) {
    x = x * 2 + xval;
    y = y * 2 + yval;
    G[x ^ 1].push_back(y);
    G[y ^ 1].push_back(x);
  }
  bool solve() {
    for (int i = 0; i < n * 2; i += 2)
      if (!mark[i] && !mark[i ^ 1])
        if (!dfs(i)) return false;
    return true;
  }
};
// $$数据结构:
// $KMP
void getFail(char *P, int m, int *f) {
  f[0] = f[1] = 0;
  for (int i = 1; i < m; i++) {
    int j = f[i];
    while (j && P[i] != P[j]) j = f[j];
    f[i + 1] = P[i] == P[j] ? j + 1 : 0;
  }
}
//T 长,P 短,f 是 fail 数组
int kmp(char *T, int n, char *P, int m, int *f) {
  getFail(P, m, f);
  int i, j = 0;
  rep(i, n) {
    while (j && P[j] != T[i]) j = f[j];
    if (P[j] == T[i]) j++;
    if (j == m) {
      return 1;
    }
  }
  return 0;
}

namespace SA {
// $后缀数组+RMQ
// sa[] 将后缀 0~n-1 放入 sa 数组按后缀字典序排序后的数组
// rank[] 后缀 i 在 sa 中的下标

// height[i] sa[i]和 sa[i-1]的最长公共前缀
// lcp(x, y) sa[x-1]和 sa[y]的最长公共前缀
const int N = 3e5 + 10;
int s[N];
int sa[N], t[N], t2[N], c[N], n;
int rank[N], height[N];
void build_sa(int m) {
  int i, *x = t, *y = t2, p;
  for (i = 0; i < m; i++) c[i] = 0;
  for (i = n - 1; i >= 0; i--) c[x[i] = s[i]]++;
  for (i = 1; i < m; i++) c[i] += c[i - 1];
  for (i = n - 1; i >= 0; i--) sa[--c[x[i]]] = i;
  for (int k = 1; k <= n; k <<= 1) {
    //x 第一关键字的排序数组 y 第二关键字的排序数组
    for (p = 0, i = n - 1; i >= n - k; i--) y[p++] = i;
    for (i = 0; i < n; i++)
      if (sa[i] >= k) y[p++] = sa[i] - k;
    //sa 数组为第一关键字
    for (i = 0; i < m; i++) c[i] = 0;
    for (i = 0; i < n; i++) c[x[y[i]]]++;
    for (i = 1; i < m; i++) c[i] += c[i - 1];
    for (i = n - 1; i >= 0; i--) sa[--c[x[y[i]]]] = y[i];
    swap(x, y);
    p = 1;
    x[sa[0]] = 0;
    for (i = 1; i < n; i++)
      x[sa[i]] = y[sa[i - 1]] == y[sa[i]] && y[sa[i - 1] + k] == y[sa[i] + k] ? p - 1 : p++;
    if (p >= n) break;
    m = p;
  }
}
void getHeight() {
  int i, k = 0;
  for (i = 0; i < n; i++) rank[sa[i]] = i;
  for (i = 0; i < n; i++) {
    if (k) k--;
    if (rank[i] == 0) continue;
    int j = sa[rank[i] - 1];
    while (s[i + k] == s[j + k]) k++;
    height[rank[i]] = k;
  }
}
int st[20][N], ln[N];
#define choose min
void initrmq(int n, int val[]) {
  int i, j, k, sk;
  ln[0] = -1;
  for (i = 1; i <= n; i++) ln[i] = ln[i >> 1] + 1;  //预处理
  for (i = 0; i < n; i++) st[0][i] = val[i];
  for (i = 1, k = 2; k < n; i++, k <<= 1)
    for (j = 0, sk = (k >> 1); sk < n; ++j, ++sk)
      st[i][j] = choose(st[i - 1][j], st[i - 1][sk]);
}
int lct(int x, int y) {  //choose the best in { val[x] ... val[y] }
  int bl = ln[y - x + 1];
  return choose(st[bl][x], st[bl][y - (1 << bl) + 1]);
}
// n = total_length; // build_sa(charset_size);
// getHeight(); // initrmq(total_length);
// using lct(x, y)
}; //  namespace SA
// =============== template ended ===============


// $AC 自动机
const int SIGMA_SIZE = 2;
const int MAXNODE = 11000;
const int MAXS = 150 + 10;
struct AhoCorasickAutomata {
  int ch[MAXNODE][SIGMA_SIZE];
  int f[MAXNODE];     // fail 函数
  int val[MAXNODE];   // 字符串结尾节点 val 非 0
  int last[MAXNODE];  // 输出链表的下一个结点
  int cnt[MAXS];
  int sz;
  void init() {
    sz = 1;
    memset(ch[0], 0, sizeof(ch[0]));
    memset(cnt, 0, sizeof(cnt));
  }
  int idx(char c) {
    return c == 'H';
  }
  // 插入字符串。v 必须非 0
  int insert(char *s, int v) {
    int u = 0, n = strlen(s);
    for (int i = 0; i < n; i++) {
      int c = idx(s[i]);
      if (!ch[u][c]) {
        memset(ch[sz], 0, sizeof(ch[sz]));
        val[sz] = 0;
        ch[u][c] = sz++;
      }
      u = ch[u][c];
    }
    val[u] = v;
    return u;
  }
  void getFail() {
    queue<int> q;
    f[0] = 0;
    // 初始化队列
    for (int c = 0; c < SIGMA_SIZE; c++) {
      int u = ch[0][c];
      if (u) {
        f[u] = 0;
        q.push(u);
        last[u] = 0;
      }
    }
    // 按 BFS 顺序计算 fail
    while (!q.empty()) {
      int r = q.front();
      q.pop();
      for (int c = 0; c < SIGMA_SIZE; c++) {
        int u = ch[r][c];
        if (!u) {
          ch[r][c] = ch[f[r]][c];
          continue;
        }
        q.push(u);
        int v = f[r];
        f[u] = ch[v][c];
        last[u] = val[f[u]] ? f[u] : last[f[u]];
      }
    }
  }
} ac;
// $Treap 平衡树
// Treap
// split(root, left, right, interval)
// 0<=interval<=len(root)
// root = merge(left, right)
namespace Treap {
const int N = 5e5 + 10;
#define len(x) (x ? x->s : 0)
#define S(x) (x ? x->sum : 0)
typedef long long int64;
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
  void up() {
    s = len(l) + len(r) + 1;
    /*sum=S(l)+S(r)+sav*s+v;*/
  }
  void dn() {
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
  p->dn();
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
  p->up();
}
// node a, b could be NULL
Node *merge(Node *a, Node *b) {
  if (!a || !b) return a ? a : b;
  a->dn();
  b->dn();
  if (a->w > b->w) {
    a->r = merge(a->r, b), a->r->fa = a;
    a->up();
    return a;
  } else {
    b->l = merge(a, b->l), b->l->fa = b;
    b->up();
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
// 递归调试 可能影响结构
void trace(Node *a) {
  if (!a) return;
  a->dn();
  trace(a->l);
  printf("%d\n", a->v);
  trace(a->r);
}
}  // namespace Treap

// $笛卡尔树
// 笛卡尔树(就是一棵 Treap 的 O(n)构造过程)
// zju 1985 一列宽为 1 高为 height[i]的木棍,求最大延伸面积,
// 利用了 T=i,H=height[i],然后父子之间夹住的区间......类似
// splay 的想法...
namespace CartesianTree {
#define maxn 110000
struct data {
  int T, H;
  data *ch[2];  //Tree,Heap
  void init(int _T, int _H) {
    T = _T;
    H = _H;
    ch[0] = ch[1] = NULL;
  }
  friend bool operator<(const data &a, const data &b) {
    return a.T < b.T;
  }
} dd[maxn], *root;
//节点 0~n-1
data *st[maxn];
int top;
void makeTree(int n) {
  top = 0;
  for (int i = 0; i < n; i++) {
    data *lson = NULL;
    while (top > 0) {
      if (dd[i].H < st[top - 1]->H) {
        lson = st[--top];
      } else
        break;
    }
    if (top == 0) {
      dd[i].ch[0] = lson;
      root = &dd[i];
    } else {
      st[top - 1]->ch[1] = &dd[i];
      dd[i].ch[0] = lson;
    }
    st[top++] = &dd[i];
  }
}
long long ans;

int dfs(data *p) {
  if (p == NULL) return 0;
  int num = 1;
  num += dfs(p->ch[0]);
  num += dfs(p->ch[1]);
  ans = max(ans, (long long)num * p->H);
  return num;
}
int main_example() {
  int n;
  while (scanf("%d", &n) && n) {
    for (int i = 0; i < n; i++) {
      int H;
      scanf("%d", &H);
      dd[i].init(i, H);
    }
    makeTree(n);
    ans = 0;
    dfs(root);
    printf("%lld\n", ans);
  }
}
};  // namespace CartesianTree
// $树链剖分
// 树链剖分模版 将树映射到链
namespace TreeToLine {
const int MAXN = 1e5 + 10;
//==================================================
int to[MAXN * 2], nxt[MAXN * 2], head[MAXN];
int son[MAXN], fa[MAXN], siz[MAXN],
    dep[MAXN], top[MAXN];
int id[MAXN], ec, cnt;  //Need initialize
void init_tree() {
  memset(head, -1, sizeof(head));
  ec = 0;
  cnt = 0;
}
inline void add(int u, int v) {
  to[ec] = v;
  nxt[ec] = head[u];
  head[u] = ec++;
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
};  // namespace TreeToLine


// $斜率优化

namespace XieLv {
// dp[p][x] = min(dp[p - 1][k - 1] +
//                (x - k) * a[x] - (sum[x] - sum[k]));
// G = dp[p][x] + sum[x] - a[x] * x;
// k = a[x];
// x = k;
// y = dp[p - 1][k] + sum[k];

const int N = 1e5 + 111;
int64 dd[N];
int64 sy[N];
int64 sum[N];
int64 dp[101][N];
int64 num[N];
struct point {
  int64 x, y;
  point(int64 x = 0, int64 y = 0) : x(x), y(y) {}
  point operator-(const point B) const {
    return point(x - B.x, y - B.y);
  }
} q[N];
inline int64 cross(point O, point A, point B) {
  point OA = A - O;
  point OB = B - O;
  return OA.x * OB.y - OB.x * OA.y;
}
int64 getx(int i) {
  return i;
}
int64 gety(int _p, int i) {
  return dp[_p - 1][i] + sum[i];
}
int64 getk(int i) {
  return sy[i];
}
int64 getg(int l, int i) {
  return q[l].y - getk(i) * q[l].x;
}
void slope(int _p, int m) {
  int l, r;
  l = r = 1;
  dp[_p - 1][0] = 0;
  q[1] = point(getx(0), gety(0, 0));
  int i, j;
  rep(i, m + 1) {
    point tmp(getx(i), gety(_p, i));
    for (j = r; j > l; j--) {
      if (cross(q[j - 1], q[j], tmp) >= 0)
        break;
      else
        r--;
    }
    q[++r] = tmp;
    for (j = l; j < r; j++) {
      if ((q[j + 1].y - q[j].y) < getk(i) * (q[j + 1].x -
                                             q[j].x))
        l++;
      else
        break;
    }
    dp[_p][i] = getg(l, i) - sum[i] + sy[i] * i;
  }
}
};  // namespace XieLv
// $各种数位计数:
// 2 行全排列 取对应的大的值的和大于等于 k 的组合数
// int dp[51][26][2501];//s,o,k s:空的列,o:空的单上和单
// 下
// 给一个串,要 a 到 b 连 c 到 d 的组合使得有 x 个连续的 g,一共多少种
// 这次是给出来的,上次是问让你染有多少种,那个是暴力 dp,观察一下后
// 化简。
// 这里是暴力枚举 a-b,后面的用一个数组存答案
// 从 n 种不同元素中取出 m 的元素 C(n+m-1,m)
// 当每一种类的个数小于 m 时比如是 a,就要减去这个种类中在上述组合中
// 的多余的为空个数。
// 考虑一个已经放了 a 个,然后有多少种,用容斥++--就可以了!
// $polya
// zoj2344
// 给定一张有 n*m 个格子的纸,每个格子有黑白两种颜色可以染。现在先将
// 纸按长边粘起来得到一个圆柱,再将纸按短边拈起来得到一个游泳圈。如
// 果两种染色方案卷起来后是一样的,那么它们同构。问不同构的的染色方
// 案。n,m <= 20
// 解题思路:Polya 的核心是找置换,本题有两类置换,一类是滚动,一类
// 是旋转。
// 我们将 n*m 的纸想成一个矩阵,那么第一行是 1,2,3...m,第二行是
// m+1...m*2,以此类推。
// 第一类滚动置换有 n*m 个。行和列是分开的,行的滚动有 n 个,列的滚动
// 有 m 个,它们相乘就是总置换数。每个置换的循环节个数我们就模拟下,
// O(n*m)找到循环节。
// 第二类旋转置换要分类讨论,当 n!= m 时,我们可以将整个矩阵旋转 180
// 度,这里是旋转而不是翻转,我写了 Rorate 函数将矩阵旋转九十度,
// Rorate 两次.当 n==m 时,不仅可以旋转 180 度,旋转 90 度和 270 度
// 都要考虑在内。
// 也就是说当 n!=m 时,置换有 2*n*m 个,当 n==m 时,置换有 4*n*m 个,
// 再套个高精度
// TAG: Polya
// 设 G 是 p 个对象的一个置换群,用 m 种颜色涂染 p 个对象,则不同染色方
// 案为
// 其中 G={g1,...gs} c(gi )为置换 gi 的循环节数(i=1...s)

// $Java 大数
// import java.math.*;
// import java.util.*;
// public class Main {
// private static BigInteger TWO = BigInteger.valueOf(2);
// private static boolean[] mark = new boolean[512];
// public static int cyclic(int[] p) {
// int ret = 0;
// Arrays.fill(mark, 0, p.length, false);
// for (int i = 0; i < p.length; ++i) {
// if (!mark[i]) {

// ++ret;
// for (int j = i; !mark[j]; j = p[j]) {
// mark[j] = true;
// }
// }
// }
// return ret;
// }
// public static void main(String[] args) {
// Scanner in = new Scanner(System.in);
// while (in.hasNextInt()) {
// int n = in.nextInt();
// int m = in.nextInt();
// int nm = n * m;
// int[][] p;
// if (n == m) {
// p = new int[4][nm];
// } else {
// p = new int[2][nm];
// }
// BigInteger num = BigInteger.ZERO;
// int den = 0;
// for (int i = 0; i < n; ++i) {
// for (int j = 0; j < m; ++j) {
// for (int k = 0; k < nm; ++k) {
// int x = ((k / m) + i) % n;
// int y = ((k % m) + j) % m;
// p[0][k] = x * m + y;
// p[1][k] = (n - 1 - x) * m + (m - 1 - y);
// // nm - 1 - p[0][k];
// if (n == m) {
// p[2][k] = y * m + (m - 1 - x);
// p[3][k] = (n - 1 - y) * m + x;
// }
// }
// for (int[] pp : p) {
// num = num.add(TWO.pow(cyclic(pp)));
// ++den;
// }
// }
// }
// System.err.println(num);
// System.err.println(den);
// num = num.divide(BigInteger.valueOf(den));
// System.out.println(num);
// }
// }
// }

// Java 排序
// import java.util.Arrays;
// class Node implements Comparable {
// int a;
// Node(int a) {
// this.a = a;
// }
// @Override
// public int compareTo(Object arg0) {
// int a = ((Node)arg0).a;
// if(this.a > a)
// return -1;
// else if(this.a < a)
// return 1;
// return 0;
// }

// @Override
// public String toString() {
// return "" + a;
// }
// }
// public class Main {
// public static void main(String[]args) {
// Node[]node=new Node[4];
// node[0]=new Node(9);
// node[1]=new Node(8);
// node[2]=new Node(100);
// node[3]=new Node(0);
// Arrays.sort(node);
// System.out.println(Arrays.toString(node));
// }
// }

// 平方和
// 1^2+2^2+3^2+......+n^2=n(n+1)(2n+1)/6
// 立方和
// 1^3+2^3+3^3+......+n^3=[n(n+1)/2]^2
// 假设 n 分解因式之后的形式是 n = a1 ^ b1 * a2 ^ b2 *...*an
// ^ bn
// 它的因子个数就是(b1+1)*(b2+1)*...*(bn+1)
// 因子和是(1+a1 ^ 1+..+a1 ^ b1)*(1 + a2 ^ 1 +..+a2 ^
// b2 )*..*(1 + an ^ 1 + ... +an ^bn)

// $FFT 高精度乘法
#define N 200005
#define pi acos(-1.0)  // PI 值
using namespace std;
struct complex {
  double r, i;
  complex(double real = 0.0, double image = 0.0) {
    r = real;
    i = image;
  }
  // 以下为三种虚数运算的定义
  complex operator+(const complex o) {
    return complex(r + o.r, i + o.i);
  }
  complex operator-(const complex o) {
    return complex(r - o.r, i - o.i);
  }
  complex operator*(const complex o) {
    return complex(r * o.r - i * o.i, r * o.i + i * o.r);
  }
} x1[N], x2[N];
char a[N / 2], b[N / 2];
int sum[N];                  // 结果存在 sum 里
void brc(complex *y, int l)  // 二进制平摊反转置换 O(logn)
{
  register int i, j, k;
  for (i = 1, j = l / 2; i < l - 1; i++) {
    if (i < j) swap(y[i], y[j]);  // 交换互为下标反转的元素
    // i<j 保证只交换一次
    k = l / 2;
    while (j >= k)  // 由最高位检索,遇 1 变 0,遇 0 变 1,跳出
    {
      j -= k;
      k /= 2;
    }
    if (j < k) j += k;
  }
}
void fft(complex *y, int l, double on)  // FFT O(nlogn)
// 其中 on==1 时为 DFT,on==-1 为 IDFT
{
  register int h, i, j, k;
  complex u, t;
  brc(y, l);                    // 调用反转置换
  for (h = 2; h <= l; h <<= 1)  // 控制层数
  {
    // 初始化单位复根
    complex wn(cos(on * 2 * pi / h), sin(on * 2 * pi / h));
    for (j = 0; j < l; j += h)  // 控制起始下标
    {
      complex w(1, 0);                 // 初始化螺旋因子
      for (k = j; k < j + h / 2; k++)  // 配对
      {
        u = y[k];
        t = w * y[k + h / 2];
        y[k] = u + t;
        y[k + h / 2] = u - t;
        w = w * wn;  // 更新螺旋因子
      }              // 据说上面的操作叫蝴蝶操作...
    }
  }
  if (on == -1)
    for (i = 0; i < l; i++) y[i].r /= l;  // IDFT
}
int main(void) {
  int l1, l2, l;
  register int i;
  while (scanf("%s%s", a, b) != EOF) {
    l1 = strlen(a);
    l2 = strlen(b);
    l = 1;
    while (l < l1 * 2 || l < l2 * 2) l <<= 1;  // 将次数界变成 2^n
    // 配合二分与反转置换
    for (i = 0; i < l1; i++)  // 倒置存入
    {
      x1[i].r = a[l1 - i - 1] - '0';
      x1[i].i = 0.0;
    }
    for (; i < l; i++) x1[i].r = x1[i].i = 0.0;
    // 将多余次数界初始化为 0
    for (i = 0; i < l2; i++) {
      x2[i].r = b[l2 - i - 1] - '0';
      x2[i].i = 0.0;
    }
    for (; i < l; i++) x2[i].r = x2[i].i = 0.0;
    fft(x1, l, 1);                                   // DFT(a)
    fft(x2, l, 1);                                   // DFT(b)
    for (i = 0; i < l; i++) x1[i] = x1[i] * x2[i];   // 点乘结果存入 a
    fft(x1, l, -1);                                  // IDFT(a*b)
    for (i = 0; i < l; i++) sum[i] = x1[i].r + 0.5;  // 四舍五入
    for (i = 0; i < l; i++)                          // 进位
    {
      sum[i + 1] += sum[i] / 10;
      sum[i] %= 10;
    }
    l = l1 + l2 - 1;
    while (sum[l] <= 0 && l > 0) l--;                // 检索最高位
    for (i = l; i >= 0; i--) putchar(sum[i] + '0');  // 倒序输出
    putchar('\n');
  }

  return 0;
}
// $高次幂取模

// $树上分治:
const int N = 1e3 + 111;
const int M = N * 3;
struct edge {
  int next, v;
} e[M];
int head[N], tot;
void init() {
  FILL(head, 0);
  tot = 1;
}
inline void add(int u, int v) {
  e[tot].v = v;
  e[tot].next = head[u];
  head[u] = tot++;
}
bool f[N];
int n, m;
//```````````````
int mx[N], size[N], mi, root;
void dfssize(int u, int fa) {
  size[u] = 1;
  mx[u] = 0;
  for (int i = head[u]; i; i = e[i].next) {
    int v = e[i].v;
    if (v != fa && !f[v]) {
      dfssize(v, u);
      size[u] += size[v];
      if (size[v] > mx[u])
        mx[u] = size[v];
    }
  }
}
void dfsroot(int r, int u, int fa) {
  if (size[r] - size[u] > mx[u])
    mx[u] = size[r] - size[u];
  if (mx[u] < mi)
    mi = mx[u], root = u;
  for (int i = head[u]; i; i = e[i].next) {
    int v = e[i].v;
    if (v != fa && !f[v])
      dfsroot(r, v, u);
  }
}
void getroot(int st) {
  mi = n;
  dfssize(st, 0);
  dfsroot(st, st, 0);
}
//````````
int ans[N];

vector<int> k;
void dfs2(int st, int now, int pa) {
  k.push_back(now);
  for (int i = head[st]; i; i = e[i].next) {
    int v = e[i].v;
    if (v != pa && !f[v]) {
      dfs2(v, now + 1, st);
    }
  }
}
int cao[N];
void cal(int st) {
  int i, j;
  FILL(cao, 0);
  cao[0] = 1;
  for (i = head[st]; i; i = e[i].next) {
    int v = e[i].v;
    if (!f[v]) {
      k.clear();
      dfs2(v, 1, st);
      rep(j, k.size()) {
        int tp = k[j];
        for (int t = tp; t < n; t++) {
          ans[t] += cao[t - tp];
        }
      }
      rep(j, k.size()) {
        cao[k[j]]++;
      }
    }
  }
}
void dfs(int st) {
  getroot(st);
  f[root] = 1;
  cal(root);
  for (int i = head[root]; i; i = e[i].next) {
    int v = e[i].v;
    if (!f[v]) {
      dfs(v);
    }
  }
}
void work() {
  int a, b, i;
  init();
  rep(i, n - 1) {
    scanf("%d%d", &a, &b);
    add(a, b);
    add(b, a);
  }
  FILL(f, 0);
  FILL(ans, 0);
  dfs(1);
  while (m--) {
    scanf("%d", &a);
    printf("%d\n", ans[a]);
  }
}

// $滚动 hash
// 给一副 1000*1000 的图
// 100 个恒定一样大小的测试图
// 问几个测试图是出现在图中的
// 先行 hash
// 再列 hash
// 这里使用 multiset 存 100 个然后找要快
int N, M, T, P, Q;
char s[1001][1001], a[1001][1001];
multiset<ull> st;
ull hash[1001][1001], tmp[1001][1001];
void computer_hash(char a[1001][1001], int n, int m) {
  const ull B1 = 9973;
  const ull B2 = 100000007;
  ull t1 = 1;
  int i, j;
  rep(j, Q) t1 *= B1;
  rep(i, n) {
    ull e = 0;
    rep(j, Q) e = e * B1 + a[i][j];
    for (j = 0; j + Q <= m; j++) {
      tmp[i][j] = e;
      e = e * B1 - t1 * a[i][j] + a[i][j + Q];
    }
  }
  ullt2 = 1;
  rep(i, P) t2 *= B2;
  rep(j, m) {
    ull e = 0;
    rep(i, P) e = e * B2 + tmp[i][j];
    for (i = 0; i + P <= n; i++) {
      hash[i][j] = e;
      e = e * B2 - t2 * tmp[i][j] + tmp[i + P][j];
    }
  }
}
void work() {
  st.clear();
  int i, j;
  rep(i, N) {
    scanf("%s", s[i]);
  }
  rep(i, T) {
    rep(j, P) {
      scanf("%s", a[j]);
    }
    computer_hash(a, P, Q);
    st.insert(hash[0][0]);
  }
  computer_hash(s, N, M);
  for (i = 0; i + P <= N; i++) {
    for (j = 0; j + Q <= M; j++) {
      st.erase(hash[i][j]);
    }
  }
  cout << T - st.size() << endl;
}
int main() {
  int cas = 1;
  while (cin >> N >> M >> T >> P >> Q && N) {
    printf("Case %d: ", cas++);
    work();
  }
  return 0;
}
// $$计算几何:
const double eps = 1e-7;
const double PI = acos(-1.0);
inline int sign(double x) {
  //中间过程大的话 eps 可适当乘 500 左右
  return (x > eps) - (x < -eps);
}
#define rep(i, n) for (i = 0; i < n; i++)
struct Point {
  double x, y;
  Point() {
  }
  Point(double _x, double _y) {
    x = _x;
    y = _y;
  }
  void in() {
    scanf("%lf%lf", &x, &y);
  }
  bool operator==(const Point B) const {
    return sign(x - B.x) == 0 && sign(y - B.y) == 0;
  }
  Point operator-(const Point B) const {
    return Point(x - B.x, y - B.y);
  }
  Point operator+(const Point B) const {
    return Point(x + B.x, y + B.y);
  }
  Point operator/(const double B) const {
    return Point(x / B, y / B);
  }
  Point operator*(const double B) const {
    return Point(x * B, y * B);
  }
};
typedef Point Vector;
bool operator<(const Point &a, const Point &b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}
int dcmp(double x) {
  if (fabs(x) < eps) {
    return 0;
  } else
    return x < 0 ? -1 : 1;
}
double Dot(Vector a, Vector b) {
  //计算点积
  return a.x * b.x + a.y * b.y;
}
inline double Dot(Point O, Point A, Point B) {
  return Dot(A - O, B - O);
}
inline double dist(Point A, Point B) {
  return sqrt(Dot(A, B, B));
}
double Cross(Vector a, Vector b) {
  //计算叉积
  return a.x * b.y - a.y * b.x;
}
double Cross(Point a, Point b, Point c) {
  return Cross(b - a, c - a);
}
double Length(Vector a) {
  //计算长度
  return sqrt(Dot(a, a));
}
double Angle(Vector A, Vector B) {
  //计算 AB 夹角
  return acos(Dot(A, B) / Length(A) / Length(B));
}
Vector Rotate(Vector a, double rad) {
  //向量逆时针旋转 rad 角
  //rad 弧度
  return Vector(a.x * cos(rad) - a.y * sin(rad),
                a.x * sin(rad) + a.y * cos(rad));
}
Vector Normal(Vector a) {
  //计算向量的单位法线,即左转 90 度,长度归一化。
  //用前确保 A 不是零向量
  double l = Length(a);
  return Vector(-a.y / l, a.x / l);
}
struct Segment {
  Point a, b;
  Segment() {
  }
  Segment(Point _a, Point _b) {
    a = _a;
    b = _b;
  }
};
struct Line {
  //有向直线,左边就是对应的半平面
  Vector v;
  Point p, p2;
  double ang;
  double ia, ib, ic;  //iax + iby + ic = 0
  Line() {}
  Line(Point pp, Vector vv) {
    p = pp;
    v = vv;
    ang = atan2(v.y, v.x);
    p2 = p + v;
  }
  Point point(double a) {
    return v + (p - v) * a;
  }
  bool operator<(const Line &L) const {
    return ang < L.ang;
  }
  void make()  //(y2-y1)x - (x2-x1)y + x2y1 - x1y2 = 0
  {
    ia = p2.y - p.y;
    ib = p.x - p2.x;
    ic = p2.x * p.y - p.x * p2.y;
    //if(sign(ia)<0)ia = -ia , ib = -ib , ic = -ic;
  }
};
int inSegment(Point pt, Segment seg) {
  //cout<<Cross(pt, seg.a, seg.b)<<endl;
  if (sign(Cross(pt, seg.a, seg.b)))
    return 0;
  int ans = sign(Dot(pt, seg.a, seg.b));
  return ans == 0 ? 2 : ans < 0;
}
//线段和线段相交情况 :
// 0:不相交, 1: 规范相交 , 2:交于端点, 3:有重合部分
// 快速排斥实验+ 相互跨立(是相互)
int aCross(Segment AB, Segment CD) {
  if (AB.b < AB.a) swap(AB.b, AB.a);
  if (CD.b < CD.a) swap(CD.b, CD.a);
  if (AB.a == CD.a && AB.b == CD.b) return 3;
  if (max(AB.a.x, AB.b.x) < min(CD.a.x, CD.b.x) || max(AB.a.y, AB.b.y) < min(CD.a.y, CD.b.y) || max(CD.a.x, CD.b.x) < min(AB.a.x, AB.b.x) || max(CD.a.y, CD.b.y) < min(AB.a.y, AB.b.y))
    return 0;
  int AB_CD = sign(Cross(AB.a, AB.b, CD.a)) * sign(Cross(AB.a, AB.b, CD.b));
  int CD_AB = sign(Cross(CD.a, CD.b, AB.a)) * sign(Cross(CD.a, CD.b, AB.b));
  //有重合部分
  if (AB_CD == 0 && CD_AB == 0 && (inSegment(AB.a, CD) == 1 || inSegment(AB.b, CD) == 1 || inSegment(CD.a, AB) == 1 || inSegment(CD.b, AB) == 1))
    return 3;
  if (AB_CD < 0)  //CD 跨立 AB
  {
    return CD_AB == 0 ? 2 : CD_AB < 0;
  } else if (AB_CD == 0)
    return CD_AB <= 0 ? 2 : 0;
  return 0;
}
bool in_line(Point a, Line d) {
  Point b = d.p, c = d.p2;
  return fabs((b.y - a.y) * (b.x - c.x) - (b.y - c.y) * (b.x - a.x)) < eps;
}
//线段和直线相交情况 :
// 0:不相交, 1: 规范相交 ,
// 2:不规范相交(交于端点或重合)
// 判断线段与直线相交只需判断线段是否跨立直线即可
int aCross(Segment seg, Line line) {
  //if(in_line (seg. a,line) &&

  // in_line(seg.b ,line)) return 0;// 如果重合返回
  int ans = sign(Cross(line.p, line.p2, seg.a)) * sign(Cross(line.p, line.p2, seg.b));
  return ans == 0 ? 2 : ans < 0;
}
//直线和直线相交情况 :
// 0:不相交( 平行), 1:规范相交 , 2:不规范相交(重合)
int aCross(Line AB, Line CD) {
  if (sign(Cross(Point(0, 0),
                 AB.p2 - AB.p, CD.p2 - CD.p)))
    return 1;
  return sign(Cross(AB.p, AB.p2, CD.p)) == 0 ? 2 : 0;
}
//计算交点前需先判断是否相交直线可以是平行于 y 轴的
Point intersect(Line LA, Line LB) {
  LA.make();
  LB.make();
  double a1 = LA.ia, b1 = LA.ib, c1 = LA.ic;
  double a2 = LB.ia, b2 = LB.ib, c2 = LB.ic;
  double x = (c1 * b2 - b1 * c2) / (a2 * b1 - a1 * b2);
  double y;
  if (sign(b1))
    y = -(a1 * x + c1) / b1;
  else
    y = -(a2 * x + c2) / b2;
  return Point(x, y);
}
Point intersect(Segment SA, Segment SB) {
  if (SA.a == SB.a || SA.a == SB.b)
    return SA.a;
  if (SA.b == SB.a || SA.b == SB.b)
    return SA.b;
  double AB_C = Cross(SA.a, SA.b, SB.a);
  double AB_D = Cross(SA.a, SA.b, SB.b);
  double x = (AB_D * SB.a.x - AB_C * SB.b.x) / (AB_D - AB_C);
  double y = (AB_D * SB.a.y - AB_C * SB.b.y) / (AB_D - AB_C);
  return Point(x, y);
}
Point intersect(Segment seg, Line line) {
  Line _line = Line(seg.a, seg.b);
  return intersect(_line, line);
}
Point Intersection(Point a1, Point a2,
                   Point b1, Point b2) {
  Point a = a2 - a1;
  Point b = b2 - b1;
  Point s = b1 - a1;
  return a1 + a * (Cross(b, s) / Cross(b, a));
}
Point GetLineProjection(Point p, Point a, Point b) {
  //求点在直线上投影
  //AB 直线:A+tv , 投影 Q 的参数为 t0
  //由 Dot(v,P-(A+t0*v))==0 得出
  Vector v = b - a;
  return a + v * (Dot(v, p - a) / Dot(v, v));
}
//这 2 个函数是实现半平米交的
Point GetLineIntersection(Line a, Line b) {
  //二直线相交,假设交点唯一存在
  Vector u = a.p - b.p;
  double t = Cross(b.v, u) / Cross(a.v, b.v);

  return a.p + a.v * t;
}
bool SegmentProperIntersection(Point a1, Point a2,
                               Point b1, Point b2) {
  //线段相交判定
  double c1 = Cross(a2 - a1, b1 - a1);
  double c2 = Cross(a2 - a1, b2 - a1);
  double c3 = Cross(b2 - b1, a1 - b1);
  double c4 = Cross(b2 - b1, a2 - b1);
  return dcmp(c1) * dcmp(c2) < 0 && dcmp(c3) * dcmp(c4) < 0;
}
double DistanceToLine(Point p, Point a, Point b) {
  Point v1 = b - a;
  Point v2 = p - a;
  return fabs(Cross(v1, v2)) / Length(v1);
}
double DistanceToSegment(Point p, Point a, Point b) {
  //计算点 P 到线段 AB 距离
  if (a == b) {
    return Length(p - a);
  }
  Vector v1 = b - a, v2 = p - a, v3 = p - b;
  if (dcmp(Dot(v1, v2)) < 0) {
    return Length(v2);
  } else if (dcmp(Dot(v1, v3)) > 0) {
    return Length(v3);
  } else {
    return fabs(Cross(v1, v2)) / Length(v1);
  }
}
bool OnSegment(Point p, Point a1, Point a2) {
  //判断一个点在线段上
  return dcmp(Cross(a1 - p, a2 - p)) == 0 && dcmp(Dot(a1 - p, a2 - p)) < 0;
}
double ConvexPolygonArea(Point *p, int n) {
  double area = 0;
  for (int i = 1; i < n - 1; i++) {
    area += Cross(p[i] - p[0], p[i + 1] - p[0]);
  }
  return area / 2;
}
//-------------------------------------------------
// $Circle
struct Circle {
  Point o;
  double R;
  Circle() {}
  Circle(Point cc, double rr) {
    o = cc;
    R = rr;
  }

  Point point(double a) {
    //通过圆心角来求坐标
    return Point(o.x + cos(a) * R, o.y + sin(a) * R);
  }
  //圆弧面积,a 圆心
  double s_y(Point a, Point b, Point c) {
    double d1 = dist(a, b), d2 = dist(a, c),
           d3 = dist(b, c);
    return fabs(
        R * R / 2 *
        acos((d1 * d1 + d2 * d2 - d3 * d3) / (2 * d1 * d2)));
  }
  void getCross(Point b, Point c, vector<Point> &k) {
    Point a = o;
    double A = (c.x - b.x) * (c.x - b.x) + (c.y - b.y) * (c.y - b.y),
           B = 2 * ((c.x - b.x) * (b.x - a.x) + (c.y - b.y) * (b.y - a.y)),
           C = a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y - 2 * (a.x * b.x + a.y * b.y) - R * R;
    double die = B * B - 4 * A * C + eps;
    if (die < eps)
      return;
    double u1, u2;
    u1 = (-B - sqrt(die)) / 2 / A;
    u2 = (-B + sqrt(die)) / 2 / A;
    Point s1(b.x + u1 * (c.x - b.x),
             b.y + u1 * (c.y - b.y));
    Point s2(b.x + u2 * (c.x - b.x),
             b.y + u2 * (c.y - b.y));
    if (s1.y < o.y + eps && inSegment(s1, Segment(b, c)))
      k.push_back(s1);
    if (s2.y < o.y + eps && inSegment(s2, Segment(b, c)))
      k.push_back(s2);
  }
  // 三角形 abc 和 a 点 R 半径的圆的有向面积交
  double s_Cross_y(Point a, Point b, Point c) {
    double f = sign(Cross(a, b, c));
    double d1 = dist(a, b), d2 = dist(a, c),
           d3 = dist(b, c),
           d4 = fabs(Cross(a, b, c)) / d3;
    if (sign(d1) == 0 || sign(d2) == 0)
      return 0;
    bool ff = 0;
    if (d1 <= R && d2 <= R)
      return f * fabs(Cross(a, b, c)) / 2;
    if (d1 >= R && d2 >= R) {
      if (d4 < R)
        ff = 1;  //标记
      else
        return f * s_y(a, b, c);
    }
    if (d1 > d2) {
      swap(d1, d2);
      swap(c, b);
    }
    //计算 bc 和圆交点
    double A = (c.x - b.x) * (c.x - b.x)

               + (c.y - b.y) * (c.y - b.y),
           B = 2 * ((c.x - b.x) * (b.x - a.x) + (c.y - b.y) * (b.y - a.y)),
           C = a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y - 2 * (a.x * b.x + a.y * b.y) - R * R;
    double die = B * B - 4 * A * C + eps;
    double u1, u2;
    u1 = (-B - sqrt(die)) / 2 / A;
    u2 = (-B + sqrt(die)) / 2 / A;
    Point s1(b.x + u1 * (c.x - b.x),
             b.y + u1 * (c.y - b.y));
    Point s2(b.x + u2 * (c.x - b.x),
             b.y + u2 * (c.y - b.y));
    if (ff) {
      if (sign(Cross(a, s1, b)) * sign(Cross(a, s1, c)) >= 0 &&
          sign(Cross(a, s2, b)) * sign(Cross(a, s2, c)) >= 0)  //特判
        return f * s_y(a, b, c);
      double t1 = s_y(a, s1, b), t2 = s_y(a, s1, c);
      if (t1 > t2)
        return f * (t2 + fabs(Cross(a, s1, s2) / 2) + s_y(a, s2, b));
      else
        return f * (t1 + fabs(Cross(a, s1, s2) / 2) + s_y(a, s2, c));
    }
    double a1 = a.x, a2 = b.x, b1 = a.y, b2 = b.y;
    if (a1 > a2)
      swap(a1, a2);
    if (b1 > b2)
      swap(b1, b2);
    if (s1.x >= a1 && s1.x <= a2 && s1.y >= b1 && s1.y <= b2)
      return f * (fabs(Cross(a, b, s1)) / 2 + s_y(a, s1,
                                                  c));
    else
      return f * (fabs(Cross(a, b, s2)) / 2 + s_y(a, s2,
                                                  c));
  }
};
// $Polygon
struct Polygon {
  int n;
  Point p[50000];
  Polygon() {
    n = 0;
  }
  Polygon(int nn, Point *pp) {
    n = nn;
    for (int q = 0; q < n; q++) {
      p[q] = pp[q];
    }
  }
  Point operator[](int a) const {
    return p[a];
  }
};
bool OnLeft(Line l, Point pp) {
  //点P在有向直线的左边,线上不算
  return Cross(l.v, pp - l.p) > 0;
}
double angle(Vector v) {
  //计算向量的极角
  return atan2(v.y, v.x);
}
int getTangents(Point p, Circle c, Vector *v) {
  //过点 p 到圆 C 的切线,
  //v[i]是第 i 条切线的向量,返回切线条数
  Vector u = c.o - p;
  double dist = Length(u);
  if (dist < c.R)
    return 0;
  else if (dcmp(dist - c.R) == 0) {
    v[0] = Rotate(u, PI / 2);
    return 1;
  } else {
    double ang = asin(c.R / dist);
    v[0] = Rotate(u, -ang);
    v[1] = Rotate(u, ang);
    return 2;
  }
}
int getTangents(Circle a, Circle b, Point *aa, Point *bb) {
  //返回切线条数,-1 表示无穷条切线
  //a[i]和 b[i]分别表示第 i 条切线在圆 A 和圆 B 上的切点
  int cnt = 0;
  if (a.R < b.R) {
    swap(a, b);
    swap(aa, bb);
  }
  double d2 = (a.o.x - b.o.x) * (a.o.x - b.o.x) + (a.o.y - b.o.y) * (a.o.y - b.o.y);
  double rdiff = a.R - b.R;
  double rsum = a.R + b.R;
  if (d2 < rdiff * rdiff) {
    return 0;  //内含
  }
  double base = atan2(b.o.y - a.o.y, b.o.x - a.o.x);
  if (d2 == 0 && a.R == b.R)
    return -1;  //无穷多条切线
  if (d2 == rdiff * rdiff) {
    aa[cnt] = a.point(base);
    bb[cnt] = b.point(base);
    cnt++;

    return 1;  //内切,1 条
  }
  //有外公切线
  double ang = acos((a.R - b.R) / sqrt(d2));
  aa[cnt] = a.point(base + ang);
  bb[cnt] = b.point(base + ang);
  cnt++;
  aa[cnt] = a.point(base - ang);
  bb[cnt] = b.point(base - ang);
  cnt++;
  if (d2 == rsum * rsum) {
    aa[cnt] = a.point(base);
    bb[cnt] = b.point(base + PI);
    cnt++;  //一条内公切线
  } else if (d2 > rsum * rsum) {
    double ang = acos((a.R + b.R) / sqrt(d2));
    aa[cnt] = a.point(base + ang);
    bb[cnt] = b.point(PI + base + ang);
    cnt++;
    aa[cnt] = a.point(base - ang);
    bb[cnt] = b.point(PI + base - ang);
    cnt++;  //两条内公切线
  }
  return cnt;
}
Circle CircumscribedCircle(Point p1, Point p2, Point p3) {
  //三角形外接圆
  double bx = p2.x - p1.x;
  double by = p2.y - p1.y;
  double cx = p3.x - p1.x;
  double cy = p3.y - p1.y;
  double d = 2 * (bx * cy - cx * by);
  double ccx = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d + p1.x;
  double ccy = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d + p1.y;
  Point p = Point(ccx, ccy);
  return Circle(p, Length(p1 - p));
}
Circle InscribedCircle(Point p1, Point p2, Point p3) {
  //三角形内切圆
  double a = Length(p2 - p3);
  double b = Length(p3 - p1);
  double c = Length(p1 - p2);
  Point p = (p1 * a + p1 * b + p3 * c) / (a + b + c);
  return Circle(p, DistanceToLine(p, p1, p2));
}
double torad(double deg) {
  //角度转弧度
  return deg / 180 * PI;
}
void get_coord(double r, double lat, double lng,
               double &x, double &y, double &z) {
  //经纬度(角度)转化为空间坐标
  lat = torad(lat);
  lng = torad(lng);

  x = r * cos(lat) * cos(lng);
  y = r * cos(lat) * sin(lng);
  z = r * sin(lat);
}
//求球面上两点最短距离,走大圆弧
//半径 r,弦长 d,圆心角 2arcsin(d/2r)
int isPointInPolygon(Point p, Polygon po) {
  //假想有一条向右的射线,
  //统计多边形穿过这条射线正反多少次,把这个数记为绕数 S
  //逆时穿过加 1,顺时针穿过减 1
  int wn = 0;
  int n = po.n;
  for (int i = 0; i < n; i++) {
    if (OnSegment(p, po[i], po[(i + 1) % n]))
      return -1;  //在边界上
    int k = dcmp(Cross(po[(i + 1) % n] - po[i], p - po[i]));
    int d1 = dcmp(po[i].y - p.y);
    int d2 = dcmp(po[(i + 1) % n].y - p.y);
    if (k > 0 && d1 <= 0 && d2 > 0)
      wn++;
    if (k < 0 && d2 <= 0 && d1 > 0)
      wn--;
  }
  if (wn != 0)
    return 1;  //内部
  return 0;    //外部
}
//点在凸多边形内判定更简单,只要判断是否在所有边的左边
//各顶点按逆时针顺序
inline bool cmpp(Point a, Point b) {
  return a.y < b.y || (a.y == b.y && a.x < b.x);
}
int grahamConvex(Point points[], int n, Point done[]) {
  int top, len;
  sort(points, points + n, cmpp);
  n = unique(points, points + n) - points;
  for (int i = 0; i <= 1; i++)
    done[i] = points[i];
  top = 1;
  for (int i = 2; i < n; i++)  //建立右链
  {
    while (top && dcmp(Cross(points[i], done[top - 1],
                             done[top])) >= 0)
      top--;
    done[++top] = points[i];
  }
  len = top;
  for (int i = n - 2; i >= 0; i--)  //建立左链
  {
    while (top != len && dcmp(Cross(points[i], done[top - 1],
                                    done[top])) >= 0)  // 防止右链中的数据跳出
      top--;
    done[++top] = points[i];
  }
  return top;
}
void RC(Point p[], int top)  //旋转卡壳
{
  int q;
  double ans = 0;
  q = 1;
  ans = Length(p[0] - p[1]);

  for (int i = 0; i != top; i++) {
    while (fabs(Cross(p[(i + 1) % top],
                      p[i % top], p[(q + 1) % top])) > fabs(Cross(p[(i + 1) % top],
                                                                  p[i % top], p[q % top]))) {
      q = (q + 1) % top;
    }
    ans = max(ans, max(Length(p[(i + 1) % top] - p[q]),
                       Length(p[i % top] - p[q])));
  }
  printf("%.lf\n", ans * ans);
}
//结果面积为 0 输出不了点,要把 onleft 改成>=
int HalfplaneIntersection(Line *l, int n, Point *poly) {
  sort(l, l + n);
  int first, last;
  Point *p = new Point[n];
  Line *q = new Line[n];
  q[first = last = 0] = l[0];
  for (int i = 1; i < n; i++) {
    while (first < last && !OnLeft(l[i], p[last - 1]))
      last--;
    while (first < last && !OnLeft(l[i], p[first]))
      first++;
    q[++last] = l[i];
    if (fabs(Cross(q[last].v, q[last - 1].v)) < eps) {
      last--;
      if (OnLeft(q[last], l[i].p))
        q[last] = l[i];
    }
    if (first < last)
      p[last - 1] = GetLineIntersection(q[last - 1], q[last]);
  }
  while (first < last && !OnLeft(q[first], p[last - 1]))
    last--;
  if (last - first <= 1)
    return 0;
  p[last] = GetLineIntersection(q[last], q[first]);
  int m = 0;
  for (int i = first; i <= last; i++)
    poly[m++] = p[i];
  return m;
}
//要确保 as 是直线,切出的是 as 左手边的
Polygon cut(Polygon &po, Point a, Point s) {
  Point x, c;
  Polygon ans;
  int sum = 0;
  for (int q = 0; q < po.n; q++) {
    x = po[q];
    c = po[(q + 1) % po.n];
    //OnSegment(x,a,s)
    if (OnLeft(Line(a, s - a), x) || in_line(x, Line(a, s - a))) {
      ans.p[sum] = x;
      sum++;
    }
    //因为这里要判断是否相交!,x 和 c 不能取

    if (aCross(Segment(x, c), Line(a, s - a)) == 1) {  //SegmentProperIntersection(x,c,a,s)
      ans.p[sum] = GetLineIntersection(Line(x, c - x),
                                       Line(a, s - a));
      sum++;
    }
  }
  ans.n = sum;
  return ans;
}

// $模拟退火
// 求一个平面的 ax+by+cz = 0 的 a,b,c,要这个平面到每个点的距离
// 的平方的和最小
int sgn(double x) {
  return (x < -eps) ? -1 : (x > eps);
}
struct PT {
  double x, y, z;
  PT() {
  }
  PT(double x, double y, double z) : x(x), y(y), z(z) {
  }
};
PT pt[105];
double X2, Y2, Z2, XY, XZ, YZ;
double work(PT pp) {
  double a = pp.x, a2 = a * a;
  double b = pp.y, b2 = b * b;
  double c = pp.z, c2 = c * c;
  double now = a2 * X2 + b2 * Y2 + c2 * Z2 +
               a * b * XY * 2 + a * c * XZ * 2 + b * c * YZ * 2;
  return now / (a2 + b2 + c2);
}
int main() {
  int cas, T;
  srand(330598937);
  for (cas = scanf("%d", &T); cas <= T; cas++) {
    int n;
    scanf("%d", &n);
    X2 = Y2 = Z2 = XY = XZ = YZ = 0;
    for (int i = 0; i < n; i++) {
      double x, y, z;
      scanf("% lf %lf %lf", &x, &y, &z);
      X2 += x * x;
      Y2 += y * y;
      Z2 += z * z;
      XY += x * y;
      XZ += x * z;
      YZ += y * z;
    }
    PT bt(1, 1, 1);
    double ans = work(bt);
    for (double step = 1e6; step > eps; step /= 2) {
      for (int i = 1; i <= 500; i++) {
        double th1 = rand();

        double th2 = rand();
        PT p(bt.x + step * sin(th1), bt.y + step * cos(th1),
             bt.z + step * cos(th2));
        double tmp = work(p);
        if (tmp < ans)
          ans = tmp, bt = p;
      }
    }
    printf("Case %d: ", cas);
    printf("%.8f %.8f %.8f\n", bt.x, bt.y, bt.z);
  }
  return 0;
}
// 求到很多点的距离的平方的和最小的点,反正就是求导为 0,如果在直线
// 上就把 y 代成 x 形式后求导。

// $欧拉定理:
// 平面的顶点数,边数,面数分别为 V,E,F,则 V + F - E = 2
// n 维最大曼哈顿距离
// 一维就是 Max (x) - Min(x)
// 就是对于 二维的 x - y 和 x + y 做两个集合。 答案肯定会在
// Max( x - y) - Min( x - y) 或者 是 Max(x + y) - Min(x
// + y)
// 而三维就是 Max(x + y + z) - Min(x + y + z) 或者是 Max(x +
// y - z) - Min(x + y - z) 略。
// 就是 a +- b +- c +- d +- e
// 2 ^ (k - 1)个集合里。

// $三维几何
int dcmp(double x) {
  if (fabs(x) < eps)
    return 0;
  else
    return x < 0 ? -1 : 1;
}
struct Point {
  double x, y;
  Point(double _x = 0, double _y = 0) : x(_x), y(_y) {}
};
struct Point3 {
  double x, y, z;
  Point3(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z) {}
  void in() {
    scanf("%lf%lf%lf", &x, &y, &z);
  }
  void show() {
    cout << "cao" << x << ' ' << y << ' ' << z << endl;
  }
};
typedef Point Vector;
typedef Point3 Vector3;
Vector operator+(Vector A, Vector B) {
  return Vector(A.x + B.x, A.y + B.y);
}
Vector operator-(Vector A, Vector B) {
  return Vector(A.x - B.x, A.y - B.y);
}
Vector operator*(Vector A, double p) {
  return Vector(A.x * p, A.y * p);
}
Vector operator/(Vector A, double p) {
  return Vector(A.x / p, A.y / p);
}
bool operator<(const Point &a, const Point &b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}
bool operator==(const Point &a, const Point &b) {
  return dcmp(a.x - b.x) == 0 && dcmp(a.x - b.y) == 0;
}
Vector3 operator+(Vector3 A, Vector3 B) {
  return Vector3(A.x + B.x, A.y + B.y, A.z + B.z);
}
Vector3 operator-(Vector3 A, Vector3 B) {
  return Vector3(A.x - B.x, A.y - B.y, A.z - B.z);
}
Vector3 operator*(Vector3 A, double p) {
  return Vector3(A.x * p, A.y * p, A.z * p);
}
Vector3 operator/(Vector3 A, double p) {
  return Vector3(A.x / p, A.y / p, A.z / p);
}
bool operator==(const Point3 &a, const Point3 &b) {
  return dcmp(a.x - b.x) == 0 &&
         dcmp(a.x - b.y) == 0 && dcmp(a.z - b.z) == 0;
}
Vector3 Cross(Vector3 A, Vector3 B) {
  return Vector3(A.y * B.z - A.z * B.y,
                 A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x);
}
double Dot(Vector3 A, Vector3 B) {
  return A.x * B.x + A.y * B.y + A.z * B.z;
}
double Length(Vector3 A) {
  return sqrt(Dot(A, A));
}
double Area2(Point3 A, Point3 B, Point3 C) {
  return Length(Cross(B - A, C - A));
}
double Angle(Vector3 A, Vector3 B) {
  return acos(Dot(A, B) / Length(A) / Length(B));
}
//p 到平面距离
double distoplane(Point3 p, Point3 p0, Vector3 n) {
  return fabs(Dot(p - p0, n) / Length(n));
}
//p 在平面的投影
Point3 getplaneprojection(Point3 p, Point3 p0,
                          Vector3 n) {
  double d = Dot(p - p0, n) / Length(n);
  return p - n * d / Length(n);
}
//直线与平面交点,没有考虑平行包含
Point3 lineplaneintersection(Point3 p1, Point3 p2,
                             Point3 p0, Vector3 n) {
  Vector3 v = p2 - p1;
  double t = (Dot(n, p0 - p1) / Dot(n, p2 - p1));
  return p1 + v * t;
}
double distoline(Point3 P, Point3 A, Point3 B) {
  Vector3 v1 = B - A, v2 = P - A;
  return Length(Cross(v1, v2)) / Length(v1);
}
double distoseg(Point3 P, Point3 A, Point3 B) {
  if (A == B) return Length(P - A);

  Vector3 v1 = B - A, v2 = P - A, v3 = P - B;
  if (dcmp(Dot(v1, v2) < 0))
    return Length(v2);
  else if (dcmp(Dot(v1, v3) > 0))
    return Length(v3);
  else
    return Length(Cross(v1, v2)) / Length(v1);
}
//直线与直线最近距离
double dislinetoline(Point3 a, Point3 b, Point3 c, Point3 d) {
  Vector3 t = Cross(a - b, c - d);
  if (Length(t) < eps) {
    return distoline(a, c, d);
  }
  return distoplane(a, c, t);
}
double dislinetoline2(Point3 a, Point3 b, Point3 c, Point3 d) {
  Vector3 t = Cross(a - b, c - d);
  if (Length(t) < eps) {
    return distoline(a, c, d);
  }
  Point3 x = getplaneprojection(a, c, t);
  return Length(a - x);
}
// $三维凸包
#define EPS 1e-5
#define N 505
#define PI acos(-1)
struct _3DCH {
  struct pt {
    double x, y, z;
    pt() {
    }
    pt(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {
    }
    pt operator-(const pt p1) {
      return pt(x - p1.x, y - p1.y, z - p1.z);
    }
    pt operator*(pt p) {
      return pt(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
    }  //叉乘
    double operator^(pt p) {
      return x * p.x + y * p.y + z * p.z;
    }  //点乘
  };
  struct fac {
    int a, b, c;  //表示凸包一个面上三个点的编号
    bool ok;      //表示该面是否属于最终凸包中的面
  };
  int n;
  pt P[N];
  int cnt;       //凸包表面的三角形数
  fac F[N * 8];  //凸包表面的三角形
  int to[N][N];
  double vlen(pt a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
  }  //向量长度
  double area(pt a, pt b, pt c) {
    return vlen((b - a) * (c - a));
  }  //三角形面积*2
  double volume(pt a, pt b, pt c, pt d) {
    return (b - a) * (c - a) ^ (d - a);

  }                           //四面体有向体积*6
  double ptof(pt &p, fac &f)  //正:点在面同向
  {
    pt m = P[f.b] - P[f.a],
       n = P[f.c] - P[f.a], t = p - P[f.a];
    return (m * n) ^ t;
  }
  void deal(int p, int a, int b) {
    int f = to[a][b];
    fac add;
    if (F[f].ok)
      if (ptof(P[p], F[f]) > EPS)
        dfs(p, f);
      else {
        add.a = b, add.b = a, add.c = p, add.ok = 1;
        to[p][b] = to[a][p] = to[b][a] = cnt;
        F[cnt++] = add;
      }
  }
  void dfs(int p, int cur) {
    F[cur].ok = 0;
    deal(p, F[cur].b, F[cur].a);
    deal(p, F[cur].c, F[cur].b);
    deal(p, F[cur].a, F[cur].c);
  }
  bool same(int s, int t) {
    pt &a = P[F[s].a], &b = P[F[s].b], &c = P[F[s].c];
    return fabs(volume(a, b, c, P[F[t].a])) < EPS && fabs(volume(a, b, c, P[F[t].b])) < EPS && fabs(volume(a, b, c, P[F[t].c])) < EPS;
  }
  void construct()  //构建三维凸包
  {
    cnt = 0;
    if (n < 4)
      return;
    /* 此段是为了保证前四个点不共面,若已保证,可去掉 */
    bool sb = 1;  //使前两点不共点
    for (int i = 1; i < n; i++)
      if (vlen(P[0] - P[i]) > EPS) {
        swap(P[1], P[i]);
        sb = 0;
        break;
      }
    if (sb)
      return;
    sb = 1;  //使前三点不共线
    for (int i = 2; i < n; i++)
      if (vlen((P[0] - P[1]) * (P[1] - P[i])) > EPS) {
        swap(P[2], P[i]);
        sb = 0;
        break;
      }
    if (sb)
      return;
    sb = 1;  //使前四点不共面
    for (int i = 3; i < n; i++)
      if (fabs((P[0] - P[1]) *
                   (P[1] - P[2]) ^
               (P[0] - P[i])) > EPS) {
        swap(P[3], P[i]);
        sb = 0;
        break;
      }
    if (sb)
      return;
    fac add;

    for (int i = 0; i < 4; i++) {
      add.a = (i + 1) % 4, add.b = (i + 2) % 4;
      add.c = (i + 3) % 4, add.ok = 1;
      if (ptof(P[i], add) > 0)
        swap(add.b, add.c);
      to[add.a][add.b] = to[add.b][add.c] = to[add.c][add.a] = cnt;
      F[cnt++] = add;
    }
    for (int i = 4; i < n; i++)
      for (int j = 0; j < cnt; j++)
        if (F[j].ok && ptof(P[i], F[j]) > EPS) {
          dfs(i, j);
          break;
        }
    int tmp = cnt;
    cnt = 0;
    for (int i = 0; i < tmp; i++)
      if (F[i].ok)
        F[cnt++] = F[i];
  }
  double volume()  //体积
  {
    pt O(0, 0, 0);
    double ret = 0.0;
    for (int i = 0; i < cnt; i++)
      ret += volume(O, P[F[i].a], P[F[i].b], P[F[i].c]);
    return fabs(ret / 6.0);
  }
  pt Fc[N * 8];
  double V[N * 8];
  double ptoface(pt p, int i) {  //点到面的距离
    return fabs(
        volume(P[F[i].a], P[F[i].b], P[F[i].c], p) / vlen(
                                                         (P[F[i].b] - P[F[i].a]) * (P[F[i].c] - P[F[i].a])));
  }
  void get_panel(double &a, double &b, double &c,
                 double &d, pt p1, pt p2, pt p3) {
    a = (p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y);
    b = (p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z);
    c = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    d = -(a * p1.x + b * p1.y + c * p1.z);
  }
  inline double dist(pt p1, pt p2) {
    return sqrt(
        (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
  }
} hull;
inline int sign(double d) {
  if (d > EPS)
    return 1;
  if (d < -EPS)
    return -1;
  return 0;
}
struct point {
  double x, y, z;
  point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {
  }
  void read() {
    scanf("%lf%lf%lf", &x, &y, &z);
  }
  point operator-(point tp) {
    return point(x - tp.x, y - tp.y, z - tp.z);
  }
  point operator+(point tp) {
    return point(x + tp.x, y + tp.y, z + tp.z);
  }
  point operator/(double len) {
    return point(x / len, y / len, z / len);
  }
  point operator*(double len) {
    return point(x * len, y * len, z * len);
  }
  point operator*(point tp) {
    return point(x * tp.x, y * tp.y, z * tp.z);
  }
  double getLen() {  //计算向量长度
    return sqrt(x * x + y * y + z * z);
  }
} con[N];
point ps[N], org;
double a, b, c, d;
int n;
//求点 tp 在直线(st,ed)上的垂足
inline point get_point(point st, point ed, point tp) {
  double a, b, t;
  a = (ed.x - st.x) * (ed.x - st.x) + (ed.y - st.y) * (ed.y - st.y) + (ed.z - st.z) * (ed.z - st.z);
  b = (tp.x - st.x) * (ed.x - st.x) + (tp.y - st.y) * (ed.y - st.y) + (tp.z - st.z) * (ed.z - st.z);
  t = b / a;
  point ans;
  ans.x = ed.x * t + (1 - t) * st.x;
  ans.y = ed.y * t + (1 - t) * st.y;
  ans.z = ed.z * t + (1 - t) * st.z;
  return ans;
}
//向量(0,0)->p1 X (0,0)->p2 得到的向量的终点坐标
point xmul(point p1, point p2) {
  point ans;
  ans.x = p1.y * p2.z - p1.z * p2.y;
  ans.y = p1.z * p2.x - p1.x * p2.z;
  ans.z = p1.x * p2.y - p1.y * p2.x;
  return ans;
}
inline double dist(point p1, point p2) {
  return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}
//求 tp 沿(st,ed)右旋转 ang 角度之后的点
point rotate(point st, point ed, point tp,
             double ang) {
  point root = get_point(st, ed, tp), e, r;
  point ans;
  e = (ed - st) / dist(st, ed);
  r = tp - root;
  e = xmul(e, r);
  ans = r * cos(ang) + e * sin(ang) + root;
  return ans;
}

double inter_pro(point st1, point ed1,
                 point st2, point ed2) {
  return (ed1.x - st1.x) * (ed2.y - st2.y) - (ed1.y - st1.y) * (ed2.x - st2.x);
}
bool cmp(point d1, point d2) {
  return d1.y < d2.y || (d1.y == d2.y && d1.x < d2.x);
}
// $$其它:
// $Dancing Links
// 精确覆盖
void remove(int &c) {
  L[R[c]] = L[c];
  R[L[c]] = R[c];
  for (int i = D[c]; i != c; i = D[i]) {
    for (int j = R[i]; j != i; j = R[j]) {
      U[D[j]] = U[j];
      D[U[j]] = D[j];
      --S[Col[j]];
    }
  }
}
void resume(int &c) {
  for (int i = U[c]; i != c; i = U[i]) {
    for (int j = L[i]; j != i; j = L[j]) {
      ++S[Col[j]];
      U[D[j]] = j;
      D[U[j]] = j;
    }
  }
  L[R[c]] = c;
  R[L[c]] = c;
}
bool dfs() {
  if (R[0] == 0) {
    return true;
  }
  int i, j;
  int idx, minnum = 999999;
  for (i = R[0]; i != 0; i = R[i]) {
    if (S[i] < minnum) {
      minnum = S[i];
      idx = i;
    }
  }
  remove(idx);
  for (i = D[idx]; i != idx; i = D[i]) {
    ans[deep++] = Row[i];
    for (j = R[i]; j != i; j = R[j]) {
      remove(Col[j]);
    }
    if (dfs()) {
      return true;
    }
    deep--;
    for (j = L[i]; j != i; j = L[j]) {
      resume(Col[j]);
    }
  }
  resume(idx);
  return false;
}
// 重复覆盖

void remove(int &c) {
  for (int i = D[c]; i != c; i = D[i]) {
    L[R[i]] = L[i];
    R[L[i]] = R[i];
  }
}
void resume(int &c) {
  for (int i = U[c]; i != c; i = U[i]) {
    L[R[i]] = i;
    R[L[i]] = i;
  }
}
int h() {
  bool hash[51];
  memset(hash, false, sizeof(hash));
  int ret = 0;
  for (int c = R[0]; c != 0; c = R[c]) {
    if (!hash[c]) {
      ret++;
      hash[c] = true;
      for (int i = D[c]; i != c; i = D[i]) {
        for (int j = R[i]; j != i; j = R[j]) {
          hash[Col[j]] = true;
        }
      }
    }
  }
  return ret;
}
bool dfs(int deep, int lim) {
  if (deep + h() > lim) {
    return false;
  }
  if (R[0] == 0) {
    return true;
  }
  int idx, i, j, minnum = 99999;
  for (i = R[0]; i != 0; i = R[i]) {
    if (S[i] < minnum) {
      minnum = S[i];
      idx = i;
    }
  }
  for (i = D[idx]; i != idx; i = D[i]) {
    remove(i);
    for (j = R[i]; j != i; j = R[j]) {
      remove(j);
    }
    if (dfs(deep + 1, lim)) {
      return true;
    }
    for (j = L[i]; j != i; j = L[j]) {
      resume(j);
    }
    resume(i);
  }
  return false;
}
// $表达式计算
typedef pair<int, int> result;
#define value first
#define p second
struct Equation {
  result equation(const string &s, int p) {
    result r = factor(s, p);
    while (s[r.p] == '+' || s[r.p] == '-') {
      result r_ = factor(s, r.p + 1);
      if (s[r.p] == '+') r.value += r_.value;
      if (s[r.p] == '-') r.value -= r_.value;
      r.p = r_.p;
    }
    return r;
  }
  result factor(const string &s, int p) {
    result r = term(s, p);
    while (s[r.p] == '*' || s[r.p] == '/') {
      result r_ = term(s, r.p + 1);
      if (s[r.p] == '*') r.value *= r_.value;
      if (s[r.p] == '/') r.value /= r_.value;
      r.p = r_.p;
    }
    return r;
  }
  result term(const string &s, int p) {
    if (s[p] == '(') {
      result r = equation(s, p + 1);
      r.p += 1;  // skip ')'
      return r;
    } else {
      int value = 0;
      while (isdigit(s[p])) {
        value = value * 10 + (s[p++] - '0');
      }
      return result(value, p);
    }
  }
} eq;
// $日期函数
//日期函数
int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
struct Date {
  int year, month, day;
};
//判闰年
inline int leap(int year) {
  return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}
//返回指定日期是星期几
int weekday(Date a) {
  int tm = a.month >= 3 ? (a.month - 2) : (a.month + 10);
  int ty = a.month >= 3 ? a.year : (a.year - 1);
  return (ty + ty / 4 - ty / 100 + ty / 400 + (int)(2.6 * tm - 0.2) + a.day) % 7;
}
//日期转天数偏移
int date2int(Date a) {
  int ret = a.year * 365 + (a.year - 1) / 4 - (a.year - 1) / 100 + (a.year - 1) / 400;
  days[1] += leap(a.year);
  for (int i = 0; i < a.month - 1; ret += days[i++])
    ;
  days[1] = 28;
  return ret + a.day;
}
//天数偏移转日期
Date int2date(int a) {
  Date ret;
  ret.year = a / 146097 * 400;
  for (a %= 146097; a >= 365 + leap(ret.year);
       a -= 365 + leap(ret.year), ret.year++)
    ;

  days[1] += leap(ret.year);
  for (ret.month = 1; a >= days[ret.month - 1];
       a -= days[ret.month - 1], ret.month++)
    ;
  days[1] = 28;
  ret.day = a + 1;
  return ret;
}
// $离线求动态区间完全被包含次数
// 2011 年吉隆坡赛区 Regional B 题 离线求动态区间完全被包含次数
// 方法归并 + 树状数组 离线搞二维矩形查询
// 每次+ a b 插入线段 a b 或? a b 询问 a b 被完全包含次数
// -10^9<=a,b<=10^9, query<=5*10^5
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
using namespace std;
const int MAXN = 512345;
struct data {
  int type, id;  //q01s23
  int qid, pid;
  int x;
  friend bool operator<(const data &a, const data &b) {
    return a.x == b.x ? a.type < b.type : a.x < b.x;
  }
} d[2 * MAXN];
int who[MAXN], who2[MAXN], ans[MAXN],
    out[MAXN], seq[MAXN];
const int MAXQ = 500001;
//调用上面的二维离线求和
int all[MAXN * 2], IDX;
int main() {
  int n;
  while (scanf("%d", &n) != EOF) {
    for (int i = 0; i < 2 * n; i++) d[i].id = i + 1;
    IDX = 0;
    int pnum = 0, qnum = 0;
    for (int i = 0; i < n; i++) {
      getchar();
      char ch = getchar();
      scanf("%d %d", &d[2 * i].x, &d[2 * i + 1].x);
      all[IDX++] = d[2 * i].x;
      all[IDX++] = d[2 * i + 1].x;
      if (ch == '+') {
        d[2 * i].type = 2;
        d[2 * i + 1].type = 3;
        d[2 * i].pid = pnum++;
        d[2 * i + 1].pid = pnum++;
      } else {
        d[2 * i].type = 0;
        d[2 * i + 1].type = 1;
        d[2 * i].qid = qnum++;
        d[2 * i + 1].qid = qnum++;
      }
    }
    sort(all, all + IDX);
    IDX = unique(all, all + IDX) - all;
    for (int i = 0; i < 2 * n; i++) {
      d[i].x = lower_bound(all, all + IDX, d[i].x) - all + 1;
    }
    sort(d, d + 2 * n);

    int num = 0;
    idx = qid = 0;
    for (int i = 0; i < 2 * n; i++) {
      if (d[i].type == 0) {
        who[d[i].qid] = d[i].x;
      } else if (d[i].type == 1) {
        int x1 = 1, y1 = 1,
            x2 = who[d[i].qid], y2 = d[i].id - 1;
        //printf("Ask %d %d %d %d\n", x1, y1, x2, y2);
        q[idx].x = x1 - 1, q[idx].y = y1 - 1,
        q[idx].id = qid, q[idx++].type = 1;
        q[idx].x = x1 - 1, q[idx].y = y2,
        q[idx].id = qid, q[idx++].type = 2;
        q[idx].x = x2, q[idx].y = y1 - 1,
        q[idx].id = qid, q[idx++].type = 3;
        q[idx].x = x2, q[idx].y = y2,
        q[idx].id = qid, q[idx++].type = 4;
        seq[qid] = d[i].qid;
        ans[qid++] = 0;
      } else if (d[i].type == 2) {
        //printf("Add %d %d %d\n", d[i].x, d[i].id, 1);
        q[idx].type = 0;
        q[idx].x = d[i].x;
        q[idx].y = d[i].id;
        q[idx].id = 1;
        who2[d[i].pid] = i;
        idx++;
      } else if (d[i].type == 3) {
        //printf("Del %d %d %d\n", d[i].x, d[i].id, 1);
        q[idx].type = 0;
        q[idx].x = d[who2[d[i].pid]].x;
        q[idx].y = d[who2[d[i].pid]].id;
        q[idx].id = -1;
        idx++;
      }
    }
    M.init(2 * n + 2, 2 * n + 2);
    M.work(q, 0, idx - 1);
    for (int i = 0; i < qid; i++) {
      out[seq[i]] = ans[i];
    }
    for (int i = 0; i < qnum; i++) {
      printf("%d\n", out[i]);
    }
  }
}
// $后缀自动机
namespace SAM {
struct Node {
  Node *ch[11], *fa;
  int val;
  int way;
  int sum;
  Node() : val(0), fa(NULL), way(0), sum(0) {
    memset(ch, 0, sizeof(ch));
  }
} pool[N * 2 + 5], *last, *root;
vector<Node *> vec[N];
namespace SAM {
int cnt;
void init() {
  if (cnt)
    for (int i = 0; i < cnt; i++)
      pool[i] = Node();
  cnt = 1;
  root = &pool[0];
  last = root;
}
void add(int c) {
  Node *p = last, *np = &pool[cnt++];
  last = np;
  np->val = p->val + 1;
  for (; p && !p->ch[c]; p = p->fa)

    p->ch[c] = np;
  if (!p) {
    np->fa = root;
  } else {
    Node *q = p->ch[c];
    if (p->val + 1 == q->val) {
      np->fa = q;
    } else {
      Node *nq = &pool[cnt++];
      *nq = *q;
      nq->val = p->val + 1;
      q->fa = nq;
      np->fa = nq;
      for (; p && p->ch[c] == q; p = p->fa)
        p->ch[c] = nq;
    }
  }
}
}  // namespace SAM
bool cmp(int i, int j) {
  return pool[i].val < pool[j].val;
}
int n, m;
char S[N], buf[N];
void calc() {
  int ans = 0;
  vector<int> vec;
  for (int i = 0; i < SAM::cnt; i++) {
    vec.push_back(i);
    pool[i].way = pool[i].sum = 0;
  }
  sort(vec.begin(), vec.end(), cmp);
  root->way = 1;
  root->sum = 0;
  foreach (it, vec) {
    int i = *it;
    Node *p = &pool[i];
    for (int c = i == 0 ? 1 : 0; c < 10; c++) {
      if (p->ch[c]) {
        p->ch[c]->way += p->way;
        p->ch[c]->way %= Mod;
        p->ch[c]->sum += p->sum * 10 + p->way * c;
        p->ch[c]->sum %= Mod;
      }
    }
    ans += p->sum;
    ans %= Mod;
  }
  printf("%d\n", ans);
}
int main() {
  while (scanf("%d", &n) != EOF) {
    m = 0;
    SAM::init();
    for (int i = 0; i < n; i++) {
      scanf("%s", buf);
      int len = strlen(buf);
      for (int j = 0; j < len; j++) {
        SAM::add(buf[j] - '0');
      }
      SAM::add(10);
    }
    calc();
  }
}
