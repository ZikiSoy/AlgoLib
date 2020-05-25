//  https://www.codechef.com/submit/CHEFSHIP

namespace SA {
// $后缀数组+RMQ
// sa[] 将后缀 0~n-1 放入 sa 数组按后缀字典序排序后的数组
// rk[] 后缀 i 在 sa 中的下标

// height[i] sa[i]和sa[i-1]的最长公共前缀
// lcp(x, y) s[x]和s[y]的最长公共前缀
const int N = 3e5 + 10;
int s[N];
int sa[N], t[N], t2[N], c[N], n;
int rk[N], height[N];
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
  for (i = 0; i < n; i++) height[i] = -1;
  for (i = 0; i < n; i++) rk[sa[i]] = i;
  for (i = 0; i < n; i++) {
    if (k) k--;
    if (rk[i] == 0) continue;
    int j = sa[rk[i] - 1];
    while (s[i + k] == s[j + k] && i+k<n && j+k<n) k++;
    height[rk[i]] = k;
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
int rmq(int x, int y) {  //choose the best in { val[x] ... val[y] }
  int bl = ln[y - x + 1];
  return choose(st[bl][x], st[bl][y - (1 << bl) + 1]);
}
template<typename T> 
void init(T str, int len) {
    int min_char = 1e9, max_char = -1;
    for (int i = 0; i < len; i++) {
        s[i] = str[i];
        if (min_char > s[i]) min_char = s[i];
        if (max_char < s[i]) max_char = s[i];
    }
    for (int i = 0; i < len; i++) s[i] -= min_char;
    n = len;
    build_sa(max_char-min_char+1);
    getHeight();
    initrmq(len, height);
}
// longest common prefix between suffix x and suffix y
// 0 based. x and y should be different.
int lcp(int x, int y) {
    // if (x == y) { dout << "XY is equal " << x << ' ' << y << endl; }
    return rmq(min(rk[x], rk[y])+1, max(rk[x], rk[y]));
}
}  //namespace SA
