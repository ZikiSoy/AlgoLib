int idx(char c) { 
    return c - 'a';
}

struct AhoCorasick {
  static const int maxnode = 26 * 100005;
  static const int sigma_size = 26;

  int ch[maxnode][sigma_size]; //child
  int f[maxnode]; //fail 
  int64 val[maxnode]; //value
  int match[maxnode]; //max match
  int len[maxnode]; //length of a node
  int sz; //size
  void reset() { sz = 1; memset(ch[0], 0, sizeof(ch[0])); }

  void insert(const char *s, int v) {
      int u = 0, n = strlen(s);
      for (int i = 0; i < n; i++) {
          int c = idx(s[i]);
          if (!ch[u][c]) {
              memset(ch[sz], 0, sizeof(ch[sz]));
              match[sz] = val[sz] = 0;
              len[sz] = len[u] + 1;
              dout << "insert " << s[i] << ' ' <<sz << endl;
              ch[u][c] = sz++;
          }
          u = ch[u][c];
      }

      match[u] = n;
      val[u] = v;
      dout << "value " << u << ' ' << v << endl;
  }

  void getFail() {
      queue<int> q;
      f[0] = 0;

      for (int c = 0; c < sigma_size; c++) {
          int u = ch[0][c];
          if (u) { f[u] = 0; q.push(u); }
      }

      while (!q.empty()) {
          int r = q.front(); q.pop();
          for (int c = 0; c < sigma_size; c++) {
              int &v = ch[r][c];
              if (!v) { 
                  v = ch[f[r]][c]; 
              } else {
                  q.push(v);
                  f[v] = ch[f[r]][c];
                  // How to combine value with same suffix.
                  // val[v] = val[v]|val[f[v]];
                  val[v] = val[v] + val[f[v]];
                  dout << v << ' ' << val[v] << endl;
                  match[v] = max(match[v], match[f[v]]);
              }
          }
      }
  }
  // Find string starting from node j.
  int64 find(const char* T, int j = 0) {
      int n = strlen(T);
      int64 ret = 0;
      for (int i = 0; i < n; i++) {
          int c = idx(T[i]);
          while (j && !ch[j][c]) j = f[j];
          j = ch[j][c];
          ret += val[j];
      }
      return ret;
  }
} ac;
