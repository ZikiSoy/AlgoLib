#include <queue>
#include <algorithm>

using namespace std;

int idx(char c) { 
	if (c == 'A') return 0;
	if (c == 'T') return 1;
	if (c == 'C') return 2;
	if (c == 'G') return 3;
}

struct AhoCorasick {
  static const int maxnode = 200005;
  static const int sigma_size = 26;
  
	int ch[maxnode][sigma_size]; //child
	int f[maxnode]; //fail 
	int val[maxnode]; //value
	int match[maxnode]; //max match
	int len[maxnode]; //length of a node
	int sz; //size
	void reset() { sz = 1; memset(ch[0], 0, sizeof(ch[0])); }

	void insert(char *s, int v) {
		int u = 0, n = strlen(s);
		for (int i = 0; i < n; i++) {
			int c = idx(s[i]);
			if (!ch[u][c]) {
				memset(ch[sz], 0, sizeof(ch[sz]));
				match[sz] = val[sz] = 0;
				len[sz] = len[u] + 1;
				ch[u][c] = sz++;
			}
			u = ch[u][c];
		}

		match[u] = n;
		val[u] = v;
	}

	int getFail() {
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
					val[v] = val[v]|val[f[v]];
					match[v] = max(match[v], match[f[v]]);
				}
			}
		}
	}
} ac;

