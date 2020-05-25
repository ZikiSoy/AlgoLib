namespace NumberTheory {

// Generate primes[], phi[Z]
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

// Greatest Common Divisor
int64 gcd(int64 a, int64 b) {
  return b ? gcd(b, a % b) : a;
}
int64 f(int64 x, int64 n) {
  return (mul(x, x, n) + 1) % n;
}

// Test whether n is a prime number.
// 判定n是否是质数
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

// a * x + b * y = gcd(a, b)
int64 extGcd(int64 a, int64 b, int64& x, int64& y) {
  if (b == 0) {
    x = 1;
    y = 0;
    return a;
  } else {
    int g = extGcd(b, a % b, y, x);
    y -= a / b * x;
    return g;
  }
}

// ASSUME: gcd(a, m) == 1
int64 modInv(int64 a, int64 m) {
  int64 x, y;
  extGcd(a, m, x, y);
  return (x % m + m) % m;
}

// Get a dividor.
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

// Get a prime dividor for a.
int64 PrimeDividor(int64 a) {
  if (Miller_Rabin(a)) return 0;
  int64 t = Pollard_Rho(a);
  int64 p = Prime(t);
  if (p) return p;
  return t;
}

// 质因数分解
// 所有素约数保存于 factor[]中, m 为其个数, 返回分解后的素因子个数
int factorize(int64 x, int64 factor[]) {
  int m = 0;
  while (x > 1) {
    if (Miller_Rabin(x)) break;
    int64 t = PrimeDividor(x);
    factor[m++] = t;
    x /= t;
  }
  if (x > 0) factor[m++] = x;
  //所有素约数保存于 factor[]中,m 为其个数
  return m;
}

// Chinese Remainder Theory

// Lucas Theory
// Combination % mod  [mod is a prime number]
// ** O(mod) memory + O(mod) preprocessing
int64 pow(int64 n, int64 p, int64 mod) {

}
// n < mod && m < mod
int64 comb(int64 n, int64 m, int64 mod) {

}
int64 lucas(int64 n, int64 m, int64 mod) {
  return m ? lucas(n/p, m/p, p) * comb(n%p, m%p, p) % p : 1;
}

};