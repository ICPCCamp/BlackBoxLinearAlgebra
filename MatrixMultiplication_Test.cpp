#include<vector>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<algorithm>
using namespace std;

const int MOD = 1000000007, M = 52, LOG = 63;

int m; // dimension

struct Vector {
	int a[M];

	Vector() {
		memset(a, 0, sizeof(a));
	}
	
	int& operator[] (const int &i) { return a[i]; }
	const int operator[] (const int &i) const { return a[i]; }

	int operator * (const Vector &b) {
		int ret = 0;
		for (int i = 0; i < m; ++i) {
			if ((ret += (long long)a[i] * b[i] % MOD) >= MOD) {
				ret -= MOD;
			}
		}
		return ret;
	}

	Vector operator + (const Vector &b) {
		Vector ret;
		for (int i = 0; i < m; ++i) {
			if ((ret[i] = a[i] + b[i]) >= MOD) {
				ret[i] -= MOD;
			}
		}
		return ret;
	}

};

Vector operator * (int k, const Vector &b) {
	Vector ret;
	for (int i = 0; i < m; ++i) {
		ret[i] = (long long)k * b[i] % MOD;
	}
	return ret;
}

// Reimplement this structure to support sparse matrix
struct Matrix {
	int a[M][M];

	int* operator[]	(const int &i) { return a[i]; }
	const int* operator[] (const int &i) const { return a[i]; }

	Vector operator * (const Vector &b) {
		Vector ret;
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < m; ++j) {
				if ((ret[i] += (long long)a[i][j] * b[j] % MOD) >= MOD) {
					ret[i] -= MOD;
				}
			}
		}
		return ret;
	}
};

// Berlekamp-Massey Algorithm

int inverse(int a) {
	return a == 1 ? 1 : (long long)(MOD - MOD / a) * inverse(MOD % a) % MOD;
}

struct Poly {
	vector<int> a;

	Poly() { a.clear(); }

	Poly(vector<int> &a): a(a) {}

	int length() const { return a.size(); }

	Poly move(int d) {
		vector<int> na(d, 0);
		na.insert(na.end(), a.begin(), a.end());
		return Poly(na);
	}

	int calc(vector<int> &d, int pos) {
		int ret = 0;
		for (int i = 0; i < (int)a.size(); ++i) {
			if ((ret += (long long)d[pos - i] * a[i] % MOD) >= MOD) {
				ret -= MOD;
			}	
		}
		return ret;
	}

	Poly operator - (const Poly &b) {
		vector<int> na(max(this->length(), b.length()));
		for (int i = 0; i < (int)na.size(); ++i) {	
			int aa = i < this->length() ? this->a[i] : 0,
				bb = i < b.length() ? b.a[i] : 0;
			na[i] = (aa + MOD - bb) % MOD;
		}
		return Poly(na);
	}
};

Poly operator * (const int &c, const Poly &p) {
	vector<int> na(p.length());
	for (int i = 0; i < (int)na.size(); ++i) {
		na[i] = (long long)c * p.a[i] % MOD;
	}
	return na;
}

vector<int> Berlekamp(vector<int> a) {
	int n = a.size();
	Poly s, b;
	s.a.push_back(1), b.a.push_back(1);
	for (int i = 1, j = 0, ld = a[0]; i < n; ++i) {
		int d = s.calc(a, i);
		if (d) {
			if ((s.length() - 1) * 2 <= i) {
				Poly ob = b;
				b = s;
				s = s - (long long)d * inverse(ld) % MOD * ob.move(i - j);
				j = i;
				ld = d;
			} else {
				s = s - (long long)d * inverse(ld) % MOD * b.move(i - j);
			}
		}
	}
	return s.a;
}

// Calculating kth term of linear recurrence sequence
// Modified for application

struct LinearRec {

	int n;
	vector<Vector> first; 
	vector<int> trans;
	vector<vector<int> > bin;

	vector<int> add(vector<int> &a, vector<int> &b) {
		vector<int> result(n * 2 + 1, 0);
		//You can apply constant optimization here to get a ~10x speedup
		for (int i = 0; i <= n; ++i) {
			for (int j = 0; j <= n; ++j) {
				if ((result[i + j] += (long long)a[i] * b[j] % MOD) >= MOD) {
					result[i + j] -= MOD;
				}
			}
		}
		for (int i = 2 * n; i > n; --i) {
			for (int j = 0; j < n; ++j) {
				if ((result[i - 1 - j] += (long long)result[i] * trans[j] % MOD) >= MOD) {
					result[i - 1 - j] -= MOD;
				}
			}
			result[i] = 0;
		}
		result.erase(result.begin() + n + 1, result.end());
		return result;
	}

	LinearRec(vector<Vector> &first, vector<int> &trans):first(first), trans(trans) {
		n = first.size();
		vector<int> a(n + 1, 0);
		a[1] = 1;
		bin.push_back(a);
		for (int i = 1; i < LOG; ++i) {
			bin.push_back(add(bin[i - 1], bin[i - 1]));
		}
	}

	Vector calc(long long k) {
		vector<int> a(n + 1, 0);
		a[0] = 1;
		for (int i = 0; i < LOG; ++i) {
			if (k >> i & 1) {
				a = add(a, bin[i]);
			}
		}
		Vector ret;
		for (int i = 0; i < n; ++i) {
			ret = ret + a[i + 1] * first[i];
		}
		return ret;
	}
};

Vector solve(Matrix &A, long long k, Vector &b) {
	vector<Vector> vs;
	vs.push_back(A * b);
	for (int i = 1; i < m * 2; ++i) {
		vs.push_back(A * vs.back());
	}
	if (k == 0) {
		return b;
	} else if (k <= m * 2) {
		return vs[k - 1];
	}
	Vector x;
	for (int i = 0; i < m; ++i) {
		x[i] = rand() % MOD;
	}
	vector<int> a(m * 2);
	for (int i = 0; i < m * 2; ++i) {
		a[i] = 	vs[i] * x;
	}
	vector<int> s = Berlekamp(a);
	s.erase(s.begin());
	for (int i = 0; i < s.size(); ++i) {
		s[i] = (MOD - s[i]) % MOD;
	}
	vector<Vector> vf(vs.begin(), vs.begin() + s.size());
	LinearRec f(vf, s);
	return f.calc(k);
}

//tested on CF 222E - Decoding Genome
int n;

long long k;

int main() {
	scanf("%lld %d %d", &k, &m, &n);
	Matrix A;
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			A[i][j] = 1;
		}
	}
	for (int i = 0; i < n; ++i) {
		char s[3];
		scanf("%s", s);
		int t1 = 'a' <= s[0] && s[0] <= 'z' ? t1 = s[0] - 'a' : t1 = s[0] - 'A' + 26;
		int t2 = 'a' <= s[1] && s[1] <= 'z' ? t2 = s[1] - 'a' : t2 = s[1] - 'A' + 26;
		A[t1][t2] = 0;
	}
	Vector b;
	for (int i = 0; i < m; ++i) {
		b[i] = 1;
	}
	int ans = solve(A, k - 1, b) * b;
	printf("%d\n", ans);
	return 0;
}
