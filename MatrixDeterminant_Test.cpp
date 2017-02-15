#include<vector>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<algorithm>
using namespace std;

const int MOD = 1000000007, N = 10005, M = N * 11;

const int BAR = 10;

struct Vector {
	int n, a[N];

	Vector(int n):n(n) {
		memset(a, 0, sizeof(a));
	}
	
	int& operator[] (const int &i) { return a[i]; }
	const int operator[] (const int &i) const { return a[i]; }

	int operator * (const Vector &b) {
		unsigned long long ret = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < BAR && i < n; ++j, ++i) {
				ret = ret + (unsigned long long)a[i] * b[i];
			}
			--i;
			ret %= MOD;
		}
		return ret;
	}

};

struct Matrix {
	int n, m;
	int x[M], y[M], a[M];

	Matrix(int n):n(n) {
		m = 0;
		memset(a, 0, sizeof(a));
	}

	void reshuffle() {
		vector<pair<int, pair<int, int> > > v(m);
		for (int i = 0; i < m; ++i) {
			v[i].first = x[i], v[i].second.first = y[i], v[i].second.second = a[i];
		}
		sort(v.begin(), v.end());
		for (int i = 0; i < m; ++i) {
			x[i] = v[i].first, y[i] = v[i].second.first, a[i] = v[i].second.second;
		}
	}

	Vector operator * (const Vector &b) const {
		Vector ret(n);
		for (int i = 0; i < m; ++i) {
			if ((ret[x[i]] += (unsigned long long)a[i] * b[y[i]] % MOD) >= MOD) {
				ret[x[i]] -= MOD;
			}
		}
		return ret;
	}
};

unsigned long long buf[N];

void mul(const Matrix &A, Vector &b) { //to save memory
	int n = A.n;
	memset(buf, 0, sizeof(unsigned long long) * n);

	for (int i = 0; i < A.m; ++i) {
		buf[A.x[i]] += (unsigned long long)A.a[i] * b[A.y[i]];
		if (i % BAR == 0) {
			buf[A.x[i]] %= MOD;
		}
	}

	for (int i = 0; i < A.n; ++i) {
		b[i] = buf[i] % MOD;
	}
}

// Berlekamp-Massey Algorithm

int inverse(int a) {
	return a == 1 ? 1 : (long long)(MOD - MOD / a) * inverse(MOD % a) % MOD;
}

vector<int> na;

struct Poly {
	vector<int> a;

	Poly() { a.clear(); }

	Poly(vector<int> &a): a(a) {}

	int length() const { return a.size(); }

	Poly move(int d) {
		na.resize(d + a.size());
		for (int i = 0; i < d + a.size(); ++i) {
			na[i] = i < d ? 0 : a[i - d];
		}
		return na;
	}

	int calc(vector<int> &d, int pos) {
		unsigned long long ret = 0;
		for (int i = 0; i < (int)a.size(); ++i) {
			for (int j = 0; j < BAR && i < (int)a.size(); ++j, ++i) {
				ret = ret + (unsigned long long)d[pos - i] * a[i];
			}
			--i;
			ret %= MOD;
		}
		return ret;
	}

	Poly operator - (const Poly &b) {
		na.resize(max(this->length(), b.length()));
		for (int i = 0; i < (int)na.size(); ++i) {	
			int aa = i < this->length() ? this->a[i] : 0,
				bb = i < b.length() ? b.a[i] : 0;
			na[i] = aa >= bb ? aa - bb : aa + MOD - bb;
		}
		return Poly(na);
	}
};

Poly operator * (const int &c, const Poly &p) {
	na.resize(p.length());
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

Vector getRandomVector(int n) {
	Vector ret(n);
	for (int i = 0; i < n; ++i) {
		ret[i] = rand() % MOD;
	}
	return ret;
}

int solve(Matrix &A) {
	Vector d = getRandomVector(A.n), x = getRandomVector(A.n), y = getRandomVector(A.n);
	for (int i = 0; i < A.m; ++i) {
		A.a[i] = (long long)A.a[i] * d[A.x[i]] % MOD;
	}
	vector<int> a;
	for (int i = 0; i < A.n * 2 + 1; ++i) {
		mul(A, x); //x = A * x;
		a.push_back(x * y);
	}
	vector<int> s = Berlekamp(a);
	int ret = s.back();
	if (A.n & 1) {
		ret = (MOD - ret) % MOD;
	}
	for (int i = 0; i < A.n; ++i) {
		ret = (long long)ret * inverse(d[i]) % MOD;
	}
	return ret;
}

//tested on CF 348F - Little Artem and Graph

int n, k;

void initMatrix(Matrix &A) {
	A.m = n - 1;
	for (int i = 0; i < n - 1; ++i) {
		A.x[i] = A.y[i] = i;
		A.a[i] = 0;
	}
}

void addEdge(Matrix &A, int u, int v) {
	if (u < A.n && v < A.n) {
		A.x[A.m] = u, A.y[A.m] = v, A.a[A.m] = MOD - 1, ++A.m;
		A.x[A.m] = v, A.y[A.m] = u, A.a[A.m] = MOD - 1, ++A.m;
	}
	if (u < A.n) {
		++A.a[u];
	}
	if (v < A.n) {
		++A.a[v];
	}
}

int main() {
	scanf("%d%d", &n, &k);
	Matrix A(n - 1);
	initMatrix(A);
	for (int i = 0; i < k; ++i) {
		for (int j = i + 1; j < k; ++j) {
			addEdge(A, i, j);
		}
	}
	for (int i = k; i < n; ++i) {
		int u = i, v;
		for (int j = 0; j < k; ++j) {
			scanf("%d", &v);
			--v;
			addEdge(A, u, v);
		}
	}
	A.reshuffle();
	int ans = solve(A);
	printf("%d\n", ans);
	return 0;
}
