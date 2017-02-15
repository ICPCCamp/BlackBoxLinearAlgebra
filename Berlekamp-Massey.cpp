// Berlekamp-Massey Algorithm
// Complexity: O(n^2)
// Requirement: const MOD, inverse(int)
// Input: vector<int> - the first elements of the sequence
// Output: vector<int> - the recursive equation of the given sequence
// Example:  In: {1, 1, 2, 3}  Out: {1, 1000000006, 1000000006} (MOD = 1e9+7)

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

vector<int> solve(vector<int> a) {
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
	//Caution: s.a might be shorter than expected
	return s.a;
}
