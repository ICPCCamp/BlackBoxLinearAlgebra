// Calculating kth term of linear recurrence sequence
// Complexity: init O(n^2log) query O(n^2logk)
// Requirement: const LOG const MOD
// Input(constructor): vector<int> - first n terms
//                     vector<int> - transition function
// Output(calc(k)): int - the kth term mod MOD
// Example: In: {1, 1} {2, 1} an = 2an-1 + an-2
// 			Out: calc(3) = 3, calc(10007) = 71480733 (MOD 1e9+7)

struct LinearRec {

	int n;
	vector<int> first, trans;
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

	LinearRec(vector<int> &first, vector<int> &trans):first(first), trans(trans) {
		n = first.size();
		vector<int> a(n + 1, 0);
		a[1] = 1;
		bin.push_back(a);
		for (int i = 1; i < LOG; ++i) {
			bin.push_back(add(bin[i - 1], bin[i - 1]));
		}
	}

	int calc(int k) {
		vector<int> a(n + 1, 0);
		a[0] = 1;
		for (int i = 0; i < LOG; ++i) {
			if (k >> i & 1) {
				a = add(a, bin[i]);
			}
		}
		int ret = 0;
		for (int i = 0; i < n; ++i) {
			if ((ret += (long long)a[i + 1] * first[i] % MOD) >= MOD) {
				ret -= MOD;
			}
		}
		return ret;
	}
};
