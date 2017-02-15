#include<vector>
#include<cstdio>
#include<cstring>
#include<iostream>
#include<algorithm>
using namespace std;

const unsigned LOG = 31;

// Calculating kth term of linear recurrence sequence
// Modified for testing

struct LinearRec {

	int n;
	vector<unsigned> first, trans;
	vector<vector<unsigned> > bin;

	vector<unsigned> add(vector<unsigned> &a, vector<unsigned> &b) {
		vector<unsigned> result(n * 2 + 1, 0);
		for (int i = 0; i <= n; ++i) {
			for (int j = 0; j <= n; ++j) {
				result[i + j] ^= a[i] & b[j];
			}
		}
		for (int i = 2 * n; i > n; --i) {
			for (int j = 0; j < n; ++j) {
				result[i - 1 - j] ^= result[i] & trans[j];
			}
			result[i] = 0;
		}
		result.erase(result.begin() + n + 1, result.end());
		return result;
	}

	LinearRec(vector<unsigned> &first, vector<unsigned> &trans):first(first), trans(trans) {
		n = first.size();
		vector<unsigned> a(n + 1, 0);
		a[1] = ~0u;
		bin.push_back(a);
		for (int i = 1; i < LOG; ++i) {
			bin.push_back(add(bin[i - 1], bin[i - 1]));
		}
	}

	unsigned calc(int k) {
		vector<unsigned> a(n + 1, 0);
		a[0] = ~0u;
		for (int i = 0; i < LOG; ++i) {
			if (k >> i & 1) {
				a = add(a, bin[i]);
			}
		}
		unsigned ret = 0;
		for (int i = 0; i < n; ++i) {
			ret = ret ^ (a[i + 1] & first[i]);
		}
		return ret;
	}
};

//end of template

//test on http://abc009.contest.atcoder.jp/tasks/abc009_4

int n, k;

int main() {
	scanf("%d%d", &n, &k);
	vector<unsigned> a(n), b(n);
	for (int i = 0; i < n; ++i) {
		scanf("%u", &a[i]);
	}
	for (int i = 0; i < n; ++i) {
		scanf("%u", &b[i]);
	}
	LinearRec f(a, b);
	printf("%u\n", f.calc(k));
	return 0;
}
