/*
 * RandomGen.h
 *
 *  Created on: Mar 11, 2014
 *      Author: marius
 */

#ifndef RANDOMGEN_H_
#define RANDOMGEN_H_
#include <string>
#include <cstdlib>

using namespace std;

class RandomGen {
public:
	RandomGen(int seed) {
		srand(seed);
	}

	char randomChar(string& alphabet) {
		return alphabet[rand() % alphabet.size()];
	}

	int randomPos(int length) {
		return rand() % length;
	}

	void generateDistinctPositions(int q, int n, bool *posFlag) {
		for (int i = 0; i < n; ++i)
			posFlag[i] = false;
		for (int i = 0; i < q;) {
			int r = randomPos(n);
			if (!posFlag[r]) {
				posFlag[r] = true;
				++i;
			}
		}
	}

	string generateString(int length, string& alphabet) {
		char s[length];
		for (int i = 0; i < length; ++i)
			s[i] = randomChar(alphabet);
		return string(s, length);
	}

	char mutate(char c, string& alphabet) {
		char d;
		do {
			d = randomChar(alphabet);
		} while (d == c);
		return d;
	}

	string mutate(string s, int d, string& alphabet) {
		bool modify[s.length()];
		generateDistinctPositions(d, s.length(), modify);
		for (unsigned i = 0; i < s.length(); ++i)
			if (modify[i]) {
				s[i] = mutate(s[i], alphabet);
			}
		return s;
	}

private:
	template<class T>
	bool find(T key, T *array, int length) {
		for (int i = 0; i < length; ++i)
			if (key == array[i])
				return true;
		return false;
	}
};

#endif /* RANDOMGEN_H_ */
