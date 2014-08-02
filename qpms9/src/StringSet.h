/*
 * StringSet.h
 *
 *  Created on: Jul 4, 2012
 *      Author: marius
 */

#ifndef STRINGSET_H_
#define STRINGSET_H_
#include <string>
#include <vector>
using namespace std;
#include "utils.h"

class StringSet {
public:
	StringSet(int n, int L, int *length, char *memStart) :
			n(n) {
		this->length = new int[n];
		memcpy(this->length, length, n * sizeof(*length));

		this->totalLength = sum(length, n);
		this->memStart = new char[totalLength];
		memcpy(this->memStart, memStart, totalLength * sizeof(*memStart));

		s = new char*[n];
		char *m = this->memStart;
		range = new int[n];
		for (int i = 0; i < n; ++i) {
			s[i] = m;
			m += length[i];
			range[i] = max(0, length[i] - L + 1);
		}
	}

	~StringSet() {
		delete[] length;
		delete[] memStart;
		delete[] range;
		delete[] s;
	}

	int n;
	int totalLength;
	char *memStart;
	int *length;
	char **s;
	int *range;
};

#endif /* STRINGSET_H_ */
