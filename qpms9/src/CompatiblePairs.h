/*
 * CompatiblePairs.h
 *
 *  Created on: Nov 19, 2012
 *      Author: marius
 */

#ifndef COMPATIBLEPAIRS_H_
#define COMPATIBLEPAIRS_H_
#include "CompressedLmers.h"
#include "utils.h"
#include <cstring>
#include <vector>
#include <map>
using namespace std;

class CompatiblePairs {
private:
	const CompressedLmers *compressedLmers;
	const int nLmers, hdLimit, nBytesInRow;
	int reservedSize;

public:
	CompatiblePairs(CompressedLmers *compressedLmers, int nLmers, int hdLimit) :
			compressedLmers(compressedLmers), nLmers(nLmers), hdLimit(hdLimit), nBytesInRow(
					(nLmers / (8 * sizeof(**pairOk)) + 1) * sizeof(**pairOk)) {
		allocateDataStructures(nLmers);
		init();
	}

	~CompatiblePairs() {
		deleteDataStructures();
	}

	static long totalBytesRequired(int nLmers) {
		long nRows = nextPow2(nLmers);
		long rowSize = nRows / (8 * sizeof(uint32)) + 1;
		return nRows * rowSize * sizeof(uint32);
	}

	uint32 **pairOk;

private:

	void init() {
		for (int i = 0; i < nLmers; ++i) {
			memset(pairOk[i], 0, nBytesInRow);
			comprWordType *a = compressedLmers->getCompressedLmer(i);
			uint32 *p = pairOk[i];
			for (int j = 0; j < i; ++j)
				if (compressedLmers->HamDist(a, j) <= hdLimit) {
					setBit(p, j);
					setBit(pairOk[j], i);
				}
		}
	}

	void allocateDataStructures(int nLmers) {
		reservedSize = nextPow2(nLmers);
		int rowSize = reservedSize / (8 * sizeof(**pairOk)) + 1;
		allocate(pairOk, reservedSize, rowSize);
	}

	void deleteDataStructures() {
		deAllocate(pairOk, reservedSize);
	}

};

#endif /* COMPATIBLEPAIRS_H_ */
