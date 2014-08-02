/*
 * PairDist.h
 *
 *  Created on: Mar 04, 2014
 *      Author: marius
 */

#ifndef PAIRDIST_H_
#define PAIRDIST_H_
#include "CompressedLmers.h"
#include "utils.h"
#include <cstring>

template<class indexType>
class PairDist {
public:
	PairDist(int nRows, int *rSize, indexType **rowItem, int L, int twoD,
			char *minLmerPtr, int nLmers, CompressedLmers *compressedLmers) :
			nLmers(nLmers) {
		if (compressedLmers != NULL)
			prepareCompatiblePairs<true>(nRows, rSize, rowItem, L, twoD,
					minLmerPtr, compressedLmers, nLmers);
		else
			prepareCompatiblePairs<false>(nRows, rSize, rowItem, L, twoD,
					minLmerPtr, compressedLmers, nLmers);
		refDist = NULL;
	}

	virtual ~PairDist() {
		deleteDataStructures();
	}

	int dist(int aInd, int bInd) {
		return pairDist[aInd][bInd];
	}
	/**
	 * Remembers refIndex for future calls to distToReference
	 */
	void setReference(int refInd) {
		refDist = pairDist[refInd];
	}

	/**
	 * Returns the distance between otherIndex and the reference index set with setReference
	 */

	int distToReference(int otherIndex) {
		return refDist[otherIndex];
	}

private:
	char **pairDist;
	char *refDist;

	template<bool areLmersCompressed>
	void precomputeCompatPairs(int *rowSize, indexType **rowItem, int n, int L,
			int twoD, char *minLmerPtr, CompressedLmers *compressedLmers) {
		for (int i = 0; i < n; ++i) {
			int nItemsA = rowSize[i];
			indexType *rowA = rowItem[i];
			for (int j = i + 1; j < n; ++j) {
				int nItemsB = rowSize[j];
				indexType *rowB = rowItem[j];
				for (int k = 0; k < nItemsA; ++k) {
					int aInd = rowA[k];
					for (int l = 0; l < nItemsB; ++l) {
						int bInd = rowB[l];
						int h;
						if (areLmersCompressed) {
							h = compressedLmers->HamDist(aInd, bInd);
						} else {
							h = HammingDist(minLmerPtr + aInd, minLmerPtr + bInd,
									L);
						}
						pairDist[aInd][bInd] = pairDist[bInd][aInd] = h;
					}
				}
			}
		}
	}

	template<bool areLmersCompressed>
	void prepareCompatiblePairs(int n, int *rSize, indexType **rowItem, int L,
			int twoD, char *minLmerPtr, CompressedLmers *compressedLmer,
			int nLmers) {
		clock_t start = clock();
		allocate(pairDist, nLmers, nLmers);
		for (int i = 0; i < nLmers; ++i)
			memset(pairDist[i], 0, nLmers * sizeof(**pairDist));
		precomputeCompatPairs<areLmersCompressed>(rSize, rowItem, n, L, twoD,
				minLmerPtr, compressedLmer);
		clock_t end = clock();
		long time = (end - start) * 1000 / CLOCKS_PER_SEC;
		cerr << "Time to preprocess pairs " << time << "ms" << endl;
	}

	void deleteDataStructures() {
		deAllocate(pairDist, nLmers);
	}

	int nLmers;
};

#endif
