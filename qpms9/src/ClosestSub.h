/*
 * closestString.h
 *
 *  Created on: Jun 21, 2012
 *      Author: marius
 */

#ifndef CLOSESTSUB_H_
#define CLOSESTSUB_H_

#include "utils.h"
#include "MyString.h"
#include "StringSet.h"
#include "CompressedLmers.h"
#include "CompatiblePairs.h"
#include <cstdlib>
#include <cstring>
using namespace std;

//#define DEBUG
#ifdef DEBUG
#define debug(code) code
#else
#define debug(code) do { } while (0)
#endif

#define LARGE_VALUE 0x7fffffff

template<class indexType>
class ClosestSub {
public:

	ClosestSub(int n, int minQuorum, MotifConfig& mc, char *lmerStart,
			int *rSize, indexType **rItem, CompressedLmers *comprLmers,
			CompatiblePairs *compatPairs) :
			n(n), L(mc.L), d(mc.d), quorumTolerance(n - minQuorum), smallN(
					min(n, mc.nPrime)), minStackSize(
					min(n, max(3, mc.minStackSize))), sigmaLen(mc.sigmaLen), twoD(
					2 * d), LminusD(L - d), stackSz(0), compressedLmers(
					comprLmers), compatPairs(compatPairs) {

		initLmerMatrices(rSize, rItem);
		lmers = lmerStart;

		cLmer = new comprWordType[L];

		int cols = smallN * (smallN - 1) / 2;
		allocate(dist, L + 1, cols);
		fill(0, dist[L], cols);

		allocate(lmer, L);

		allocate(r, L + 1, smallN);
		fill(d, r[0], smallN);

		allocate(colMf, smallN + 1, L);
		memset(colMf[0], 0, L * sizeof(*colMf[0]));
		colMaxFreq = colMf[0];

		allocate(colFreq, L, sigmaLen);
		fill(0, colFreq, L, sigmaLen);

		allocate(colS, L, smallN);

		allocate(tempArray, n);

		allocate(maxLB, n);
		for (int i = 0; i < n; ++i) {
			maxLB[i] = i * LminusD;
		}
	}

	~ClosestSub() {
		deleteLmerMatrices();
		deAllocate(cLmer);
		deAllocate(dist, L + 1);
		deAllocate(lmer);
		deAllocate(r, L + 1);
		deAllocate(colMf, smallN + 1);
		deAllocate(colFreq, L);
		deAllocate(colS, L);
		deAllocate(tempArray);
		deAllocate(maxLB);
	}

	char *toNumerical(char *s, char *d, int n) {
		string sigma("ACGT");
		for (int i = 0; i < n; ++i)
			d[i] = sigma.find(s[i]);
		return d;
	}

	bool checkMotifCompr(const comprWordType *m, int nItems, indexType *item,
			CompressedLmers *compLmers) {
		for (int j = 0; j < nItems; ++j) {
			if (compLmers->HamDist(m, item[j]) <= d)
				return true;
		}
		return false;
	}

	bool checkMotifCompressed(comprWordType *cLmer) {
		int incompat = globalNumberOfRows - stackSz;
		for (int i = globalNumberOfRows; i < smallN; ++i)
			if (!checkMotifCompr(cLmer, rowSize[i][stackSz - 1], rowItem[i],
					compressedLmers)) {
				incompat++;
				if (incompat > quorumTolerance)
					return false;
			}
		for (int i = smallN; i < n; ++i)
			if (!checkMotifCompr(cLmer, rowSize[i][2], rowItem[i],
					compressedLmers)) {
				incompat++;
				if (incompat > quorumTolerance)
					return false;
			}
		return true;
	}

	void handle(const char *sol) {
		compressedLmers->compressLmer(sol, cLmer);
		if (checkMotifCompressed(cLmer)) {
			motifs.insert(MyString(sol, L));
		}
	}

	set<MyString> getMotifs() {
		return motifs;
	}

	void processIndex(int i) {
		if (compatPairs == NULL) {
			processInd<false>(i);
		} else {
			processInd<true>(i);
		}
	}

	template<bool hasPreprocessedPairs>
	void processInd(int i) {
		indexType ind = rowItem[0][i];
		addString(lmers + ind);

		if (quorumTolerance > 0) {
			if (getValid<hasPreprocessedPairs, true>(1, n, L - twoD, ind,
					quorumTolerance)) {
				enumerateSubStrings<hasPreprocessedPairs, true>(1,
						quorumTolerance);
			}
		} else {
			if (getValid<hasPreprocessedPairs, false>(1, n, L - twoD, ind, 0)) {
				enumerateSubStrings<hasPreprocessedPairs, false>(1, 0);
			}
		}
		removeLastString();
	}

private:

	template<class K, class V>
	void insert(K *aKey, V *aVal, int l1, int l2) {
		K rk = aKey[l2];
		V rv = aVal[l2];
		int j;
		for (j = l2; j > l1 && aKey[j - 1] > rk; --j) {
			aKey[j] = aKey[j - 1];
			aVal[j] = aVal[j - 1];
		}
		aKey[j] = rk;
		aVal[j] = rv;
	}

	template<bool hasPreprocessedPairs, bool useQ>
	bool getValid(int startRow, int nRows, int extraMaxFreqLB,
			indexType topLmerIndex, int remainQTolerance) {
		comprWordType *topLmer;
		uint32 *pairOk;
		if (hasPreprocessedPairs)
			pairOk = compatPairs->pairOk[topLmerIndex];
		else
			topLmer = compressedLmers->getCompressedLmer(topLmerIndex);
		for (int k = startRow; k < nRows; ++k) {
			indexType *row = rowItem[k];
			int nItems = rowSize[k][stackSz - 1];
			int i = 0;
			int sortStat = 0;
			for (int j = 0; j < nItems; ++j) {
				indexType ind = row[j];
				bool ok;
				if (hasPreprocessedPairs)
					ok = isBitSet(pairOk, ind);
				else
					ok = compressedLmers->HamDist(topLmer, ind) <= twoD;
				if (ok) {
					char *lmerPtr = lmers + ind;
					int m = extraMaxFreq(lmerPtr);
					if (m >= extraMaxFreqLB) {
						swap(row, i++, j);
						if (m > sortStat) {
							/* Take the stack of current kmers and compute for each column, the most frequent characters in that column.
							 For every kmer in a candidate row compute how many columns have the same character as the most frequent
							 character in that column of the stack of current kmers.
							 Let this number be X.
							 For every row compute the largest X. Let this number be Y.
							 Sort rows increasingly by Y.
							 Rows with smallest Y are likely to lead to small neighborhoods.*/
							sortStat = m;
						}
					}
				}
			}

			if (!i) { // empty row
				if (useQ) {
					--remainQTolerance;
					if (remainQTolerance < 0)
						return false;
				} else {
					return false;
				}
			}
			rowSize[k][stackSz] = i;
			sortStat = (sortStat << 20) | i;
			insertRow(tempArray, sortStat, k, startRow);
		}
		return true;
	}

	void insertRow(int *tempArray, int sortStat, int k, int startRow) {
		int *rs = rowSize[k];
		indexType *rv = rowItem[k];
		int j = k;
		for (; j > startRow && tempArray[j - 1] > sortStat; --j) {
			rowSize[j] = rowSize[j - 1];
			rowItem[j] = rowItem[j - 1];
			tempArray[j] = tempArray[j - 1];
		}
		rowSize[j] = rs;
		rowItem[j] = rv;
		tempArray[j] = sortStat;
	}

	bool shouldGenerateNeighborhood(int rowNumber, int itemsInRow) {
		return (stackSz + 1 >= minStackSize);
//				&& (rowNumber + 1 >= smallN || itemsInRow < bruteForceThreshold);
	}

	void bruteForceIt(int rowNumber, int nItems) {
		indexType *row = rowItem[rowNumber];
		globalNumberOfRows = rowNumber + 1;
		for (int j = 0; j < nItems; ++j) {
			int ind = row[j];
			addString(lmers + ind);
			preprocessLowerBounds();
			enumerateStrings(0, maxLB[stackSz] - addMaxFreq());
			removeLastString();
		}
	}

	template<bool hasPreprocessedPairs, bool useQ>
	void enumerateSubStrings(int rowNumber, int remainQTolerance) {
		int nItems = rowSize[rowNumber][stackSz];
		if (shouldGenerateNeighborhood(rowNumber, nItems)) {
			bruteForceIt(rowNumber, nItems);
		} else {
			indexType *row = rowItem[rowNumber];
			for (int j = 0; j < nItems; ++j) {
				indexType ind = row[j];
				addString(lmers + ind);
				preprocessLowerBounds();
				uint threshold = maxLB[stackSz] - addMaxFreq();
				if (hasSolution(0, threshold)) {
					if (getValid<hasPreprocessedPairs, useQ>(rowNumber + 1,
							(stackSz <= 2 ? n : smallN), threshold + LminusD,
							ind, remainQTolerance)) {
						enumerateSubStrings<hasPreprocessedPairs, useQ>(
								rowNumber + 1, remainQTolerance);
					}
				}
				removeLastString();
			}
		}

		if (useQ && remainQTolerance > 0) {
			enumerateSubStrings<hasPreprocessedPairs, useQ>(rowNumber + 1,
					remainQTolerance - 1);
		}
	}

	void initLmerMatrices(int *rSize, indexType **rItem) {
		allocate(rowSize, n, n);
		allocate(rowItem, n);
		for (int i = 0; i < n; ++i) {
			int nItems = rSize[i];
			rowSize[i][0] = nItems;
			allocate(rowItem[i], nItems);
			memcpy(rowItem[i], rItem[i], nItems * sizeof(indexType));
		}
	}

	void deleteLmerMatrices() {
		deAllocate(rowSize, n);
		deAllocate(rowItem, n);
	}

	void addString(const char *t) {
		int *mf = colMf[stackSz + 1];
		for (int j = 0; j < L; ++j) {
			int c = t[j];
			colS[j][stackSz] = c;
			mf[j] = colMaxFreq[j] + (colMaxFreq[j] == colFreq[j][c]++);
		}
		colMaxFreq = mf;
		++stackSz;
	}

	void removeLastString() {
		--stackSz;
		for (int j = 0; j < L; ++j)
			--colFreq[j][colS[j][stackSz]];
		colMaxFreq = colMf[stackSz];
	}

	int addMaxFreq() {
		int total = 0;
		for (int j = 0; j < L; ++j)
			total += colMaxFreq[j];
		return total;
	}

	int extraMaxFreq(const char *t) {
		int total = 0;
		for (int j = 0; j < L; ++j)
			total += (colMaxFreq[j] == colFreq[j][(unsigned char) t[j]]);
		return total;
	}

	void preprocessLowerBounds() {
		int i = stackSz - 1;
		int pairOffset = (i * (i - 1)) >> 1;
		for (int k = L; k; --k) {
			int *dsn = dist[k] + pairOffset;
			int *ds = dist[k - 1] + pairOffset;
			int *s = colS[k - 1];
			char ci = s[i];
			for (int j = 0; j < i; ++j) {
				char cj = s[j];
				*ds++ = (*dsn++) + (ci != cj);
			}
		}
	}

	bool prunePassed(int pos) {
		int *rem = r[pos];
		int *ds = dist[pos];
		for (int i = 1; i < stackSz; ++i) {
			for (int j = 0, ri = rem[i]; j < i; ++j, ++ds)
				if (*ds > ri + rem[j])
					return false;
		}
		return true;
	}

	bool updateRemainingDistances(int *colS, int a, int pos) {
		int *oldRem = r[pos];
		int *newRem = r[pos + 1];
		for (int j = 0; j < stackSz; ++j)
			if ((newRem[j] = oldRem[j] - (a != colS[j])) < 0)
				return false;
		return true;
	}

	bool hasSolution(int pos, int sumFreqLB) {
		if (pos == L)
			return true;
		int *freq = colFreq[pos];
		int *s = colS[pos];
		sumFreqLB += colMaxFreq[pos];
		for (int a = 0; a < sigmaLen; ++a) {
			int f = freq[a];
			if (f && f >= sumFreqLB)
				if (updateRemainingDistances(s, a, pos))
					if (prunePassed(pos + 1))
						if (hasSolution(pos + 1, sumFreqLB - f))
							return true;
		}
		return false;
	}

	void enumerateStrings(int pos, int sumFreqLB) {
		if (pos == L) {
			handle(lmer);
			return;
		}
		int *freq = colFreq[pos];
		int *s = colS[pos];
		sumFreqLB += colMaxFreq[pos];
		for (int a = 0; a < sigmaLen; ++a) {
			int f = freq[a];
			if (f >= sumFreqLB)
				if (updateRemainingDistances(s, a, pos))
					if (prunePassed(pos + 1)) {
						lmer[pos] = a;
						enumerateStrings(pos + 1, sumFreqLB - f);
					}
		}
	}

	const int n, L, d, quorumTolerance, smallN;
	const int minStackSize;
	const int sigmaLen;
	const int twoD, LminusD;
	int stackSz, globalNumberOfRows;
	int **colFreq;
	int *colMaxFreq;
	int **colMf;
	int **colS;

	char *lmer;
	int **r;
	int **dist; // dist[pos][pairNumber] = Hamming distance between suffixes starting at pos of a pair of strings

	int **rowSize;
	indexType **rowItem;
	int *tempArray; // array to store whatever
	int *maxLB;

	set<MyString> motifs;

	char *lmers;

	CompressedLmers * const compressedLmers;
	CompatiblePairs * const compatPairs;

	comprWordType *cLmer;

	void printRows() {
		cerr << "indexType **rItem={" << endl;
		for (int i = 0; i < n; ++i) {
			cerr << "{ ";
			for (int j = 0; j < rowSize[i][0]; ++j) {
				cerr << rowItem[i][j] << ", ";
			}
			cerr << "}," << endl;
		}
		cerr << "};" << endl;
	}

	void loadCustomRows() {
//		indexType rItem[20][580]={{},};
		indexType rItem[1][0] = { { }, };
		for (int i = 0; i < n; ++i) {
			int nItems = rowSize[i][0];
			memcpy(rowItem[i], rItem[i], nItems * sizeof(indexType));
		}
	}

}
;

#endif
