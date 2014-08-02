/*
 * CompressedLmers.h
 *
 *  Created on: Nov 19, 2012
 *      Author: marius
 */

#ifndef COMPRESSEDLMERS_H_
#define COMPRESSEDLMERS_H_
#include "utils.h"

typedef uint16 comprWordType;

class CompressedLmers {
	typedef uint16 preprocHDType;
	static const int PREPROC_HD_SIZE = 1 << (8 * sizeof(preprocHDType));

public:
	static const unsigned bitsPerComprWord = 8 * sizeof(comprWordType);

	CompressedLmers(int L, int sigmaLen, char *s = NULL, int sLength = 0) :
			bitsPerSymbol(bitsFor(sigmaLen - 1)), symbolsPerWord(
					bitsPerComprWord / bitsPerSymbol), symbolsInLastWord(
					L % symbolsPerWord == 0 ?
							symbolsPerWord : L % symbolsPerWord),
					wordsPerLmer(
					(L / symbolsPerWord) + (L % symbolsPerWord ? 1 : 0)) {

		preprocessDistances();
		reservedSize = 0;
		compressedLmer = NULL;
		if (s != NULL)
			compressAllLmers(s, sLength);
	}

	virtual ~CompressedLmers() {
		deleteCompressedLmers();
	}

	int HamDist(int indA, int indB) const {
		return HamDist(getCompressedLmer(indA), getCompressedLmer(indB));
	}

	int HamDist(const comprWordType *ca, int ind) const {
		return HamDist(ca, getCompressedLmer(ind));
	}

	int HamDist(const comprWordType *ca, const comprWordType *cb) const {
		int d = 0;
		for (unsigned i = 0; i < wordsPerLmer; ++i) {
			comprWordType a = ca[i] ^ cb[i];
			d += compressedMismatches[a];
		}
		return d;
	}

	template<typename T>
	void compressLmer(const T *l, comprWordType *cl) {
		comprWordType w;
		for (unsigned i = 1; i < wordsPerLmer; ++i) {
			w = 0;
			for (unsigned j = 0; j < symbolsPerWord; ++j)
				w = (w << bitsPerSymbol) | *l++;
			*cl++ = w;
		}
		w = 0;
		for (unsigned j = 0; j < symbolsInLastWord; ++j) {
			w = (w << bitsPerSymbol) | *l++;
		}
		*cl = w;
	}

	comprWordType *getCompressedLmer(int index) const {
		return compressedLmer + index * wordsPerLmer;
	}

	void reserveAndClear(int nLmers) {
		if (nLmers > reservedSize) {
			if (compressedLmer != NULL)
				delete[] compressedLmer;
			int reservedLen = nextPow2(nLmers * wordsPerLmer);
			reservedSize = reservedLen / wordsPerLmer;
			compressedLmer = new comprWordType[reservedLen];
		}
	}

	void copy(int dst, comprWordType *src) {
		memcpy(getCompressedLmer(dst), src,
				wordsPerLmer * sizeof(comprWordType));
	}

private:

	void compressAllLmers(char *lmers, int nLmers) {
		compressedLmer = new comprWordType[nLmers * wordsPerLmer];
		for (int i = 0, j = 0; i < nLmers; ++i, j += wordsPerLmer)
			compressLmer(lmers + i, compressedLmer + j);
		reservedSize = nLmers;
	}

	int countNonZeroGroups(preprocHDType w, preprocHDType lastSymbolMask,
			unsigned bitsPerSymbol) {
		int d = 0;
		for (; w;) {
			d += 0 != (w & lastSymbolMask);
			w >>= bitsPerSymbol;
		}
		return d;
	}

	void preprocessDistances() {
		preprocHDType rightMask = getLSOnesMask(bitsPerSymbol);
		for (int i = 0; i < PREPROC_HD_SIZE; ++i) {
			compressedMismatches[i] = countNonZeroGroups(i, rightMask,
					bitsPerSymbol);
		}
	}

	void deleteCompressedLmers() {
		if (compressedLmer != NULL)
			delete[] compressedLmer;
	}

	comprWordType getLSOnesMask(int nBits) {
		return ~((~0) << nBits);
	}

	comprWordType *compressedLmer;
	char compressedMismatches[PREPROC_HD_SIZE];
	const unsigned bitsPerSymbol;
	const unsigned symbolsPerWord;
	const unsigned symbolsInLastWord;
	const unsigned wordsPerLmer;
	int reservedSize;
}
;

#endif /* COMPRESSEDLMERS_H_ */
