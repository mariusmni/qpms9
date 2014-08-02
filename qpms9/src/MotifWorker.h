/*
 * MotifWorker.h
 *
 *  Created on: Nov 24, 2012
 *      Author: marius
 */

#ifndef MOTIFWORKER_H_
#define MOTIFWORKER_H_
#include <utility>
#include <algorithm>
#include <vector>
#include <set>
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include "utils.h"
#include "MyString.h"
#include "StringSet.h"
#include "ClosestSub.h"
#include "IndexScheduler.h"
using namespace std;

class MotifWorker {
	typedef uint32 indexType;
public:
	MotifWorker(int myRank, int nProcs, clock_t startTime, int *buffer) :
			myRank(myRank), nProcs(nProcs), startTime(startTime), scheduler(
			NULL) {
		s = NULL;
		comprLmers = NULL;
		compatPairs = NULL;
		loadInputFromBuffer(buffer);
		scheduler = new IndexScheduler(getNTasks() + nProcs);
	}

	~MotifWorker() {
		if (s != NULL)
			delete s;
		if (comprLmers != NULL)
			delete comprLmers;
		if (scheduler != NULL)
			delete scheduler;
		if (compatPairs != NULL)
			delete compatPairs;
	}

	static int *readAndEncodeInput(const char *fileName, int& totalMsgWords,
			MotifConfig& mc) {
		vector<string> strings;
		string alphabet;
		int totalLen = parseInput(fileName, strings, mc, alphabet);
		if (totalLen < 0)
			return NULL;
		else {
			return encodeData(strings, totalLen, mc, alphabet, totalMsgWords);
		}
	}

	void loadInputFromBuffer(int *b) {
		char *bc = (char*) b;
		int n = b[0];

		// read alphabet
		int sigmaLen = b[n + 1];
		int pos = (n + 2) * sizeof(int);
		sigma = string(bc + pos, sigmaLen);
		pos += sigmaLen;

		// read config
		memcpy(&mc, bc + pos, sizeof(MotifConfig));

		// read strings
		pos = (mc.n + 2) * sizeof(int) + mc.sigmaLen + sizeof(MotifConfig);
		s = new StringSet(mc.n, mc.L, b + 1, bc + pos);
	}

	void schedulerLoop() {
		scheduler->loop();
	}

	void doWork() {
		int index;
		long prevPrintTime = 0;
		prepareForWork();
		int nTasks = getNTasks();
		int q = getMinQuorum(mc);

		int *rSize = NULL;
		indexType **rItem = NULL;
		ClosestSub<indexType> *mf = NULL;

		int currentString = -1;
		int totalItemsIncludingCurrentString = 0;
		int allocatedRows = -1;

		while ((index = scheduler->requestWork(myRank)) < nTasks) {
			if (index >= totalItemsIncludingCurrentString) {
				if (rSize != NULL) {
					deAllocate(rSize);
					deAllocate(rItem, allocatedRows);
					set<MyString> m = mf->getMotifs();
					motifs.insert(m.begin(), m.end());
					delete mf;
				}

				do {
					currentString++;
					totalItemsIncludingCurrentString += s->range[currentString];
				} while (index >= totalItemsIncludingCurrentString);

				allocatedRows = mc.n - currentString;
				initLmerMatrices(*s, currentString, rSize, rItem);
				mf = new ClosestSub<indexType>(allocatedRows, q, mc,
						s->memStart, rSize, rItem, comprLmers, compatPairs);
			}

			int indexInCurrentString = index + s->range[currentString]
					- totalItemsIncludingCurrentString;
			mf->processIndex(indexInCurrentString);
			clock_t endInd = clock();
			long globalTime = (endInd - startTime) / CLOCKS_PER_SEC;
			long projected = (long) ((endInd - startTime) * nTasks
					/ ((double) (index + 1)) / CLOCKS_PER_SEC);

			if (globalTime - prevPrintTime > 1) {
				cerr << "\r[P" << myRank << "] Index " << (index + 1)
						<< " out of " << nTasks << "; current time "
						<< (int) globalTime << "s; projected " << projected
						<< "s       ";
				prevPrintTime = globalTime;
			}
		}
		if (rSize != NULL) {
			deAllocate(rSize);
			deAllocate(rItem, allocatedRows);
			set<MyString> m = mf->getMotifs();
			motifs.insert(m.begin(), m.end());
			delete mf;
		}

	}

	set<MyString> getMotifs() {
		return motifs;
	}

	int *allocateMotifBuffer(int nMotifs, int& sz) {
		sz = 1 + ceil(((double) nMotifs * mc.L) / sizeof(int));
		return new int[sz];
	}

	int *encodeMotifs(set<MyString>& motifs, int nMotifs, int& requiredMem) {
		requiredMem = ceil(((double) nMotifs * mc.L) / sizeof(int)) + 1;
		int *buf = new int[requiredMem];
		buf[0] = nMotifs;
		char *bc = ((char *) buf) + sizeof(int);
		set<MyString>::iterator it = motifs.begin();
		for (int i = 0; i < nMotifs; ++i) {
			memcpy(bc, it->s, it->L);
			bc += mc.L;
			++it;
		}
		return buf;
	}

	int decodeMotifs(int *buf, set<MyString>& output) {
		int nm = buf[0];
		char *bc = ((char*) buf) + sizeof(int);
		for (int j = 0; j < nm; ++j, bc += mc.L) {
			output.insert(MyString(bc, mc.L));
		}
		return nm;
	}

	void printMotifs(set<MyString>& motifs) {
		for (set<MyString>::iterator it = motifs.begin(); it != motifs.end();
				++it) {
			std::cout << decodeString(it->s, mc.L, sigma) << std::endl;
		}
	}
private:

	void prepareForWork() {
		int nLmers = s->totalLength - mc.L + 1;
		comprLmers = new CompressedLmers(mc.L, mc.sigmaLen, s->memStart,
				nLmers);
		if (CompatiblePairs::totalBytesRequired(nLmers) <= 0x10000000) { // not too much memory required
			clock_t start = clock();
			compatPairs = new CompatiblePairs(comprLmers, nLmers, 2 * mc.d);
			clock_t end = clock();
			long preprocTime = (end - start) * 1000 / CLOCKS_PER_SEC;
//			cerr << "Time to preprocess pairs " << preprocTime << "ms" << endl;
		}
	}

	void initLmerMatrices(StringSet& s, int startingString, int *&rowSize,
			indexType **&rowItem) {
		int nRows = s.n - startingString;
		allocate(rowSize, nRows);
		allocate(rowItem, nRows);
		for (int i = startingString; i < s.n; ++i) {
			int nItems = s.range[i];
			int j = i - startingString;
			rowSize[j] = nItems;
			indexType *item = allocate(rowItem[j], nItems);
			indexType ind = s.s[i] - s.memStart;
			for (int c = 0; c < nItems; ++c) {
				item[c] = ind + c;
			}
		}
	}

	static int getMinQuorum(MotifConfig& mc) {
		return (int) (mc.q_percent * mc.n / 100);
	}

	int getNTasks() {
		int q = getMinQuorum(mc);
		int maxExcludedStrings = s->n - q;
		int nTasks = 0;
		for (int i = 0; i <= maxExcludedStrings; ++i)
			nTasks += s->range[i];
		return nTasks;
	}

	static bool compare_pairs(pair<int, MyString> a, pair<int, MyString> b) {
		return a.first < b.first;
	}

	void getMyChunk(int originalFirstStringLen, int myRank, int nProcessors,
			int& beginIndex, int& nLmers) {
		int totalLmers = (originalFirstStringLen - mc.L + 1);
		nLmers = totalLmers / nProcessors;
		beginIndex = myRank * nLmers;
		if (myRank < totalLmers % nProcessors) {
			nLmers++;
			beginIndex += myRank;
		} else
			beginIndex += totalLmers % nProcessors;
	}

	MotifConfig mc;

	static void populateThresholds(int totalLen, MotifConfig& mc) {
		double avgStringLen = 1.0 * totalLen / mc.n;
//		mc.sigmaLen = 20;
//		for (int d = 1; d <= 50; ++d) {
//			mc.d = d;
//			mc.L = 2 * mc.d + 5;
//			mc.q_percent=50;

		double tsquare = (2 * mc.d) * log2(mc.sigmaLen) - log2(avgStringLen)
				+ 4;
		double t = round(sqrt(max(9.0, tsquare)));
		if (mc.minStackSize > 0)
			t = mc.minStackSize;
		int quorum = getMinQuorum(mc);
		mc.minStackSize = min(quorum, (int) t);

		double nPrime = t + mc.n / 4 - log2(t);
		if (mc.nPrime > 0)
			nPrime = mc.nPrime;
		int qTolerance = mc.n - quorum;
		mc.nPrime = max(mc.minStackSize + qTolerance, min(mc.n, (int) nPrime));

		cerr << "s=" << mc.minStackSize << " n'=" << mc.nPrime << endl;
//			mc.minStackSize = -1;
//			mc.nPrime = -1;
//		}
//		exit(0);
	}

	static int parseInput(const char *fileName, vector<string>& strings,
			MotifConfig& mc, string& sigma) {
		cerr << "Reading file " << fileName << endl;
		if (readFasta(fileName, strings)) {
			mc.n = strings.size();
			sigma = getAlphabet(strings);
			mc.sigmaLen = sigma.length();
			encodeStrings(strings, sigma);
			int totalLen = sumLengths(strings);
			cerr << "n=" << mc.n << " q=" << mc.q_percent << "% L=" << mc.L
					<< " d=" << mc.d << " totalStringLength=" << totalLen;
			cerr << " alphabet=" << sigma << " of size=" << sigma.length()
					<< endl;
			populateThresholds(totalLen, mc);
			return totalLen;
		} else {
			return -1;
		}
	}

	static int *encodeData(vector<string>& strings, int totalLen,
			MotifConfig& mc, string& sigma, int& totalMsgWords) {
		totalMsgWords = mc.n + 3
				+ (sigma.length() + sizeof(MotifConfig) + totalLen)
						/ sizeof(int);
		int *b = new int[totalMsgWords];

		// write number of strings and string lengths
		b[0] = mc.n;
		for (int i = 0; i < mc.n; ++i)
			b[i + 1] = strings[i].length();

		// write alphabet
		b[mc.n + 1] = sigma.length();
		char *bc = (char*) b;
		int pos = (mc.n + 2) * sizeof(int);
		memcpy(bc + pos, sigma.c_str(), sigma.length());
		pos += sigma.length();

		// write configuration
		memcpy(bc + pos, (char*) &mc, sizeof(MotifConfig));
		pos += sizeof(MotifConfig);

		// write strings
		for (int j = 0; j < mc.n; ++j) {
			int len = strings[j].length();
			memcpy(bc + pos, strings[j].c_str(), len);
			pos += len;
		}
		return b;
	}

	StringSet *s;
	string sigma;
	CompressedLmers *comprLmers;
	CompatiblePairs *compatPairs;
	int myRank;
	int nProcs;
	clock_t startTime;
	IndexScheduler *scheduler;
	set<MyString> motifs;
};

#endif /* MOTIFWORKER_H_ */
