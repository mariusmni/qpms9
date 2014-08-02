/*
 * utils.h
 *
 *  Created on: Jun 9, 2012
 *      Author: marius
 */

#ifndef UTILS_H_
#define UTILS_H_

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long uint64;

typedef char int8;
typedef short int16;
typedef int int32;
typedef long int64;

struct MotifConfig {
	int n, L, d, sigmaLen;
	float q_percent;
	int nPrime, minStackSize;
};

#include <string>
#include <vector>
#include <iostream>
using namespace std;

string getAlphabet(vector<string>& strings);

string getAlphabetC(char** s, int n, int m);

char** readInputC(char *fileName, int& n, int& L, int& d);

template<class T>
T* allocate(T*&dest, int n) {
	dest = new T[n];
	return dest;
}

template<class T>
T** allocate(T**&dest, int n, int m) {
	T **a = new T*[n];
	for (int i = 0; i < n; ++i)
		allocate(a[i], m);
	dest = a;
	return dest;
}

template<class T>
T*** allocate(T***&dest, int n, int m, int k) {
	T ***a = new T**[n];
	for (int i = 0; i < n; ++i)
		allocate(a[i], m, k);
	dest = a;
	return dest;
}

char **allocate2DBlockChar(int n, int m); // n x m matrix stored in single contiguous memory

char **allocate2DBlockCharB(int n, int m, int bufSize);// n x m matrix stored in single contiguous memory of specified size (has to fit)

int HammingDist(const char *a, const char *b, int n);

bool isHamDistWithin(char *a, char *b, int n, int d);

int charCmp(char *a, char *b, int L);

char *getCopy(char *a, int n);

/** Returns the minimum number of bits needed to represent number n */
int bitsFor(int n);

/** Returns the smallest power of two greater or equal to n */
int nextPow2(int n);

template<class T>
void deAllocate(T **a, int n) {
	if (a == NULL)
		return;
	for (int i = 0; i < n; ++i)
		delete[] a[i];
	delete[] a;
}

template<class T>
void deAllocate2DBlock(T **a) {
	if (a == NULL)
		return;
	delete[] a[0]; // delete memory
	delete[] a; // delete pointers
}

template<class T>
void deAllocate(T *a) {
	delete[] a;
}

template<class T>
void deAllocate(T ***a, int n, int m) {
	for (int i = 0; i < n; ++i)
		deAllocate(a[i], m);
	delete[] a;
}

template<class T>
void deAllocate3DWithBlock23(T ***a, int n) {
	for (int i = 0; i < n; ++i)
		deAllocate2DBlock(a[i]);
	delete[] a;
}

template<class T> void swap(T *s, int i, int j) {
	T t = s[i];
	s[i] = s[j];
	s[j] = t;
}

template<class T>
void fill(T val, T *a, int n) {
	for (int i = 0; i < n; ++i)
		a[i] = val;
}

template<class T>
void fill(T val, T **a, int n, int m) {
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			a[i][j] = val;
}

template<class T>
void fill(T val, T ***a, int n, int m, int k) {
	for (int i = 0; i < n; ++i)
		fill(val, a[i], m, k);
}

template<class T>
bool equal(T **a, T **b, int n, int m) {
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			if (a[i][j] != b[i][j])
				return false;
	return true;
}

template<class T>
void rotateRight(T *a, int c1, int c2) {
	T t = a[c2];
	for (int i = c2; i > c1; --i)
		a[i] = a[i - 1];
	a[c1] = t;
}

template<class T>
void rotateDown(T **a, int col, int r1, int r2) {
	T t = a[r2][col];
	for (int i = r2; i > r1; --i)
		a[i][col] = a[i - 1][col];
	a[r1][col] = t;
}

template<class T>
int columnCompare(T **a, int nRows, int c1, int c2) {
	for (int i = 0; i < nRows; ++i)
		if (a[i][c1] != a[i][c2])
			return a[i][c1] - a[i][c2];
	return 0;
}

template<class T>
int compare(T *a, T *b, int n) {
	for (int i = 0; i < n; ++i)
		if (a[i] != b[i])
			return a[i] - b[i];
	return 0;
}
template<class T>
void sortColumns(T **a, int n, int m) {
	for (int i = 1; i < m; ++i) {
		int j = 0;
		while (j < i && columnCompare(a, n, j, i) <= 0)
			++j;
		if (j < i) {
			for (int k = 0; k < n; ++k)
				rotateRight(a[k], j, i);
		}
	}
}

template<class T>
void sortRows(T **a, int n, int m) {
	for (int i = 1; i < n; ++i) {
		int j = 0;
		while (j < i && compare(a[j], a[i], m) <= 0)
			++j;
		if (j < i)
			for (int k = 0; k < m; ++k)
				rotateDown(a, k, j, i);
	}
}

template<class T>
void printMatrix(T **a, int n, int m) {
	for (int i = 0; i < n; ++i, cerr << endl)
		for (int j = 0; j < m; ++j)
			cerr << " " << (int) a[i][j];
}

template<class T>
int maxPos(T a, int n) {
	int m = a[0];
	int p = 0;
	for (int i = 1; i < n; ++i)
		if (a[i] > m) {
			m = a[i];
			p = i;
		}
	return p;
}

template<class T>
T maxArray(T *a, int n) {
	T m = a[0];
	for (int i = 1; i < n; ++i)
		if (a[i] > m)
			m = a[i];
	return m;
}

template<class T>
T minArray(T *a, int n) {
	T m = a[0];
	for (int i = 1; i < n; ++i)
		if (a[i] < m)
			m = a[i];
	return m;
}
void *buildMultiDArray(int dimensions, int n);

void fillMultiDArray(void *a, int dimensions, int n, int val);

template<class T>
void printArray(T *a, int n, string msg) {
	std::cerr << msg << " ";
	for (int i = 0; i < n; ++i)
		std::cerr << a[i] << " ";
	std::cerr << endl;
}

template<class T>
void permute(T *a, int *perm, int n) {
	T b[n];
	for (int i = 0; i < n; ++i)
		b[i] = a[perm[i]];
	for (int i = 0; i < n; ++i)
		a[i] = b[i];
}

bool hasNegative(int *a, int n);

void insertLeft(int *a, int pos);

void insertSort(int *a, int n);

template<class T>
T findMin(T *a, int n) {
	T m = a[0];
	for (int i = 0; i < n; ++i)
		if (a[i] < m)
			m = a[i];
	return m;
}

template<class T>
T findMax(T *a, int n) {
	T m = a[0];
	for (int i = 0; i < n; ++i)
		if (a[i] > m)
			m = a[i];
	return m;
}

template<class T>
T sum(T* data, int n) {
	T total = 0;
	for (int j = 0; j < n; ++j)
		total += data[j];
	return total;
}

void lexySmallestNeighbor(char *s, int L, int d);

void printLmer(const char *x, int L);

int maxFreq(int a, int b, int c);

int sumMaxFreq(char *a, char *b, char *c, int L);

int sumMaxFreq(char *a, char *b, char *c, char *d, int L);

void encodeString(char *s, int m, string &sigma);

void encodeStrings(char **s, int n, int m, string &sigma);

void trim(string& s);

bool readFasta(const char *fileName, vector<string>& strings);

void encodeString(string& s, string &sigma);

void encodeStrings(vector<string>& strings, string &sigma);

int sumLengths(vector<string>& v);

char **pack(vector<string>& strings);

string decodeString(const char *s, int n, const string& sigma);

bool isBitSet(int64 *a, int index, int*bitMask);

void setBit(int64 *a, int index, int*bitMask);

bool isBitSet(uint32 *a, int index);

void setBit(uint32 *a, int index);

bool isBitSet(char *a, int index, int*bitMask);

void setBit(char *a, int index, int*bitMask);

void parseIntOrExit(char option, char *argument, int *result);

void parseFloatOrExit(char option, char *argument, float *result);

//bool isMotif(const char *m, PriorityMatrix *good, int L, int d, int threshold);

/*
 template<class T>
 bool isMotif(const char *m, const char *memStart, row<T> **good, int n, int L,
 int d) {
 bool exist = true;
 for (int i = 0; (i < n) & exist; ++i) {
 exist = false;
 row<T> *r = good[i];
 for (int j = 0; (j < r->nItems) & (!exist); ++j)
 exist = HamDist(m, memStart + r->itemIndex[j], L) <= d;
 }
 return exist;
 }
 */
#endif /* UTILS_H_ */
