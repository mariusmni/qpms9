/*
 * MyString.cpp
 *
 *  Created on: Jul 4, 2012
 *      Author: marius
 */

#include "MyString.h"
#include <cstring>
#include <algorithm>

MyString::MyString(const MyString& other) :
		L(other.L) {
	s = new char[L];
	memcpy(s, other.s, L * sizeof(char));
}

MyString::MyString(const char *st, int L) :
		L(L) {
	s = new char[L];
	memcpy(s, st, L * sizeof(char));
}

MyString::MyString(const int *st, int L) :
		L(L) {
	s = new char[L];
	for (int i = 0; i < L; ++i)
		s[i] = (char) st[i];
}

MyString::~MyString() {
	delete[] s;
}

bool MyString::operator<(const MyString& b) const {
	for (int i = 0; i < L; ++i)
		if (s[i] != b.s[i])
			return s[i] < b.s[i];
	return false;
}

bool MyString::operator==(const MyString& b) const {
	if (L != b.L)
		return false;
	for (int i = 0; i < L; ++i)
		if (s[i] != b.s[i])
			return false;
	return true;
}

MyString& MyString::operator=(MyString rhs) {
	std::swap(L, rhs.L);
	std::swap(s, rhs.s);
	return *this;
}
