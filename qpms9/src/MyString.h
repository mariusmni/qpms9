/*
 * MyString.h
 *
 *  Created on: Jul 4, 2012
 *      Author: marius
 */

#ifndef MYSTRING_H_
#define MYSTRING_H_

class MyString {
public:

	MyString(const MyString& other);

	MyString(const char *s, int L);

	MyString(const int *st, int L);

	virtual ~MyString();

	bool operator<(const MyString& b) const;

	bool operator==(const MyString& b) const;

	MyString& operator=(MyString rhs);

	char *s;
	int L;
private:
};

#endif /* MYSTRING_H_ */
