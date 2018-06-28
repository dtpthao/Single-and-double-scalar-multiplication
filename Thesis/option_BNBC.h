#ifndef _OPTION_BNBC_H
#define _OPTION_BNBC_H

#include <iostream>
#include "Ellipse.h"
#include "Option_NAF.h"

struct Vector {
	char m;
	char n;
	Vector(char m = 0, char n = 0);
};

struct TripleVector {
public:
	Vector V;
	int d;
	Vector v[3];
	TripleVector();
	TripleVector(Vector V, int d);
};

typedef struct DifferentSequence {
	Vector ds1;
	Vector ds2;
} DiffSeq;

int dNext(int D, Vector V, TripleVector tV);

int v2mod4(Vector v2);

int v3parities(Vector v3, Vector v3next);

DWORD GenBNBC(big a, big b, Vector* CS, DiffSeq &DS);

void InitialState(pepoint PQ[3], pepoint P, pepoint Q, DiffSeq DS, Vector CS);

void ShamirMul_BNBC(big, pepoint, big, pepoint, pepoint);

void TestGenBinChain(csprng &);

void TestShrMulBNBC(csprng &, pepoint);

//class Vector {
//public:
//	int m;
//	int n;
//	Vector(int m = 0, int n = 0);
//	/*~Vector();
//	Vector& operator= (const Vector &v);
//	Vector& operator+= (const Vector &v);
//	bool operator== (const Vector &v);
//	friend Vector operator+ (Vector const&, Vector const&);
//	friend Vector operator+ (Vector const&, int const);
//	friend Vector operator& (Vector const&, Vector const&);
//	friend Vector operator& (Vector const&, int const);
//	friend Vector operator- (Vector const&, Vector const&);
//	friend Vector operator>> (Vector const&, int const&);*/
//};

//struct Vector {
//	int m;
//	int n;
//	Vector(int m = 0, int n = 0);
//};

#endif
