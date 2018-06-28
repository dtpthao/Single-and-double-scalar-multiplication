#include "option_BNBC.h"

Vector::Vector(char m, char n) {
	this->m = m;
	this->n = n;
}

TripleVector::TripleVector(Vector V, int d) {
	this->d = d;
	this->V = V;
	Vector Vmod2(V.m & 1, V.n & 1);
	v[0].m = V.m + ((Vmod2.m + 1) & 1);
	v[0].n = V.n + ((Vmod2.n + 1) & 1);
	v[1] = Vector(V.m + Vmod2.m, V.n + Vmod2.n);
	v[2].m = V.m + ((V.m + d) & 1);
	v[2].n = V.n + ((V.n + d + 1) & 1);
}

TripleVector::TripleVector() {
	TripleVector(Vector(0, 0), 0);
}

int dNext(int D, Vector V, TripleVector tV) {
	
	Vector cond(V.m + tV.V.m, V.n + tV.V.n);
	cond.m &= 1;
	cond.n &= 1;
	if (cond.m == 0 && cond.n == 1) return 0;
	if (cond.m == 1 && cond.n == 0) return 1;
	if (cond.m == 0 && cond.n == 0) return D;
	return 1 - D;
}

int v2mod4(Vector v2) {
	if ((v2.m & 3) != (v2.n & 3)) return 2;
	if ((v2.m ^ 3 || v2.n ^ 3) && ((v2.m & 3) != 2)) return 1;
	return 0;
}

int v3parities(Vector v3, Vector v3next) {
	return (((v3.m ^ v3next.m) | (v3.n ^ v3next.n)) & 1) ? 0 : 1;
}

///////////////////////////////////////////////////////////////////////////////////
void print3v(Vector v[3]) {
	printf("(%8x, %8x)   (%8x, %8x)   (%8x, %8x)\n", v[0].m, v[0].n,
		v[1].m, v[1].n, v[2].m, v[2].n);
}
void printVector(Vector v) {
	printf("(%2x, %2x)", v.m, v.n);
}

void printTv(TripleVector tv) {
	std::cout << "V: "; printVector(tv.V);
	std::cout << "  d: " << tv.d << ", ";
	std::cout << "(";
	printVector(tv.v[0]);
	std::cout << ", ";
	printVector(tv.v[1]);
	std::cout << ", ";
	printVector(tv.v[2]);
	std::cout << ")   "; //<< std::endl;
}
////////////////////////////////////////////////////////////////////////////
//_________________________________________________________
DWORD bitlen(big k) {
	int i = 31;
	while (!(k->w[k->len - 1] & (1 << i))) i--;
	return (k->len << 5) - (31 - i);
}

void nextStep(int &d, Vector &V, TripleVector &tV, Vector* CS, DWORD &len) {
	
	d = dNext(d, V, tV);
	TripleVector tVnext = TripleVector(V, d);
	CS[len] = Vector(v2mod4(tV.v[1]), v3parities(tV.v[2], tVnext.v[2]));
	tV = tVnext;
	len++;
}

DWORD GenBNBC(big a, big b, Vector* CS, DiffSeq &DS) {	//without using sftbit(big, n, big);

	DWORD len = 1,
		maxlen = (compare(a, b) > 0) ? bitlen(a) : bitlen(b);
	int d = a->w[0] & 1;
	Vector V(a->w[0], b->w[0]), constV = V;
	if (a->len > 1 || b->len > 1) { V.m += 4; V.n += 4; }
	TripleVector tV(V, d);
	int iw = 0, sft = 1;
	while (len < maxlen - 2) {
		if (sft < 31) {
			V.m = (a->w[iw] >> sft) | 4;
			V.n = (b->w[iw] >> sft) | 4;
			sft++;
		}
		else {
			V.m = (a->w[iw] >> 31) | (a->w[iw + 1] << 1) | 4;
			V.n = (b->w[iw] >> 31) | (b->w[iw + 1] << 1) | 4;
			sft = 0; iw++;
		}
		nextStep(d, V, tV, CS, len);
	}
	if (sft == 31) {
		V.m = (a->w[iw] >> 31) | (a->w[iw + 1] << 1);
		V.n = (b->w[iw] >> 31) | (b->w[iw + 1] << 1);
		iw++;
	}
	else {
		V.m = (a->w[iw] >> sft) & 3;
		V.n = (b->w[iw] >> sft) & 3;
	}
	nextStep(d, V, tV, CS, len);
	sft++;
	V.m = (a->w[iw] >> sft) & 1;
	V.n = (b->w[iw] >> sft) & 1;
	nextStep(d, V, tV, CS, len);
	DS.ds1 = Vector(tV.v[0].m - tV.v[1].m, tV.v[0].n - tV.v[1].n);
	DS.ds2.m = tV.v[CS[len - 1].n].m - tV.v[2].m;
	DS.ds2.n = tV.v[CS[len - 1].n].n - tV.v[2].n;
	return len;
}

void InitialState(pepoint R[3], pepoint P, pepoint Q, DiffSeq DS, Vector CS) {
	R[0] = epoint_init();
	R[1] = epoint_init();
	R[2] = epoint_init();
	epoint2_copy(P, R[0]);
	ecurve2_add(Q, R[0]);
	if (DS.ds1.m == -1 && DS.ds1.n == -1) {
		epoint2_copy(R[0], R[1]);
	}
	else if (DS.ds1.m == 1 && DS.ds1.n == -1) { 
		epoint2_copy(Q, R[1]);
	}
	else {
		epoint2_copy(P, R[1]);
	}
	ecurve2_add(R[1], R[1]);

	epoint2_copy(R[CS.n], R[2]);
	if (DS.ds2.m == 0 && DS.ds2.n == 1) 
		ecurve2_sub(Q, R[2]);
	else if (DS.ds2.m == 0 && DS.ds2.n == -1) 
		ecurve2_add(Q, R[2]);
	else if (DS.ds2.m == 1 && DS.ds2.n == 0) 
		ecurve2_sub(P, R[2]);
	else ecurve2_add(P, R[2]);
}

void ShamirMul_BNBC(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	DWORD len = ((a->len > b->len) ? a->len : b->len) << 5;
	Vector CS[MAX_M];
	DiffSeq DS;
	len = GenBNBC(a, b, CS, DS);
	pepoint Ri[3], Rtmp = epoint_init();
	InitialState(Ri, P, Q, DS, CS[len - 1]);
	for (int i = len - 1; i > 0; i--) {
		epoint2_copy(Ri[CS[i].m], Rtmp);
		ecurve2_add(Rtmp, Rtmp);
		ecurve2_add(Ri[CS[i].n], Ri[2]);
		ecurve2_add(Ri[1], Ri[0]);
		epoint2_copy(Rtmp, Ri[1]);
	}
	int amod2 = a->w[0] & 1, bmod2 = b->w[0] & 1;
	if (amod2 != bmod2)
		epoint2_copy(Ri[2], R);
	else if (amod2 == 0)
		epoint2_copy(Ri[1], R);
	else
		epoint2_copy(Ri[0], R);
	epoint_free(Ri[0]);
	epoint_free(Ri[1]);
	epoint_free(Ri[2]);
	epoint_free(Rtmp);
}

void TestShrMulBNBC(csprng &Rng, pepoint P) {
	big a = mirvar(71), b = mirvar(1);
	pepoint Q = epoint_init();
	pepoint R = epoint_init();
	pepoint R1 = epoint_init();
	ecurve2_mult(a, P, Q);
	int count = 0;
	for (int i = 0; i < 10000; i++) {
		strong_bigdig(&Rng, 50, 10, a);
		strong_bigdig(&Rng, 50, 10, b);
		//std::cout << "a: "; cotnum(a, stdout);
		//std::cout << "b: "; cotnum(b, stdout);
		ShamirMul_BNBC(a, P, b, Q, R); //std::cout << "R: "; cotnumEp(R);
		ecurve2_mult2(a, P, b, Q, R1); //std::cout << "R3: "; cotnumEp(R3);
		count += epoint2_comp(R, R1);
		//ecurve2_sub(R, R3);
		//std::cout << "R3: "; cotnumEp(R3);
	}
	std::cout << "Cmp: " << count << std::endl;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//Vector::~Vector() {}

//Vector& Vector::operator= (const Vector &v) {
//	if (this == &v) return *this;
//	m = v.m;
//	n = v.n;
//	return *this;
//}

//bool Vector::operator== (const Vector &v) {
//	return (this->m == v.m) && (this->n == v.n);
//}
//
//Vector operator+ (Vector const &v1, Vector const &v2) {
//	return Vector(v1.m + v2.m, v1.n + v2.n);
//}
//
//Vector operator+ (Vector const &v1, int const x) {
//	return Vector(v1.m + x, v1.n + x);
//}
//
//Vector& Vector::operator+= (const Vector &v) {
//	m += v.m;
//	n += v.n;
//	return *this;
//}
//
//Vector operator- (Vector const &v1, Vector const &v2) {
//	return Vector(v1.m - v2.m, v1.n - v2.n);
//}
//
//Vector operator& (Vector const &v1, Vector const &v2) {
//	return Vector(v1.m & v2.m, v1.n & v2.n);
//}
//
//Vector operator& (Vector const &v1, int const x) {
//	return Vector(v1.m & x, v1.n & x);
//}
//
//Vector operator| (Vector const&v1, Vector const& v2) {
//	return Vector(v1.m | v2.m, v1.n | v2.n);
//}
//
//Vector operator| (Vector const&v, int const& x) {
//	return Vector(v.m | x, v.n | x);
//}
//
//Vector operator>> (Vector const&v, int const& x) {
//	return Vector(v.m >> x, v.n >> x);
//}
//
//Vector operator<< (Vector const&v, int const& x) {
//	return Vector(v.m << x, v.n << x);
//}
//
//DWORD GenBNBC2(big a, big b, Vector* CS, DiffSeq &DS) {	//using sftbit(big, n, big);
//	big a2 = mirvar(1), b2 = mirvar(1);
//	copy(a, a2); copy(b, b2);
//	DWORD len = 1;
//	int d = a2->w[0] & 1;
//	Vector V(a2->w[0], b2->w[0]);
//	if (a->len > 1 || b->len > 1) V += 4;
//	TripleVector tV(V, d), tVnext;
//	sftbit(a2, -1, a2);
//	sftbit(b2, -1, b2);
//	while (a2->len != 0 || b2->len != 0) {
//		printf("%2d   ", len); printTv(tV);							////////
//		V = Vector(a2->w[0], b2->w[0]);
//		if (a->len > 1 || b->len > 1) V += 4;
//		d = dNext(d, V, tV);
//		tVnext = TripleVector(V, d);
//		CS[len] = Vector(v2mod4(tV.v[1]), v3parities(tV.v[2], tVnext.v[2]));
//		printVector(CS[len]); std::cout << std::endl;//////////////////////
//		tV = tVnext;
//		len++;
//		sftbit(a2, -1, a2);
//		sftbit(b2, -1, b2);
//	}
//	DS.ds1 = tV.v[0] - tV.v[1];
//	DS.ds2 = tV.v[CS[len - 1].n] - tV.v[2];
//	printf("%2d   ", len); printTv(tV); std::cout << std::endl;//////////////
//	printVector(DS.ds1); printVector(DS.ds2);
//	mirkill(a2); mirkill(b2);
//	return len;
//}
//
//DWORD GenBNBC3(big a, big b, Vector* CS, DiffSeq &DS) {		//without using sftbit(big, n, big);
//
//	DWORD len = 1,
//		maxlen = (compare(a, b) > 0) ? bitlen(a) : bitlen(b);
//	int d = a->w[0] & 1;
//	Vector V(a->w[0], b->w[0]), constV = V;
//	if (a->len > 1 || b->len > 1) V += 4;
//	TripleVector tV(V, d), tVnext;
//	int iw = 0, sft = 1;
//	while (len < maxlen - 2) {
//		if (sft < 31) {
//			V = (constV >> sft) | 4;
//			sft++;
//		}
//		else {
//			V = (constV >> 31) & 1;
//			constV = Vector(a->w[++iw], b->w[iw]);
//			V = V | (constV << 1) | 4;
//			sft = 0;
//		}
//		nextStep(d, V, tV, CS, len);
//	}
//	if (sft == 31) {
//		V = (constV >> 31) & 1;
//		constV = Vector(a->w[++iw], b->w[iw]);
//		V = V | (constV << 1);
//		sft = 0;
//	}
//	else {
//		V = (constV >> sft) & 3;
//		sft++;
//	}
//	nextStep(d, V, tV, CS, len);
//
//	V = (constV >> sft) & 1;
//	nextStep(d, V, tV, CS, len);
//
//	DS.ds1 = tV.v[0] - tV.v[1];
//	DS.ds2 = tV.v[CS[len - 1].n] - tV.v[2];
//	return len;
//}

//Vector CheckGenBinChain(Vector* CS, DiffSeq DS, int len, big a, big b) {
//	Vector V[3];
//	V[0] = Vector(1, 1);
//	V[1] = V[0] - DS.ds1;
//	V[2] = V[CS[len - 1].n] - DS.ds2;
//	/*printVector(DS.ds1); std::cout << ", ";
//	printVector(DS.ds2); std::cout << std::endl;
//	printf("%2d: ", len - 1); std::cout << "CS: ";
//	printVector(Vector(0, 0)); std::cout << "   ";
//	print3v(V);*/
//	Vector Vtmp;
//	for (int i = len - 1; i > 0; i--) {
//		/*printf("%2d: ", i); std::cout << "CS: ";
//		printVector(CS[i]); std::cout << "   ";*/
//		Vtmp = V[CS[i].m] + V[CS[i].m];
//		V[2] += V[CS[i].n];
//		V[0] += V[1];
//		V[1] = Vtmp;
//		//print3v(V);
//	}
//	int amod2 = a->w[0] & 1, bmod2 = b->w[0] & 1;
//	if (amod2 != bmod2) return V[2];
//	else if (amod2 == 0) return V[1];
//	else return V[0];
//}
//
//void TestGenBinChain(csprng &Rng) {
//	big a = mirvar(1), b = mirvar(1);
//	DWORD len;
//	Vector* CS = new Vector[200];
//	//DiffSeq* DS = new DiffSeq[200];
//	DiffSeq DS;
//	for (int i = 0; i < 1; i++) {
//		strong_bigdig(&Rng, 9, 16, a);
//		strong_bigdig(&Rng, 9, 16, b);
//		//a->w[0] = 71;// 0x83701758;
//		//a->len = 2; a->w[1] = 1;
//		//b->len = 2; b->w[0] = 93;// 0xE2A082C8;
//		//b->w[1] = 1;
//		len = GenBNBC(a, b, CS, DS);
//		std::cout << std::endl << std::endl;
//		Vector v = CheckGenBinChain(CS, DS, len, a, b);
//
//		printf("(%8x, %8x)\n", v.m, v.n);
//		std::cout << "a: " << a->len << "  "; cotnum(a, stdout);
//		std::cout << "b: " << b->len << "  "; cotnum(b, stdout);
//		
//		/*std::cout << std::endl << "~~~~~~~~" << std::endl;
//		len = GenBNBC2(a, b, CS, DS);
//		std::cout << std::endl << std::endl;
//		v = CheckGenBinChain(CS, DS, len, a, b);
//
//		printf("(%8x, %8x)\n", v.m, v.n);
//		std::cout << "a: " << a->len << "  "; cotnum(a, stdout);
//		std::cout << "b: " << b->len << "  "; cotnum(b, stdout);*/
//	}
//}
