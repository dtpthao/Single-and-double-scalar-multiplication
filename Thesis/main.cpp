#include <iostream>
#include <math.h>
#include <iomanip>

extern "C" {
#include "miracl.h"
	FILE _iob[] = { *stdin, *stdout, *stderr };
	extern "C" FILE * __cdecl __iob_func(void) { return _iob; }
}
using namespace std;

#include "Ellipse.h"
#include "Test.h"
#include "PrintTest.h"

int main()
{
	srand(time(NULL));
	miracl *M = mirsys(100, 0);
	M->IOBASE = 16;
	csprng Rng; InitStrongRNG(&Rng);
	big a = mirvar(1), x = mirvar(1), y = mirvar(1), n = mirvar(1);
	big b = mirvar(0xe);
	EC_CONSTANTS_F2m_POLY EC = {};

	pepoint P = epoint_init();
	big k = mirvar(1);
	int m[5] = { 163,233,283,409,571 };
	double tmin[6][LIB + 1];
	double per[6][LIB + 1] = { 0 };
	unsigned int corr[6][LIB] = { 0 };

	double tminShr[6][LIB + 1];
	double perShr[6][LIB + 1] = { 0 };
	unsigned int corrShr[6][LIB] = { 0 };

	double tminScl[6][LIB + 1];
	double perScl[6][LIB + 1] = { 0 };
	unsigned int corrScl[6][LIB] = { 0 };

	double tminwNAF[6][7];
	double perwNAF[6][7] = { 0 };
	unsigned int corrwNAF[6][7] = { 0 };

	double tminGenForm[5][LIB + 1];
	DWORD tickGF[5][LIB + 1] = { 0 };

	for (int i = 0; i < 5; i++) {
		GetConstainsEC(EC, m[i]);
		if (!GenEC(EC, a, b, P, x, y, n)) return 1;
		TestGenForm(k, Rng, n, tminGenForm[i], tickGF[i]);
		TestAllScalar(k, P, Rng, n, tmin[i], per[i], corr[i]);
		TestShamir(k, P, Rng, n, tminShr[i], perShr[i], corrShr[i]);
		TestScalarwNAF(k, P, Rng, n, tminwNAF[i], perwNAF[i], corrwNAF[i]);
	}
	//cout << "wNAF: " << wNAF << endl;
	printGenFormScl(m, tminGenForm, tickGF);
	printGenFormShr(m, tminGenForm, tickGF);
	Line; Line;
	printAll_m(m, tmin, per, corr);
	Line;Line;
	printAll_Scl(m, tmin, per, corr);
	Line;Line;
	printAll_Shr(m, tminShr, perShr, corrShr);
	Line; Line;
	printAll_SclwNAF(m, tminwNAF, perwNAF, corrwNAF);

	/*GetConstainsEC(EC, 163);
	if (!GenEC(EC, a, b, P, x, y, n)) return 1;
	TestShrMulBNBC(Rng, P);*/
	//TestShrMulAKDAC(Rng, P);
	//Test_ShrMul_JT(Rng, P);
	//TestGenBinChain(Rng);

	epoint_free(P);
	mirkill(a); mirkill(b); mirkill(k);
	mirexit();
	system("pause");
	return 0;
}

//__________________________________________________

