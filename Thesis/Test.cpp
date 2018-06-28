#include "Test.h"


void TestAllScalar(big k, pepoint P, csprng &Rng, big n,
	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB])
{
	pepoint R = epoint_init();
	pepoint R1 = epoint_init();
	for (int i = 0; i <= LIB; i++) {
		tmin[i] = 100;
	}
	for (int i = 0; i < TESTALL; i++) {
		strong_bigrand(&Rng, n, k);
		//strong_bigdig(&Rng, 7, 16, k);
		SclDuration(k, P, R, ecurve2_mult, tmin[LIB]);

		SclDuration(k, P, R1, ScalarMul_Bin, tmin[SclBin]);
		correct[SclBin] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_Bin, tmin[ShrBin]);
		correct[ShrBin] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_JSF, tmin[ShrJSF]);
		correct[ShrJSF] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_NAF, tmin[SclNAF]);
		correct[SclNAF] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_NAF, tmin[ShrNAF]);
		correct[ShrNAF] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_wNAF, tmin[SclwNAF]);
		correct[SclwNAF] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_wNAF, tmin[ShrwNAF]);//*********
		correct[ShrwNAF] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_JT, tmin[SclJT]);
		correct[SclJT] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_JT, tmin[ShrJT]);
		correct[ShrJT] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_BNBC, tmin[ShrBNBC]);
		correct[ShrBNBC] += epoint2_comp(R1, R);

		ShrDuration(k, P, R1, ShamirMul_AKDAC, tmin[ShrAKDAC]);
		correct[ShrAKDAC] += epoint2_comp(R1, R);
	}
	for (int i = 0; i <= LIB; i++) {
		per[i] = (tmin[i] / tmin[0]) * 100;
	}
	epoint_free(R); epoint_free(R1);
}

void TestScalar(big k, pepoint P, csprng &Rng, big n,
	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB])
{
	pepoint R = epoint_init(),
			R1 = epoint_init();
	for (int i = 0; i <= LIB; i++) {
		tmin[i] = 100;
	}
	for (int i = 0; i < TESTSCL; i++) {
		strong_bigrand(&Rng, n, k);
		//strong_bigdig(&Rng, 7, 16, k);
		SclDuration(k, P, R, ecurve2_mult, tmin[LIB]);

		SclDuration(k, P, R1, ScalarMul_Bin, tmin[SclBin]);
		correct[SclBin] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_NAF, tmin[SclNAF]);
		correct[SclNAF] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_wNAF, tmin[SclwNAF]);
		correct[SclwNAF] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_JT, tmin[SclJT]);
		correct[SclJT] += epoint2_comp(R1, R);
	}
	for (int i = 0; i <= LIB; i++) {
		per[i] = (tmin[i] / tmin[0]) * 100;
	}
	epoint_free(R); epoint_free(R1);
}

void TestScalarwNAF(big k, pepoint P, csprng &Rng, big n,
	double tmin[7], double per[7], unsigned int correct[7])
{
	pepoint R = epoint_init(),
		R1 = epoint_init();
	for (int i = 0; i <= 6; i++) tmin[i] = 100;
	for (int i = 0; i < TESTwNAF; i++) {
		strong_bigrand(&Rng, n, k);
		SclDuration(k, P, R, ecurve2_mult, tmin[6]);

		SclDuration(k, P, R1, ScalarMul_Bin, tmin[SclBin]);
		correct[SclBin] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_w3NAF, tmin[1]);
		correct[1] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_w4NAF, tmin[2]);
		correct[2] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_w5NAF, tmin[3]);
		correct[3] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_w6NAF, tmin[4]);
		correct[4] += epoint2_comp(R1, R);

		SclDuration(k, P, R1, ScalarMul_w7NAF, tmin[5]);
		correct[5] += epoint2_comp(R1, R);
	}
	for (int i = 0; i <= 6; i++) {
		per[i] = (tmin[i] / tmin[0]) * 100;
	}
	epoint_free(R); epoint_free(R1);
}

void TestShrwNAF(big k, pepoint P, csprng &Rng, big n,
	double tmin[7], double per[7], unsigned int correct[7])
{
	pepoint R = epoint_init(),
		R1 = epoint_init(),
		Q = epoint_init();
	big a = mirvar(1);
	big b = mirvar(1);
	for (int i = 0; i <= 6; i++) tmin[i] = 100;
	for (int i = 0; i < TESTwNAF; i++) {
		strong_bigrand(&Rng, n, k);
		strong_bigrand(&Rng, n, a);
		strong_bigrand(&Rng, n, b);
		ecurve2_mult(k, P, Q);
		ShrDuration(a, P, b, Q, R, ecurve2_mult2, tmin[6]);
		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_Bin, tmin[0]);
		correct[0] += epoint2_comp(R1, R);
		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_w3NAF, tmin[1]);//********
		correct[1] += epoint2_comp(R1, R);

		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_w4NAF, tmin[2]);//********
		correct[2] += epoint2_comp(R1, R);

		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_w5NAF, tmin[3]);//********
		correct[3] += epoint2_comp(R1, R);

		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_w6NAF, tmin[4]);//********
		correct[4] += epoint2_comp(R1, R);

		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_w7NAF, tmin[5]);//********
		correct[5] += epoint2_comp(R1, R);
	}
	for (int i = 0; i <= 6; i++) {
		per[i] = (tmin[i] / tmin[0]) * 100;
	}
	epoint_free(R); epoint_free(R1);
}

void TestShamir(big k, pepoint P, csprng &Rng, big n,
	double tmin[LIB + 1], double per[LIB + 1], unsigned int correct[LIB]) 
{
	pepoint R = epoint_init(),
		R1 = epoint_init(),
		Q = epoint_init();
	big a = mirvar(1), b = mirvar(1);
	for (int i = 0; i <= LIB; i++) tmin[i] = 100;
	
	for (int i = 0; i < TESTSHR; i++) {
		strong_bigrand(&Rng, n, k);
		strong_bigrand(&Rng, n, a);
		strong_bigrand(&Rng, n, b);
		//ShamirMul(k, P, a, Q, b);
		ecurve2_mult(k, P, Q);
		ShrDuration(a, P, b, Q, R, ecurve2_mult2, tmin[LIB]);
		/*printf("i: %d\t", i);
		printf("bin\t");*/
		ShrDuration(a, P, b, Q, R1, ShamirMul_Bin, tmin[ShrBin]);
		correct[ShrBin] += epoint2_comp(R1, R);

		//printf("JSF\t");
		ShrDuration(a, P, b, Q, R1, ShamirMul_JSF, tmin[ShrJSF]);
		correct[ShrJSF] += epoint2_comp(R1, R);

		//printf("NAF\t");
		epoint2_set(0, 0, 1, R1);
		ShrDuration(a, P, b, Q, R1, ShamirMul_NAF, tmin[ShrNAF]);
		correct[ShrNAF] += epoint2_comp(R1, R);

		//printf("wNAF\t");
		ShrDuration(a, P, b, Q, R1, ShamirMul_wNAF, tmin[ShrwNAF]);//*******
		correct[ShrwNAF] += epoint2_comp(R1, R);

		//printf("JT\t");
		ShrDuration(a, P, b, Q, R1, ShamirMul_JT, tmin[ShrJT]);
		correct[ShrJT] += epoint2_comp(R1, R);

		//printf("BNBC\t");
		ShrDuration(a, P, b, Q, R1, ShamirMul_BNBC, tmin[ShrBNBC]);
		correct[ShrBNBC] += epoint2_comp(R1, R);

		//printf("AKDAC\n");
		ShrDuration(a, P, b, Q, R1, ShamirMul_AKDAC, tmin[ShrAKDAC]);
		correct[ShrAKDAC] += epoint2_comp(R1, R);
	}
	for (int i = ShrBin; i <= LIB; i++) {
		per[i] = (tmin[i] / tmin[ShrBin]) * 100;
	}

	epoint_free(Q);
	mirkill(a); mirkill(b);
	epoint_free(R); epoint_free(R1);
}

void TestGenForm(big k, csprng &Rng, big n, double tmin[LIB + 1], DWORD tick[LIB + 1]) {

	big a = mirvar(1), b = mirvar(1);
	double dur;
	DWORD qPart;
	stopWatch timer;
	Vector* CS;
	DiffSeq DS;
	char *JSFa = new char[(MAX_M + 1) / 2];
	char *JSFb = new char[(MAX_M + 1) / 2];
	for (int i = 0; i <= LIB; i++) {
		tmin[i] = 100; tick[i] = 1000;
	}
	for (int i = 0; i < TESTGF; i++) {
		strong_bigrand(&Rng, n, k);
		ShamirMul(k, a, b);
		GenFormSclDuration(k, GenNAF, tmin[SclNAF], tick[SclNAF]);
		GenFormSclDuration(k, GenwNAF, tmin[SclwNAF], tick[SclwNAF]);
		GenFormSclDuration(k, GenJT, tmin[SclJT], tick[SclJT]);
		
		GenFormShrDuration(a, b, GenNAF, tmin[ShrNAF], tick[ShrNAF]);
		GenFormShrDuration(a, b, GenwNAF, tmin[ShrwNAF], tick[ShrwNAF]);
		GenFormShrDuration(a, b, GenJT, tmin[ShrJT], tick[ShrJT]);

		startTimer(&timer);
		GenJSF(a, b,JSFa, JSFb);
		stopTimer(&timer);
		dur = getElapsedTime(&timer) * 1000;
		tmin[ShrJSF] = (tmin[ShrJSF] < dur) ? tmin[ShrJSF] : dur;
		qPart = getQuadPart(&timer);
		tick[ShrJSF] = (tick[ShrJSF] < qPart) ? tick[ShrJSF] : qPart;

		CS = new Vector[k->len <<4];
		startTimer(&timer);
		GenBNBC(a, b, CS, DS);
		stopTimer(&timer);
		dur = getElapsedTime(&timer) * 1000;
		tmin[ShrBNBC] = (tmin[ShrBNBC] < dur) ? tmin[ShrBNBC] : dur;
		qPart = getQuadPart(&timer);
		tick[ShrBNBC] = (tick[ShrBNBC] < qPart) ? tick[ShrBNBC] : qPart;
	}
	mirkill(a); mirkill(b);
}