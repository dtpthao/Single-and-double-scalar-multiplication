#include "option_JT.h"

DWORD GenJT(big k, char* rJT) {
	DWORD i = 0, j = 0, len = 0;
	while (i < k->len - 1) {
		for (j = 0; j < 30; j += 2) {
			rJT[len] = (k->w[i] >> j) & 7;
			rJT[len] += (len != 0 && rJT[len - 1] < 0) ? 1 : 0;
			rJT[len++] -= 4;
		}
		rJT[len] = (k->w[i] >> 30) & 3;
		rJT[len] += (k->w[i + 1] & 1) << 2;
		rJT[len] += (rJT[len - 1] < 0) ? 1 : 0;
		rJT[len++] -= 4;
		i++;
	}
	j = 0;
	while ((k->w[i] >> j) && j < 32) {
		rJT[len] = (k->w[i] >> j) & 7;
		rJT[len] += (len != 0 && rJT[len - 1] < 0) ? 1 : 0;
		rJT[len++] -= 4;
		j += 2;
	}
	if (len != 0) rJT[len - 1] += 4;
	return len;
}

int checkGenJT(char* rJT, DWORD len, big k) {

	big kCheck = mirvar(0);
	for (int i = len - 1; i >= 0; i--) {
		sftbit(kCheck, 2, kCheck); 
		incr(kCheck, rJT[i], kCheck);
	}
	return compare(k, kCheck);
}

void PreMulScl_JT(pepoint P, pepoint *listP) //{-4P, -3P, -2P, -P, , P, 2P, 3P}
{
	for (int i = 0; i < 8; i++) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[5]);			//listP[5] = P
	epoint2_copy(P, listP[6]);
	ecurve2_add(listP[6], listP[6]);	//listP[6] = 2P
	epoint2_copy(listP[6], listP[7]);
	ecurve2_add(P, listP[7]);			//listP[7] = 3P
	ecurve2_sub(listP[5], listP[3]);	//listP[3] = -P
	ecurve2_sub(listP[6], listP[2]);	//listP[2] = -2P
	ecurve2_sub(listP[7], listP[1]);	//listP[1] = -3P
	epoint2_copy(listP[2], listP[0]);
	ecurve2_add(listP[0], listP[0]);	//listP[0] = -4P

	for (int i = 0; i < 8; i++) {
		epoint2_norm(listP[i]);
	}
}

void ScalarMul_JT(big k, pepoint P, pepoint R)
{
	char rJT[300];
	DWORD lenJT;
	pepoint listP[8];
	for (int i = 0; i < 8; i++) {
		listP[i] = epoint_init();
	}
	lenJT = GenJT(k, rJT);
	epoint2_set(0, 0, 1, R);
	PreMulScl_JT(P, listP);
	for (int i = lenJT - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(R, R);
		ecurve2_add(listP[rJT[i] + 4], R);
	}
	for (int i = 0; i < 8; i++) {
		epoint_free(listP[i]);
	}
}

void PreMul_JT_PQ(pepoint P, pepoint listP[8], 
	pepoint Q, pepoint listQ[8], pepoint listPQ[16]) 
{
	for (int i = 0; i < 8; i++) {
		listPQ[i] = epoint_init();
		listPQ[i + 8] = epoint_init();
	}
	int index;
	for (int i = -3; i < 4; i += 2) {
		for (int j = -3; j < 4; j += 2) {
			index = ((i + 3) << 1) + ((j + 3) >> 1);
			epoint2_copy(listP[i + 4], listPQ[index]);
			ecurve2_add(listQ[j + 4], listPQ[index]);
			epoint2_norm(listPQ[index]);

		}
	}
}

void fourfold_R(pepoint R) {
	ecurve2_add(R, R);
	ecurve2_add(R, R);
}

void ShrMulJT_expStep(pepoint R, pepoint* listPQ, 
	char* rJT, DWORD begin, DWORD end) 
{
	for (DWORD i = begin - 1; i >= end; i--) {
		fourfold_R(R);
		ecurve2_add(listPQ[rJT[i] + 4], R);
	}
}

void freeListP(pepoint* listP, pepoint* listQ, pepoint* listPQ) {
	for (int i = 0; i < 8; i++) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
		epoint_free(listPQ[i]);
		epoint_free(listPQ[i + 8]);
	}
}

void ShamirMul_JT(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char rJTa[300] = { 0 },
		 rJTb[300] = { 0 };
	DWORD lenJTa, lenJTb;
	pepoint listP[8], listQ[8], listPQ[16];
	lenJTa = GenJT(a, rJTa);
	lenJTb = GenJT(b, rJTb);
	PreMulScl_JT(P, listP);
	PreMulScl_JT(Q, listQ);
	PreMul_JT_PQ(P, listP, Q, listQ, listPQ);

	epoint2_set(0, 0, 1, R);
	DWORD len = lenJTa >= lenJTb ? lenJTb : lenJTa;
	if (lenJTa > lenJTb) {
		ShrMulJT_expStep(R, listP, rJTa, lenJTa, len);
	}
	else {
		ShrMulJT_expStep(R, listQ, rJTb, lenJTb, len);
	}
	int index;
	for (int i = len - 1; i > 0; i--) {
		fourfold_R(R);
		index = ((rJTa[i] + 3)<<1) + ((rJTb[i] + 3)>>1);
		ecurve2_add(listPQ[index], R);
	}
	fourfold_R(R);
	ecurve2_add(listP[rJTa[0] + 4], R);
	ecurve2_add(listQ[rJTb[0] + 4], R);
	freeListP(listP, listQ, listPQ);
}

void Test_SclMul_JT(csprng &Rng, pepoint P) {
	big k = mirvar(1);
	pepoint R = epoint_init();
	pepoint R1 = epoint_init();
	DWORD count = 0;
	for (int i = 0; i < 1000; i++) {
		strong_bigdig(&Rng, 10, 10, k);
		//cout << "k: "; cotnum(k, stdout);
	/*	len = GenJT(k, rJT);
		count += !checkGenJT(rJT, len, k);*/
		ecurve2_mult(k, P, R);
		
		ScalarMul_JT(k, P, R1);
	
		count += epoint2_comp(R1, R);
		if (!epoint2_comp(R1, R)) {
			cout << "R: " << endl;
			cotnumEp(R);
			cout << "R1: " << endl;
			cotnumEp(R1);
			break;
		}
		//cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}
}

void Test_ShrMul_JT(csprng &Rng, pepoint P) {
	pepoint Q = epoint_init(),
			R = epoint_init(),
			R1 = epoint_init();
	double dur1 = 10, dur2 = 10, dur;
	big a = mirvar(3), b = mirvar(2);
	ecurve2_mult(a, P, Q);
	cout << "P: " << endl;
	cotnumEp(P);
	cout << "Q: " << endl;
	cotnumEp(Q);
	DWORD len, count = 0, count1 = 0;
	stopWatch timer;
	for (int i = 0; i < 10000; i++) {
		strong_bigdig(&Rng, 20, 10, a);
		strong_bigdig(&Rng, 20, 10, b);
		/*a->w[0] = 0x5f;
		b->w[0] = 0x14;*/
		////if (i == 0) {
		//	cout << "a: "; cotnum(a, stdout);
		//	cout << "b: "; cotnum(b, stdout);
		////}
		ecurve2_mult2(a, P, b, Q, R);

		startTimer(&timer);
		ShamirMul_JT(a, P, b, Q, R1);
		stopTimer(&timer);
		dur = getElapsedTime(&timer) * 1000;
		dur1 = (dur > dur1) ? dur1 : dur;
		count += epoint2_comp(R1, R);
		if (!epoint2_comp(R1, R)) {
			cout << "Variant 1" << endl;
			cout << "a: "; cotnum(a, stdout);
			cout << "b: "; cotnum(b, stdout);
			cout << "P: " << endl;
			cotnumEp(P);
			cout << "Q: " << endl;
			cotnumEp(Q);
			cout << "R: " << endl;
			cotnumEp(R);
			cout << "R1: " << endl;
			cotnumEp(R1);
			//ShamirMul_JT(a, P, b, Q, R1);
			break;
		}

		/*startTimer(&timer);
		ShamirMul_JT2(a, P, b, Q, R1);
		stopTimer(&timer);
		dur = getElapsedTime(&timer) * 1000;
		dur2 = (dur > dur2) ? dur2 : dur;

		count1 += epoint2_comp(R1, R);
		if (!epoint2_comp(R1, R)) {
			cout << "R: " << endl;
			cotnumEp(R);
			cout << "R1: " << endl;
			cotnumEp(R1);
			break;
		}*/
		//cout << "______________________________________________" << endl;
	}

	cout << "dur_ShamirMul_JT: " << dur1 << endl;
	cout << "dur_ShamirMul_JT2: " << dur2 << endl;
	//Enum(P, R1);
	cout << "Correct: " << count << "/10000" << endl;
	cout << "Correct1: " << count1 << "/10000" << endl;
}
