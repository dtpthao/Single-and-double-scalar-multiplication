#include <Windows.h>
#include "option_NAF.h"

DWORD GenNAF(big k, char *rNAF)
{
	DWORD i = 0;
	big d = mirvar(1);
	copy(k, d);
	while (d->len != 1 || d->w[0] != 1) {
		rNAF[i] = (2 - (d->w[0] & 3)) %2;
		if (rNAF[i] == -1) {
			incr(d, 1, d);
		}
		sftbit(d, -1, d);
		i++;
	}
	rNAF[i] = 1;
	mirkill(d);
	return ++i;
}

//ScalarMul
void ScalarMul_NAF(big k, epoint* P, epoint* R)
{
	char SBR[600];
	DWORD lenSBR;
	pepoint plist[3];
	plist[0] = epoint_init();
	plist[1] = epoint_init();
	plist[2] = epoint_init();
	epoint2_copy(P, plist[2]);
	ecurve2_sub(P, plist[0]);
	epoint2_norm(plist[0]);
	epoint2_norm(plist[2]);
	lenSBR = GenNAF(k, SBR);
	epoint2_copy(P, R);
	for (int i = lenSBR - 2; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(plist[SBR[i] + 1], R);
	}
}

//ShamirMul
void PreMul(pepoint P, pepoint Q, pepoint *listpoint)
{
	for (int i = 0; i < 9; i++) {
		listpoint[i] = epoint_init();
	}
	epoint2_set(0, 0, 1, listpoint[4]);
	epoint2_copy(P, listpoint[0]); ecurve2_add(Q, listpoint[0]);
	epoint2_copy(P, listpoint[2]); ecurve2_sub(Q, listpoint[2]);
	epoint2_copy(Q, listpoint[6]); ecurve2_sub(P, listpoint[6]);
	ecurve2_sub(listpoint[0], listpoint[8]);
	epoint2_copy(P, listpoint[1]);
	epoint2_copy(Q, listpoint[3]);
	ecurve2_sub(Q, listpoint[5]);
	ecurve2_sub(P, listpoint[7]);
	for (int i = 0; i < 9; i++) {
		epoint2_norm(listpoint[i]);
	}
}

void ShamirMul_NAF(big a, pepoint P, big b, pepoint Q, pepoint R)   //SBR = NAF
{
	pepoint listpoint[9];
	char NAFa[MAX_M + 1] = { 0 };
	char NAFb[MAX_M + 1] = { 0 };
	DWORD lenNAFa, lenNAFb, lenNAF;
	PreMul(P, Q, listpoint);
	lenNAFa = GenNAF(a, NAFa);
	lenNAFb = GenNAF(b, NAFb);
	lenNAF = max(lenNAFa, lenNAFb);
	epoint2_set(0, 0, 1, R);
	for (int i = lenNAF - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listpoint[4 - 3 * NAFa[i] - NAFb[i]], R);
	}
	for (int i = 0; i < 9; i++) epoint_free(listpoint[i]);
}

DWORD GenNAF2(big k, char *rNAF)
{
	DWORD i = 0;
	big d = mirvar(1);
	copy(k, d);
	while (d->len != 1 || d->w[0] != 1) {
		if (d->w[0] & 1) {
			if ((d->w[0] & 3) == 3) {
				rNAF[i] = -1;
				incr(d, 1, d);
			}
			else rNAF[i] = 1;
		}
		else rNAF[i] = 0;
		sftbit(d, -1, d);
		i++;
	}
	rNAF[i] = 1;
	mirkill(d);
	return ++i;
}

DWORD GenNAF1(big k, char *rNAF)
{
	DWORD i = 0;
	big d = mirvar(1);
	copy(k, d);
	while (d->len != 1 || d->w[0] != 1) {
		if ((d->w[0] & 3) == 3) {
			rNAF[i] = -1;
			incr(d, 1, d);
		}
		else if ((d->w[0] & 3) == 1) {
			rNAF[i] = 1;
		}
		else rNAF[i] = 0;
		sftbit(d, -1, d);
		i++;
	}
	rNAF[i] = 1;
	mirkill(d);
	return ++i;
}