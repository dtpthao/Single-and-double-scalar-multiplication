#include "option_wNAF.h"

DWORD GenwNAF(big k, char *SDR)
{
	int eSDR[Pow2wNAF];
	for (int i = 1; i < (Pow2wNAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2wNAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + (y << 5)), d);
		SDR[lenSDR++] = eSDR[d->w[0] & MAXeSDR];
		sftbit(d, -wNAF, d);
		if (SDR[lenSDR - 1] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < wNAF; i++) SDR[lenSDR++] = 0;
	}
	mirkill(d);
	return lenSDR;
}

DWORD GenwNAF(big k, char *SDR, char *eSDRout, DWORD &lene)
{
	int eSDR[Pow2wNAF];
	int tmp[Pow2wNAF] = { 0 };
	for (int i = 1; i < (Pow2wNAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2wNAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	lene = 0;
	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + y * 32), d);
		SDR[lenSDR] = eSDR[d->w[0] & MAXeSDR];
		tmp[d->w[0] & MAXeSDR]++;
		sftbit(d, -wNAF, d);
		if (SDR[lenSDR++] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < wNAF; i++) SDR[lenSDR++] = 0;
	}
	for (int i = HalfPow2wNAF - 1; i > 0; i -= 2) {
		if (tmp[i] || tmp[Pow2wNAF - i]) {
			eSDRout[lene++] = eSDR[i];
		}
	}
	mirkill(d);
	return lenSDR;
}

//////////////////////////////////////////////////////
void PreMul_wNAF(pepoint P, pepoint *listP)
{
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	for (int i = 1; i < Pow2wNAF; i += 2) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[HalfPow2wNAF + 1]);
	ecurve2_sub(listP[HalfPow2wNAF + 1], listP[HalfPow2wNAF - 1]);
	for (int i = 3; i < HalfPow2wNAF; i += 2) {
		epoint2_copy(listP[HalfPow2wNAF + i - 2], listP[HalfPow2wNAF + i]);
		ecurve2_add(P2, listP[HalfPow2wNAF + i]);
		ecurve2_sub(listP[HalfPow2wNAF + i], listP[HalfPow2wNAF - i]);
	}
	listP[HalfPow2wNAF] = epoint_init();
	epoint_free(P2);
}

void ScalarMul_wNAF(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[Pow2wNAF];
	lenSDR = GenwNAF(k, SDR);
	PreMul_wNAF(P, listP);
	epoint2_set(0, 0, 1, R);
	//ecurve2_add(listP[HalfPow2w3NAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2wNAF + SDR[i]], R);
	}

	for (int i = 1; i < Pow2wNAF; i += 2) {
		epoint_free(listP[i]);
	}
}

void ShamirMul_wNAF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[600], SDRb[600];
	DWORD lenSDRa, lenSDRb;
	pepoint listP[Pow2wNAF], listQ[Pow2wNAF];
	lenSDRa = GenwNAF(a, SDRa);
	lenSDRb = GenwNAF(b, SDRb);
	PreMul_wNAF(P, listP);
	PreMul_wNAF(Q, listQ);
	while (lenSDRa < lenSDRb) SDRa[lenSDRa++] = 0;
	while (lenSDRb < lenSDRa) SDRb[lenSDRb++] = 0;

	epoint2_set(0, 0, 1, R);
	for (int i = lenSDRa - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2wNAF + SDRa[i]], R);
		ecurve2_add(listQ[HalfPow2wNAF + SDRb[i]], R);
	}

	for (int i = 1; i < Pow2wNAF; i += 2) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
	}
}

////////////////////////////////////////////////////////
void PreMul_wNAF(pepoint P, pepoint *listP, char *eSDR, DWORD lene)
{
	for (int i = 1; i < HalfPow2wNAF; i += 2) {
		listP[HalfPow2wNAF - i] = epoint_init();
	}
	if (eSDR[lene - 1] != 1) {
		eSDR[lene++] = 1;
	}
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	epoint_copy(P, listP[HalfPow2wNAF - 1]);
	for (int i = lene - 2, j, k; i >= 0; i--) {
		j = eSDR[i] - eSDR[i + 1];
		k = HalfPow2wNAF - eSDR[i];
		ecurve2_add(listP[k + j], listP[k]);
		do {
			ecurve2_add(P2, listP[k]);
			j -= 2;
		} while (j > 0);
	}
	epoint_free(P2);
}

void ScalarMul_wNAF1(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[HalfPow2wNAF];

	lenSDR = GenwNAF(k, SDR, eSDR, lene);
	PreMul_wNAF(P, listP, eSDR, lene);
	epoint2_set(0, 0, 1, R);
	ecurve2_add(listP[HalfPow2wNAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 2; i >= 0; i--) {
		ecurve2_add(R, R);
		if (SDR[i] > 0) {
			ecurve2_add(listP[HalfPow2wNAF - SDR[i]], R);
		}
		if (SDR[i] < 0) {
			ecurve2_sub(listP[HalfPow2wNAF + SDR[i]], R);
		}
	}
	for (int i = 1; i < HalfPow2wNAF; i += 2) {
		epoint_free(listP[HalfPow2wNAF - i]);
	}
}

void ShamirMul_wNAF1(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[MAX_M + 1] = { 0 };
	char SDRb[MAX_M + 1] = { 0 };
	char eSDRa[Pow2wNAF / 2], eSDRb[Pow2wNAF / 2];
	DWORD lenSDRa, lenSDRb, lenea, leneb, lenSDR;
	pepoint listP[HalfPow2wNAF + 1];
	pepoint listQ[HalfPow2wNAF + 1];
	listP[HalfPow2wNAF] = epoint_init();
	listQ[HalfPow2wNAF] = epoint_init();
	lenSDRa = GenwNAF(a, SDRa, eSDRa, lenea);
	lenSDRb = GenwNAF(b, SDRb, eSDRb, leneb);
	PreMul_wNAF(P, listP, eSDRa, lenea);
	PreMul_wNAF(Q, listQ, eSDRb, leneb);
	lenSDR = max(lenSDRa, lenSDRb);
	epoint2_set(0, 0, 1, R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		if (SDRa[i] > 0) ecurve2_add(listP[HalfPow2wNAF - SDRa[i]], R);
		else ecurve2_sub(listP[HalfPow2wNAF + SDRa[i]], R);
		if (SDRb[i] > 0) ecurve2_add(listQ[HalfPow2wNAF - SDRb[i]], R);
		else ecurve2_sub(listQ[HalfPow2wNAF + SDRb[i]], R);
	}
	for (int i = 1; i < HalfPow2wNAF; i += 2) {
		epoint_free(listP[HalfPow2wNAF - i]);
		epoint_free(listQ[HalfPow2wNAF - i]);
	}
}

/***************************************************************************
Ignore below, please!!!
****************************************************************************/

///////////////////////////////////////////////////
DWORD Genw3NAF(big k, char *SDR)
{
	int eSDR[Pow2w3NAF];
	for (int i = 1; i < (Pow2w3NAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2w3NAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + (y << 5)), d);
		SDR[lenSDR++] = eSDR[d->w[0] & MAXe3SDR];
		sftbit(d, -w3NAF, d);
		if (SDR[lenSDR - 1] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < w3NAF; i++) SDR[lenSDR++] = 0;
	}
	mirkill(d);
	return lenSDR;
}

void PreMul_w3NAF(pepoint P, pepoint *listP)
{
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	for (int i = 1; i < Pow2w3NAF; i += 2) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[HalfPow2w3NAF + 1]);
	ecurve2_sub(listP[HalfPow2w3NAF + 1], listP[HalfPow2w3NAF - 1]);
	for (int i = 3; i < HalfPow2w3NAF; i += 2) {
		epoint2_copy(listP[HalfPow2w3NAF + i - 2], listP[HalfPow2w3NAF + i]);
		ecurve2_add(P2, listP[HalfPow2w3NAF + i]);
		ecurve2_sub(listP[HalfPow2w3NAF + i], listP[HalfPow2w3NAF - i]);
	}
	listP[HalfPow2w3NAF] = epoint_init();
	epoint_free(P2);
}

void ScalarMul_w3NAF(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[Pow2w3NAF];
	lenSDR = Genw3NAF(k, SDR);
	PreMul_w3NAF(P, listP);
	epoint2_set(0, 0, 1, R);
	//ecurve2_add(listP[HalfPow2w3NAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w3NAF + SDR[i]], R);
	}
	
	for (int i = 1; i < Pow2w3NAF; i += 2) {
		epoint_free(listP[i]);
	}
}

void ShamirMul_w3NAF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[600], SDRb[600];
	DWORD lenSDRa, lenSDRb;
	pepoint listP[Pow2w3NAF], listQ[Pow2w3NAF];
	lenSDRa = Genw3NAF(a, SDRa);
	lenSDRb = Genw3NAF(b, SDRb);
	PreMul_w3NAF(P, listP);
	PreMul_w3NAF(Q, listQ);
	while (lenSDRa < lenSDRb) SDRa[lenSDRa++] = 0;
	while (lenSDRb < lenSDRa) SDRb[lenSDRb++] = 0;

	epoint2_set(0, 0, 1, R);
	for (int i = lenSDRa - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w3NAF + SDRa[i]], R);
		ecurve2_add(listQ[HalfPow2w3NAF + SDRb[i]], R);
	}

	for (int i = 1; i < Pow2w3NAF; i += 2) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
	}
}

////////////////////////////////////////////////////
DWORD Genw4NAF(big k, char *SDR)
{
	int eSDR[Pow2w4NAF];
	for (int i = 1; i < (Pow2w4NAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2w4NAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + (y << 5)), d);
		SDR[lenSDR++] = eSDR[d->w[0] & MAXe4SDR];
		sftbit(d, -w4NAF, d);
		if (SDR[lenSDR - 1] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < w4NAF; i++) SDR[lenSDR++] = 0;
	}
	mirkill(d);
	return lenSDR;
}

void PreMul_w4NAF(pepoint P, pepoint *listP)
{
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	for (int i = 1; i < Pow2w4NAF; i += 2) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[HalfPow2w4NAF + 1]);
	ecurve2_sub(listP[HalfPow2w4NAF + 1], listP[HalfPow2w4NAF - 1]);
	for (int i = 3; i < HalfPow2w4NAF; i += 2) {
		epoint2_copy(listP[HalfPow2w4NAF + i - 2], listP[HalfPow2w4NAF + i]);
		ecurve2_add(P2, listP[HalfPow2w4NAF + i]);
		ecurve2_sub(listP[HalfPow2w4NAF + i], listP[HalfPow2w4NAF - i]);
	}
	listP[HalfPow2w4NAF] = epoint_init();
	epoint_free(P2);
}

void ScalarMul_w4NAF(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[Pow2w4NAF];
	lenSDR = Genw4NAF(k, SDR);
	PreMul_w4NAF(P, listP);
	epoint2_set(0, 0, 1, R);
	//ecurve2_add(listP[HalfPow2w3NAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w4NAF + SDR[i]], R);
	}
	for (int i = 1; i < Pow2w4NAF; i += 2) {
		epoint_free(listP[i]);
	}
}

void ShamirMul_w4NAF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[600], SDRb[600];
	DWORD lenSDRa, lenSDRb;
	pepoint listP[Pow2w4NAF], listQ[Pow2w4NAF];
	lenSDRa = Genw4NAF(a, SDRa);
	lenSDRb = Genw4NAF(b, SDRb);
	PreMul_w4NAF(P, listP);
	PreMul_w4NAF(Q, listQ);
	while (lenSDRa < lenSDRb) SDRa[lenSDRa++] = 0;
	while (lenSDRb < lenSDRa) SDRb[lenSDRb++] = 0;

	epoint2_set(0, 0, 1, R);
	for (int i = lenSDRa - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w4NAF + SDRa[i]], R);
		ecurve2_add(listQ[HalfPow2w4NAF + SDRb[i]], R);
	}

	for (int i = 1; i < Pow2w4NAF; i += 2) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
	}
}

////////////////////////////////////////////////////////
DWORD Genw5NAF(big k, char *SDR)
{
	int eSDR[Pow2w5NAF];
	for (int i = 1; i < (Pow2w5NAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2w5NAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + (y << 5)), d);
		SDR[lenSDR++] = eSDR[d->w[0] & MAXe5SDR];
		sftbit(d, -w5NAF, d);
		if (SDR[lenSDR - 1] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < w5NAF; i++) SDR[lenSDR++] = 0;
	}
	mirkill(d);
	return lenSDR;
}

void PreMul_w5NAF(pepoint P, pepoint *listP)
{
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	for (int i = 1; i < Pow2w5NAF; i += 2) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[HalfPow2w5NAF + 1]);
	ecurve2_sub(listP[HalfPow2w5NAF + 1], listP[HalfPow2w5NAF - 1]);
	for (int i = 3; i < HalfPow2w5NAF; i += 2) {
		epoint2_copy(listP[HalfPow2w5NAF + i - 2], listP[HalfPow2w5NAF + i]);
		ecurve2_add(P2, listP[HalfPow2w5NAF + i]);
		ecurve2_sub(listP[HalfPow2w5NAF + i], listP[HalfPow2w5NAF - i]);
	}
	listP[HalfPow2w5NAF] = epoint_init();
	epoint_free(P2);
}

void ScalarMul_w5NAF(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[Pow2w5NAF];
	lenSDR = Genw5NAF(k, SDR);
	PreMul_w5NAF(P, listP);
	epoint2_set(0, 0, 1, R);
	//ecurve2_add(listP[HalfPow2w3NAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w5NAF + SDR[i]], R);
	}
	for (int i = 1; i < Pow2w5NAF; i += 2) {
		epoint_free(listP[i]);
	}
}

void ShamirMul_w5NAF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[600], SDRb[600];
	DWORD lenSDRa, lenSDRb;
	pepoint listP[Pow2w5NAF], listQ[Pow2w5NAF];
	lenSDRa = Genw5NAF(a, SDRa);
	lenSDRb = Genw5NAF(b, SDRb);
	PreMul_w5NAF(P, listP);
	PreMul_w5NAF(Q, listQ);
	while (lenSDRa < lenSDRb) SDRa[lenSDRa++] = 0;
	while (lenSDRb < lenSDRa) SDRb[lenSDRb++] = 0;

	epoint2_set(0, 0, 1, R);
	for (int i = lenSDRa - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w5NAF + SDRa[i]], R);
		ecurve2_add(listQ[HalfPow2w5NAF + SDRb[i]], R);
	}

	for (int i = 1; i < Pow2w5NAF; i += 2) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
	}
}

///////////////////////////////////////////////////////////
DWORD Genw6NAF(big k, char *SDR)
{
	int eSDR[Pow2w6NAF];
	for (int i = 1; i < (Pow2w6NAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2w6NAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + (y << 5)), d);
		SDR[lenSDR++] = eSDR[d->w[0] & MAXe6SDR];
		sftbit(d, -w6NAF, d);
		if (SDR[lenSDR - 1] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < w6NAF; i++) SDR[lenSDR++] = 0;
	}
	mirkill(d);
	return lenSDR;
}

void PreMul_w6NAF(pepoint P, pepoint *listP)
{
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	for (int i = 1; i < Pow2w6NAF; i += 2) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[HalfPow2w6NAF + 1]);
	ecurve2_sub(listP[HalfPow2w6NAF + 1], listP[HalfPow2w6NAF - 1]);
	for (int i = 3; i < HalfPow2w6NAF; i += 2) {
		epoint2_copy(listP[HalfPow2w6NAF + i - 2], listP[HalfPow2w6NAF + i]);
		ecurve2_add(P2, listP[HalfPow2w6NAF + i]);
		ecurve2_sub(listP[HalfPow2w6NAF + i], listP[HalfPow2w6NAF - i]);
	}
	listP[HalfPow2w6NAF] = epoint_init();
	epoint_free(P2);
}

void ScalarMul_w6NAF(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[Pow2w6NAF];
	lenSDR = Genw6NAF(k, SDR);
	PreMul_w6NAF(P, listP);
	epoint2_set(0, 0, 1, R);
	//ecurve2_add(listP[HalfPow2w3NAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w6NAF + SDR[i]], R);
	}
	for (int i = 1; i < Pow2w6NAF; i += 2) {
		epoint_free(listP[i]);
	}
}

void ShamirMul_w6NAF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[600], SDRb[600];
	DWORD lenSDRa, lenSDRb;
	pepoint listP[Pow2w6NAF], listQ[Pow2w6NAF];
	lenSDRa = Genw6NAF(a, SDRa);
	lenSDRb = Genw6NAF(b, SDRb);
	PreMul_w6NAF(P, listP);
	PreMul_w6NAF(Q, listQ);
	while (lenSDRa < lenSDRb) SDRa[lenSDRa++] = 0;
	while (lenSDRb < lenSDRa) SDRb[lenSDRb++] = 0;

	epoint2_set(0, 0, 1, R);
	for (int i = lenSDRa - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w6NAF + SDRa[i]], R);
		ecurve2_add(listQ[HalfPow2w6NAF + SDRb[i]], R);
	}

	for (int i = 1; i < Pow2w6NAF; i += 2) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
	}
}

///////////////////////////////////////////////////////////
DWORD Genw7NAF(big k, char *SDR)
{
	int eSDR[Pow2w7NAF];
	for (int i = 1; i < (Pow2w7NAF >> 1); i += 2) {
		eSDR[i] = i; eSDR[Pow2w7NAF - i] = -i;
	}
	big d = mirvar(3); copy(k, d);
	DWORD lenSDR = 0;	int x, y;
	while (d->len) {
		x = 0; y = 0;
		while (!(d->w[y] & (1 << x))) {
			x++; if (!(x % 32)) y++;
			SDR[lenSDR++] = 0;
		}
		sftbit(d, -(x + (y << 5)), d);
		SDR[lenSDR++] = eSDR[d->w[0] & MAXe7SDR];
		sftbit(d, -w7NAF, d);
		if (SDR[lenSDR - 1] < 0) incr(d, 1, d);
		if (d->len)for (int i = 1; i < w7NAF; i++) SDR[lenSDR++] = 0;
	}
	mirkill(d);
	return lenSDR;
}

void PreMul_w7NAF(pepoint P, pepoint *listP)
{
	pepoint P2 = epoint_init();
	epoint2_copy(P, P2);
	ecurve2_add(P2, P2);
	for (int i = 1; i < Pow2w7NAF; i += 2) {
		listP[i] = epoint_init();
	}
	epoint2_copy(P, listP[HalfPow2w7NAF + 1]);
	ecurve2_sub(listP[HalfPow2w7NAF + 1], listP[HalfPow2w7NAF - 1]);
	for (int i = 3; i < HalfPow2w7NAF; i += 2) {
		epoint2_copy(listP[HalfPow2w7NAF + i - 2], listP[HalfPow2w7NAF + i]);
		ecurve2_add(P2, listP[HalfPow2w7NAF + i]);
		ecurve2_sub(listP[HalfPow2w7NAF + i], listP[HalfPow2w7NAF - i]);
	}
	listP[HalfPow2w7NAF] = epoint_init();
	epoint_free(P2);
}

void ScalarMul_w7NAF(big k, pepoint P, pepoint R)
{
	char SDR[600], eSDR[200];
	DWORD lenSDR, lene;
	pepoint listP[Pow2w7NAF];
	lenSDR = Genw7NAF(k, SDR);
	PreMul_w7NAF(P, listP);
	epoint2_set(0, 0, 1, R);
	//ecurve2_add(listP[HalfPow2w3NAF - SDR[lenSDR - 1]], R);
	for (int i = lenSDR - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w7NAF + SDR[i]], R);
	}
	for (int i = 1; i < Pow2w7NAF; i += 2) {
		epoint_free(listP[i]);
	}
}

void ShamirMul_w7NAF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	char SDRa[600], SDRb[600];
	DWORD lenSDRa, lenSDRb;
	pepoint listP[Pow2w7NAF], listQ[Pow2w7NAF];
	lenSDRa = Genw7NAF(a, SDRa);
	lenSDRb = Genw7NAF(b, SDRb);
	PreMul_w7NAF(P, listP);
	PreMul_w7NAF(Q, listQ);
	while (lenSDRa < lenSDRb) SDRa[lenSDRa++] = 0;
	while (lenSDRb < lenSDRa) SDRb[lenSDRb++] = 0;

	epoint2_set(0, 0, 1, R);
	for (int i = lenSDRa - 1; i >= 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listP[HalfPow2w7NAF + SDRa[i]], R);
		ecurve2_add(listQ[HalfPow2w7NAF + SDRb[i]], R);
	}

	for (int i = 1; i < Pow2w7NAF; i += 2) {
		epoint_free(listP[i]);
		epoint_free(listQ[i]);
	}
}