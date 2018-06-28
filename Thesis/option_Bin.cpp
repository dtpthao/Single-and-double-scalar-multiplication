#include "option_Bin.h"
#include <Windows.h>

//ScalarMul
void ScalarMul_Bin(big k, pepoint P, pepoint R)
{
	unsigned int len, i = 31;
	pepoint D = epoint_init();
	epoint2_copy(P, D);
	epoint_set(0, 0, 1, R);
	while (!(k->w[k->len - 1] & (1 << i))) i--;
	len = (k->len << 5) - (31 - i);
	for (i = 0; i < len; i++) {
		if (k->w[i / 32] & (1 << i)) {
			ecurve2_add(D, R);
		}
		ecurve2_add(D, D);
	}
	epoint_free(D);
}

//ShamirMul
void PreMul_Bin(pepoint P, pepoint Q, pepoint *plist)
{
	for (int i = 0; i < 4; i++) {
		plist[i] = epoint_init();
	}
	epoint2_set(0, 0, 1, plist[0]);
	epoint2_copy(P, plist[1]);
	epoint2_copy(Q, plist[2]);
	epoint2_copy(P, plist[3]);
	ecurve2_add(Q, plist[3]);
}

void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	pepoint plist[4];
	int bita, bitb, i, j = 0;
	DWORD shift = 1, lastw;
	lastw = (compare(a, b) >= 0) ? a->w[a->len - 1] : b->w[b->len - 1];
	//printf("lastw: %x\n", lastw);
	while (lastw >> j && j != 31) {
		j++;
	}
	PreMul_Bin(P, Q, plist);
	epoint_set(0, 0, 1, R);
	i = (lastw == a->w[a->len - 1]) ? a->len - 1 : b->len - 1;
	for (shift <<= j; i >= 0; i--, shift = 0x80000000, j = 31) {
		while (shift) {
			bita = (a->w[i] & shift) >> j;
			bitb = (b->w[i] & shift) >> j;
			ecurve2_add(R, R);
			ecurve2_add(plist[(bitb << 1) + bita], R);
			shift >>= 1;
			j--;
		}
	}
	for (int i = 0; i < 4; i++) epoint_free(plist[i]);
}

void TestShrMul_Bin(csprng &Rng, pepoint P, big n) {

}



//WORD GetBin2(big k, char*Bit)
//{
//	DWORD i = 0;
//	big d = mirvar(1);
//	copy(k, d);
//	while (d->len != 1 || d->w[0] != 1) {
//		if (d->w[0] & 1) Bit[i] = 1;
//		else Bit[i] = 0;
//		sftbit(d, -1, d);
//		i++;
//	}
//	Bit[i] = 1;
//	mirkill(d);
//	return ++i;
//}


//void ShamirMul_Bin(big a, pepoint P, big b, pepoint Q, pepoint R, DWORD len)
//{
//	pepoint plist[4];
//	//char Bita[(MAX_M + 1) / 2] = { 0 };
//	//char Bitb[(MAX_M + 1) / 2] = { 0 };
//	PreMul_Bin(P, Q, plist);
//	epoint_set(0, 0, 1, R);
//	//DWORD lenBita, lenBitb;
//	//lenBita = GetBin(a, Bita);
//	//lenBitb = GetBin(b, Bitb);
//	int bita, bitb, lena = 0, lenb = 0;
//	int i, j = 0; len = 0;// (a->len - 1) << 5;
//	if (a->len == b->len) {
//		i = max(a->w[a->len - 1], b->w[b->len - 1]);
//	}
//	else if (a->len < b->len) i = b->w[b->len - 1];
//	else i = a->w[a->len - 1];
//	//i = a->w[a->len - 1];
//
//	while (i) {
//		len++; i >>= 1;
//	}
//	/*i = b->w[b->len - 1];
//	while (i) {
//	lenb++; i >>= 1;
//	}
//	len = max(lena, lenb);*/
//	//while (j < len) {
//	//	i = j / 32;
//	//	bita = (a->w[i] & (1 << j)) ? 1 : 0;
//	//	bitb = (b->w[i] & (1 << j)) ? 1 : 0;
//	//	ecurve2_add(R, R);
//	//	ecurve2_add(plist[(bitb << 1) + bita], R);
//	//	j++;
//	//	//if (j / 32) { i++; j = 0; }
//	//}
//	DWORD shift = 1 << (len - 1);
//	for (i = a->len - 1; i >= 0; i--) {
//		while (shift) {
//			bita = (a->w[i] & shift) ? 1 : 0;
//			bitb = (b->w[i] & shift) ? 1 : 0;
//			ecurve2_add(R, R);
//			ecurve2_add(plist[(bitb << 1) + bita], R);
//			shift >>= 1;
//		}
//		shift = 0x80000000;
//	}
//	/*for (int i = lenBita - 1; i >= 0; i--) {
//	ecurve2_add(R, R);
//	ecurve2_add(plist[(Bitb[i] << 1) + Bita[i]], R);
//	}*/
//	for (int i = 0; i < 4; i++) epoint_free(plist[i]);
//}



//DWORD GetBin(big k, char*Bit)
//{
//	DWORD lenbit = 0;
//	DWORD shift = 1;
//
//	for (int i = 0; i < k->len - 1; i++, shift = 1) {
//		while (shift) {
//			Bit[lenbit++] = (k->w[i] & shift) ? 1 : 0;
//			shift <<= 1;
//		}
//	}
//	shift = k->w[k->len - 1];
//	while (shift) {
//		Bit[lenbit++] = (shift & 1) ? 1 : 0;
//		shift >>= 1;
//	}
//	return lenbit;
//}
