#include "ShamirMul.h"

void ShamirMul(big k, pepoint P, pepoint R,
	void(*func) (big, pepoint, big, pepoint, pepoint))
{
	if (k->len == 0) return;
	pepoint Q = epoint_init();
	big a = mirvar(1);
	big b = mirvar(1);
	DWORD len, i = 31;
	big tmp2l = mirvar(1);
	while (!(k->w[k->len - 1] & (1 << i))) i--;
	len = (k->len << 4) - ((31 - i) >> 1);
	sftbit(tmp2l, len, tmp2l);
	copy(k, a);
	divide(a, tmp2l, b);
	
	epoint2_copy(P, Q);
	for (i = 0; i < len; i++) {
		ecurve2_add(Q, Q);
	}
	//ecurve2_mult(tmp2l, P, Q);
	(*func) (a, P, b, Q, R);

	epoint_free(Q);
	mirkill(a); mirkill(b); mirkill(tmp2l);
}

void ShamirMul(big k, pepoint P, big a, pepoint Q, big b) {
	if (k->len == 0) return;
	DWORD len, i = 31;
	big tmp2l = mirvar(1);
	while (!(k->w[k->len - 1] & (1 << i))) i--;
	len = (k->len << 4) - ((31 - i) >> 1);
	sftbit(tmp2l, len, tmp2l);
	copy(k, a);
	divide(a, tmp2l, b);
	epoint2_copy(P, Q);
	for (i = 0; i < len; i++) {
		ecurve2_add(Q, Q);
	}
	mirkill(tmp2l);
}

void ShamirMul(big k, big a, big b) {
	if (k->len == 0) return;
	DWORD len, i = 31;
	big tmp2l = mirvar(1);
	while (!(k->w[k->len - 1] & (1 << i))) i--;
	len = (k->len << 4) - ((31 - i) >> 1);
	sftbit(tmp2l, len, tmp2l);
	copy(k, a);
	divide(a, tmp2l, b);
	mirkill(tmp2l);
}
