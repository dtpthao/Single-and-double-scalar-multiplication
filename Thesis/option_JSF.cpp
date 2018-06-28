#include "option_JSF.h"
#include "option_NAF.h"

void SubGenJSF(big Var1, big Var2, char *JSF, DWORD &lenJSF, bool &d, big RS) {
	if (!(Var1->w[0] & 1)) JSF[lenJSF] = 0;
	else {
		JSF[lenJSF] = (Var1->w[0] & 2) ? -1 : 1;
		if ((!(Var1->w[0] & 7 ^ 3) || !(Var1->w[0] & 7 ^ 5)) && ((Var2->w[0] & 3) == 2))
			JSF[lenJSF] = -JSF[lenJSF];
	}
	if (((int)d << 1) == (JSF[lenJSF] + 1)) d = 1 - d;
	sftbit(RS, -1, RS);
}

DWORD GenJSF(big R, big S, char *JSFr, char *JSFs)
{
	big L1 = mirvar(1), L2 = mirvar(1), R1 = mirvar(0), S1 = mirvar(0);
	copy(R, L1); copy(S, L2); copy(R, R1); copy(S, S1);
	bool d1 = 0, d2 = 0;
	DWORD lenJSF = 0;
	while (L1->len > 0 || L2->len > 0) {
		lenJSF++;
		SubGenJSF(L1, L2, JSFr, lenJSF, d1, R1);
		SubGenJSF(L2, L1, JSFs, lenJSF, d2, S1);
		incr(R1, d1, L1); incr(S1, d2, L2);
	}
	return lenJSF;
}

void ShamirMul_JSF(big a, pepoint P, big b, pepoint Q, pepoint R)
{
	pepoint listpoint[9];
	char JSFa[MAX_M + 1] = { 0 };
	char JSFb[MAX_M + 1] = { 0 };
	DWORD lenJSF;
	PreMul(P, Q, listpoint);
	lenJSF = GenJSF(a, b, JSFa, JSFb);
	epoint2_set(0, 0, 1, R);
	for (int i = lenJSF; i > 0; i--) {
		ecurve2_add(R, R);
		ecurve2_add(listpoint[4 - 3 * JSFa[i] - JSFb[i]], R);
	}
	for (int i = 0; i < 9; i++) epoint_free(listpoint[i]);
}