#include "option_AKDAC.h"

void R1(big d, big e, pepoint Ru, pepoint Rv) {
	subtract(d, e, d);
	sftbit(d, -1, d);
	ecurve2_add(Ru, Rv);
	ecurve2_add(Ru, Ru);
}

void R2(big d, pepoint Ru) {
	sftbit(d, -1, d);
	ecurve2_add(Ru, Ru);
}

void ShamirMul_AKDAC(big a, pepoint P, big b, pepoint Q, pepoint R) {

	big d = mirvar(1), e = mirvar(1);
	pepoint Ru, Rv;
	Ru = epoint_init();
	Rv = epoint_init();
	copy(a, d); copy(b, e);
	epoint2_copy(P, Ru);
	epoint2_copy(Q, Rv);
	int dmod2, emod2;
	while (compare(d, e)) {
		dmod2 = d->w[0] & 1;
		emod2 = e->w[0] & 1;
		if (dmod2 && emod2) {
			if (compare(d, e) > 0) R1(d, e, Ru, Rv);
			else R1(e, d, Rv, Ru);
		}
		else {
			if (!dmod2) R2(d, Ru);
			else R2(e, Rv);
		}
	}
	ecurve2_add(Rv, Ru);
	ecurve2_mult(d, Ru, R);
}

void TestShrMulAKDAC(csprng &Rng, pepoint P) {
	big a = mirvar(0x29C), b = mirvar(0xD1);
	pepoint Q = epoint_init();
	pepoint R = epoint_init();
	pepoint R1 = epoint_init();
	ecurve2_mult(a, P, Q);
	int count = 0;
	for (int i = 0; i < 10000; i++) {
		strong_bigdig(&Rng, 50, 10, a);
		strong_bigdig(&Rng, 50, 10, b);
		ShamirMul_AKDAC(a, P, b, Q, R); 
		ecurve2_mult2(a, P, b, Q, R1); 

		count += epoint2_comp(R, R1);
		if (!epoint2_comp(R, R1)) {
			std::cout << "a: "; cotnum(a, stdout);
			std::cout << "b: "; cotnum(b, stdout);
			std::cout << "R: "; cotnumEp(R);
			std::cout << "R1: "; cotnumEp(R1);
			break;
		}
	}
	std::cout << "Cmp: " << count << std::endl;
}

