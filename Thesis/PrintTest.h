#ifndef _PRINTTEST_H
#define _PRINTTEST_H

#include "Test.h"

#define Line cout << "________________________________________________\n\n"
#define Line1 cout << "_______________________________________________\n"
#define Line2 {\
cout << "_______________________________________________\t";\
cout << "_______________________________________________\n";\
}

#define PRINTALL {\
printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s\n", \
		"SclBin", "ShrBin", "ShrJSF", "ShrNAF", "ShrwNAF", \
		"ShrJT", "ShrBNBC", "ShrAKDAC", "Lib");\
}

#define PRINTSHR {\
printf("%10s %10s %10s %10s %10s %10s %10s %10s\n", \
		"ShrBin", "ShrJSF", "ShrNAF", "ShrwNAF", \
		"ShrJT", "ShrBNBC", "ShrAKDAC", "Lib");\
}

#define PRINTSCL {\
printf("%10s %10s %10s %10s %10s\n", \
		"SclBin", "SclNAF", "SclwNAF", "SclJT", "Lib");\
}

#define PRINTSCLwNAF {\
printf("%10s %10s %10s %10s %10s %10s %10s\n", \
		"SclBin", "Sclw3NAF", "Sclw4NAF", "Sclw5NAF", "Sclw6NAF", "Sclw7NAF", "Lib");\
}

#define PRINTGENFORMSCL {\
printf("%10s %10s %10s\n", \
		"GenNAF", "GenwNAF", "GenJT");\
}

#define PRINTGENFORMSHR {\
printf("%10s %10s %10s %10s %10s\n", \
		"GenNAF", "GenwNAF", "GenJSF", "GenJT", "GenBNBC");\
}

#define PRINTTMIN(tmin) printf("%10.5f ", tmin);
#define PRINTCORRECT(correct) printf("%10d ", correct);
#define PRINTPERCENT(per) printf("%9.1f%% ", per);
#define PRINTTASK(task) printf("%10d ", task);

void printAll_tmin(int* m, double tmin[6][LIB + 1]);
void printAll_per(int* m, double per[6][LIB + 1]);
void printAll_corr(int* m, unsigned int correct[6][LIB]);

/*

*/
void printAll_m(int* m, double tmin[6][LIB + 1],
	double per[6][LIB + 1], unsigned int correct[6][LIB]);

void printAll_Shr(int* m, double tmin[6][LIB + 1],
	double per[6][LIB + 1], unsigned int correct[6][LIB]);

void printAll_Scl(int* m, double tmin[6][LIB + 1],
	double per[6][LIB + 1], unsigned int correct[6][LIB]);

void printGenFormScl(int* m, double tmin[6][LIB + 1], DWORD task[6][LIB + 1]);
void printGenFormShr(int* m, double tmin[6][LIB + 1], DWORD task[6][LIB + 1]);

void printAll_SclwNAF(int* m, double tmin[6][7],
	double per[6][7], unsigned int correct[6][7]);

#endif