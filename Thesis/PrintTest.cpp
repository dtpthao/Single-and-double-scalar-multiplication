#include "PrintTest.h"

void printCmpSclShr_allm(const char* algName, int* m,
	double tmin[6][LIB + 1], double per[6][LIB + 1], int index) {
	int index2 = index + 1;
	if (index == 0) index2++;
	Line2;
	printf("_____________________%4s______________________\t", algName);
	printf("_____________________%4s______________________\n", algName);
	for (int i = 0; i < 5; i++) {
		printf("%4d: %13.5f %13.5f %13.5f\t",
			m[i], tmin[i][index], tmin[i][index2], tmin[i][LIB]);
		printf("%4d: %12.1f%% %12.1f%% %12.1f%%\n",
			m[i], per[i][index], per[i][index2], per[i][LIB]);
	}
}

void printAll_tmin(int* m, double tmin[6][LIB + 1]) {

	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for (int j = 0; j <= LIB; j += 2) PRINTTMIN(tmin[i][j]);
		cout << endl;
	}
	cout << endl;
}

void printAll_per(int* m, double per[6][LIB + 1]) {
	for (int i = 0; i < 6; i++) {
		printf("%%%4d:\t", m[i]);
		for (int j = 0; j <= LIB; j += 2) PRINTPERCENT(per[i][j]);
		cout << endl;
	}
	cout << endl;
}

void printAll_corr(int* m, unsigned int correct[6][LIB]) {
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for (int j = 0; j < LIB; j += 2) PRINTCORRECT(correct[i][j]);
		printf("    /%5d\n", TESTALL);
	}
}

void printAll_m(int* m, double tmin[6][LIB + 1],
	double per[6][LIB + 1], unsigned int correct[6][LIB])
{
	for (int i = 0; i <= LIB; i++) {
		per[5][i] = 0;
		for (int j = 0; j < 5; j++) {
			per[5][i] += per[j][i];
		}
		per[5][i] /= 5;
	}
	cout << "  m\t";
	PRINTALL;
	printAll_tmin(m, tmin);
	printAll_per(m, per);
	printAll_corr(m, correct);

	Line;
	cout << endl;
	cout << "  m " << setw(13) << "kP" << setw(16) << "aP + bQ";
	cout << setw(12) << "Lib" << endl;
	printCmpSclShr_allm("Bin", m, tmin, per, SclBin);
	printCmpSclShr_allm("NAF", m, tmin, per, SclNAF);
	printCmpSclShr_allm("wNAF", m, tmin, per, SclwNAF);
	printCmpSclShr_allm("JT", m, tmin, per, SclJT);
}

void printAll_Shr(int* m, double tmin[6][LIB + 1],
	double per[6][LIB + 1], unsigned int correct[6][LIB])
{
	cout << "Compare Shamir options" << endl;
	Line;
	cout << "  m\t";
	PRINTSHR;
	for (int i = ShrBin; i <= LIB; i++) {
		per[5][i] = 0;
		for (int j = 0; j < 5; j++) {
			per[5][i] += per[j][i];
		}
		per[5][i] /= 5;
	}
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for (int j = ShrBin; j <= LIB; j += 2) PRINTTMIN(tmin[i][j]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 6; i++) {
		printf("%%%4d:\t", m[i]);
		for (int j = ShrBin; j <= LIB; j += 2) PRINTPERCENT(per[i][j]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for (int j = ShrBin; j < LIB; j += 2) PRINTCORRECT(correct[i][j]);
		printf("    /%5d\n", TESTSHR);
	}
}

void printAll_Scl(int* m, double tmin[6][LIB + 1],
	double per[6][LIB + 1], unsigned int correct[6][LIB])
{
	cout << "Compare Scalar methods" << endl;
	Line;
	cout << "  m\t";
	PRINTSCL;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		PRINTTMIN(tmin[i][SclBin]);
		PRINTTMIN(tmin[i][SclNAF]);
		PRINTTMIN(tmin[i][SclwNAF]);
		PRINTTMIN(tmin[i][SclJT]);
		PRINTTMIN(tmin[i][LIB]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 6; i++) {
		printf("%5d:\t", m[i]);
		PRINTPERCENT(per[i][SclBin]);
		PRINTPERCENT(per[i][SclNAF]);
		PRINTPERCENT(per[i][SclwNAF]);
		PRINTPERCENT(per[i][SclJT]);
		PRINTPERCENT(per[i][LIB]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		PRINTCORRECT(correct[i][SclBin]);
		PRINTCORRECT(correct[i][SclNAF]);
		PRINTCORRECT(correct[i][SclwNAF]);
		PRINTCORRECT(correct[i][SclJT]);
		printf("    /%5d\n", TESTALL);
	}
}


void printAll_SclwNAF(int* m, double tmin[6][7],
	double per[6][7], unsigned int correct[6][7])
{
	for (int i = 0; i <= 5; i++) {
		for (int j = 0; j < 5; j++) {
			per[5][i] += per[j][i];
		}
		per[5][i] /= 5;
	}
	cout << "Compare different values of w for wNAF Scalar" << endl;
	Line;
	cout << "  m\t";
	PRINTSCLwNAF;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for(int j = 0; j< 7; j++) PRINTTMIN(tmin[i][j]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for (int j = 0; j< 7; j++) PRINTPERCENT(per[i][j]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		for (int j = 0; j< 6; j++) PRINTCORRECT(correct[i][j]);
		printf("    /%5d\n", TESTwNAF);
	}
}

void printGenFormScl(int* m, double tmin[6][LIB + 1], DWORD task[6][LIB + 1]) {

	cout << "Time of computing number representations (k)" << endl;
	Line;
	cout << "  m\t";
	PRINTGENFORMSCL;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		PRINTTMIN(tmin[i][SclNAF]);
		PRINTTMIN(tmin[i][SclwNAF]);
		PRINTTMIN(tmin[i][SclJT]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		PRINTTASK(task[i][SclNAF]);
		PRINTTASK(task[i][SclwNAF]);
		PRINTTASK(task[i][SclJT]);
		cout << endl;
	}
	cout << endl;
}

void printGenFormShr(int* m, double tmin[6][LIB + 1], DWORD task[6][LIB + 1]) {

	cout << "Time of computing number representations (a, b)" << endl;
	Line;
	cout << "  m\t";
	PRINTGENFORMSHR;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		PRINTTMIN(tmin[i][ShrNAF]);
		PRINTTMIN(tmin[i][ShrwNAF]);
		PRINTTMIN(tmin[i][ShrJSF]);
		PRINTTMIN(tmin[i][ShrJT]);
		PRINTTMIN(tmin[i][ShrBNBC]);
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 5; i++) {
		printf("%5d:\t", m[i]);
		PRINTTASK(task[i][ShrNAF]);
		PRINTTASK(task[i][ShrwNAF]);
		PRINTTASK(task[i][ShrJSF]);
		PRINTTASK(task[i][ShrJT]);
		PRINTTASK(task[i][ShrBNBC]);
		cout << endl;
	}
	cout << endl;
}
