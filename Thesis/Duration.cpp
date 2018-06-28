#include "Duration.h"
#include "ShamirMul.h"

void startTimer(stopWatch *timer) {
	QueryPerformanceCounter(&timer->start);
}

void stopTimer(stopWatch *timer) {
	QueryPerformanceCounter(&timer->stop);
}

double LIToSecs(LARGE_INTEGER *L) {
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	return ((double)L->QuadPart / (double)frequency.QuadPart);
}

double getElapsedTime(stopWatch *timer) {
	LARGE_INTEGER time;
	time.QuadPart = timer->stop.QuadPart - timer->start.QuadPart;
	return LIToSecs(&time);
}

DWORD getQuadPart(stopWatch *timer) {
	return timer->stop.QuadPart - timer->start.QuadPart;
}

void GenDuration(big k, char* r, int code,
	void(*func) (big, char*)) {
	double dur;
	stopWatch timer;
	startTimer(&timer);
	(*func)(k, r);
	stopTimer(&timer);
	dur = getElapsedTime(&timer) * 1000;
};

void SclDuration(big k, pepoint P, pepoint R, 
	void(*func) (big, pepoint, pepoint), double &min)
{
	double dur;
	stopWatch timer;
	startTimer(&timer);
	(*func)(k, P, R);
	stopTimer(&timer);
	dur = getElapsedTime(&timer) * 1000;
	min = (min < dur) ? min : dur;
}

void ShrDuration(big k, pepoint P, pepoint R, 
	void(*func) (big, pepoint, big, pepoint, pepoint), double &min)
{
	double dur;
	stopWatch timer;
	startTimer(&timer);
	big a = mirvar(1);
	big b = mirvar(1);
	pepoint Q = epoint_init();
	/*ShamirMul(k, P, a, Q, b);
	(*func) (a, P, b, Q, R);
	epoint_free(Q);
	mirkill(a); mirkill(b);*/
	ShamirMul(k, P, R, (*func));
	stopTimer(&timer);
	dur = getElapsedTime(&timer) * 1000;
	min = (min < dur) ? min : dur;
}

void ShrDuration(big a, pepoint P, big b, pepoint Q, pepoint R, 
	void(*func) (big, pepoint, big, pepoint, pepoint), double &min)
{
	double dur;
	stopWatch timer;
	startTimer(&timer);
	(*func) (a, P, b, Q, R);
	stopTimer(&timer);
	dur = getElapsedTime(&timer) * 1000;
	min = (min < dur) ? min : dur;
}

void GenFormSclDuration(big k, DWORD (*func) (big, char*), double &min, DWORD &ticks) {
	
	char form[1000] = { 0 };
	double dur;
	stopWatch timer;
	startTimer(&timer);
	(*func) (k, form);
	stopTimer(&timer);
	dur = getElapsedTime(&timer) * 1000;
	min = (min < dur) ? min : dur;
	DWORD qP = getQuadPart(&timer);
	ticks = (ticks < qP) ? ticks : qP;
}

void GenFormShrDuration(big a, big b, DWORD(*func) (big, char*), double &min, DWORD &ticks) {
	char form[1000] = { 0 };
	double dur;
	stopWatch timer;
	startTimer(&timer);
	(*func) (a, form);
	(*func) (b, form);
	stopTimer(&timer);
	dur = getElapsedTime(&timer) * 1000;
	min = (min < dur) ? min : dur;
	DWORD qP = getQuadPart(&timer);
	ticks = (ticks < qP) ? ticks : qP;
}
