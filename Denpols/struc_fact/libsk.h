#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//void calc_fouriergr(void *datav, void *outputv, int qbins, double dq, int nbins, double dr);
void calc_sk(void *qinputv, void *datav, void *outputv, int Ldata, int qbins);
void calc_sk_2(void *qinputv, void *datav, void *outputv, int Ldata, int qbins);
