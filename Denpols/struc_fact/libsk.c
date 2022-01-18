#include "libsk.h"


void calc_sk(void *qinputv, void *datav, void *outputv, int Ldata, int qbins)
{
        int i,k;
        double qx, qy, qz;
        double q;

        double * qinput = (double *) qinputv;
        double * data = (double *) datav;
        double * output = (double *) outputv;

        double tempsin, tempcos;

        for(k = 0; k < qbins; k++){

                qx = qinput[3*k]; qy = qinput[3*k+1]; qz = qinput[3*k+2];
                tempsin = 0; tempcos = 0;

                for(i=0; i < Ldata; i++){
                        tempsin += sin(qx*data[3*i] + qy*data[3*i+1] + qz*data[3*i+2]);
                        tempcos += cos(qx*data[3*i] + qy*data[3*i+1] + qz*data[3*i+2]);
                }

                output[k] += tempsin*tempsin + tempcos*tempcos;
        }

}


void calc_sk_2(void *qinputv, void *datav, void *outputv, int Ldata, int qbins)
{
        int i,k,j;
        double qx, qy, qz;
        double q;

        double * qinput = (double *) qinputv;
        double * data = (double *) datav;
        double * output = (double *) outputv;

	double dx, dy, dz;
        double temp;

	output[0] += Ldata*(Ldata-1);
        for(k = 1; k < qbins; k++){

                qx = qinput[3*k]; qy = qinput[3*k+1]; qz = qinput[3*k+2];
                temp = 0;

                for(i=0; i < Ldata-1; i++){
			for(j = i+1; j < Ldata; j++){
                        	dx = data[3*i] - data[3*j];
				dy = data[3*i+1] - data[3*j+1];
				dz = data[3*i+2] - data[3*j+2];
				temp += sin(qx*dx + qy*dy + qz*dz)/(qx*dx + qy*dy + qz*dz);
                	
				if(temp != temp){ printf("error! %d %d %.8e %.8e %.8e %.8e\n", i,j, dx, dy, dz, temp); exit(1);} 
			}
		}

                output[k] += temp;
        }

}



