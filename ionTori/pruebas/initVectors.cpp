#include "nr.h"
#include "nrutil.h"
#include "functions.h"
#include <stdio.h>
#include <math.h>



void initVectors() {
	
	extern double *r, *rCells, *ne, paso,rMax,rMin;
	extern int nR;
	FILE *fne;
	
	fne=fopen("ne.txt","r");
	if (NULL == fne) {
        perror("File not found!");
    }
	
	r=dvector(1,nR);
	rCells=dvector(0,nR);
	ne=dvector(1,nR);
	
	paso=pow(rMax/rMin,1.0/nR);
	rCells[0]=rMin;
	for(int i=1;i<=nR;i++) {
		double dens;
		rCells[i]=rCells[i-1]*paso;
		r[i]=sqrt(rCells[i]*rCells[i-1]);
		if(fscanf(fne,"%lf",&dens) > 0) {
			ne[i]=dens;
		}
	}
	fclose(fne);
}

void interpol(double rad, double *result, double r[], double ne[]) {
	
	int i=0;
	int j=1;
	extern int nR;
	if (rad < r[1] || rad >= r[nR]) {
		*result = 0.0;
	} else {
		while (i==0) {
			if (rad >= r[j] && rad < r[j+1]) i=j;
			j++;
		}
	}
	double m=log10(ne[i+1]/ne[i])/log10(r[i+1]/r[i]);
	*result = ne[i]*pow(rad/r[i],m);
}

double eDensity(double rad, double theta) {
	extern double *r, *ne;
	extern int nR;
	double result,err;
	if (rad >= r[nR] || rad <= r[1] || theta >= PI/3.0) {
		result=0.0;
	} else {
		ratint(r,ne,nR,rad,&result,&err);
		//interpol(rad,&result,r,ne);
	}
	return result;
}