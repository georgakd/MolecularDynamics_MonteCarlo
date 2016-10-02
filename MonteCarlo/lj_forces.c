/*
 *  lj_forces.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 07/03/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void lj_forces(double *x,double *y,double *z);

void lj_forces(double *x,double *y,double *z)
{ 
	double coef1,coef2,sigma6,sigma12,sigma14,sigma8,mysqrt6,mysqrt8,mysqrt12,mysqrt14;
	/////////////forces/////////////////////////////////		  

	
	sumV=0;
	for (i=0;i<npart-1;++i)
	{
		for(j=i+1;j<npart;++j)
		{

			sigma=0.5*(mysigma[i]+mysigma[j]);
			epsilon=sqrt(myepsilon[i]*myepsilon[j]);
			rc=2.5*sigma;   // rc >= 2.5sigma  cutoff distance
			
			xforce=x[j];yforce=y[j];zforce=z[j];
			if (x[j]-x[i]>0.5*Lx) {xforce=x[j]-Lx;}
			if (x[j]-x[i]<-0.5*Lx) {xforce=x[j]+Lx;}
			if (y[j]-y[i]>0.5*Ly) {yforce=y[j]-Ly;}
			if (y[j]-y[i]<-0.5*Ly) {yforce=y[j]+Ly;}
		
			if (z[j]-z[i]>0.5*Lz) {zforce=z[j]-Lz;}
			if (z[j]-z[i]<-0.5*Lz) {zforce=z[j]+Lz;}
	
			r2=(xforce-x[i])*(xforce-x[i])+(yforce-y[i])*(yforce-y[i])+(zforce-z[i])*(zforce-z[i]);
		
			
			r=sqrt(r2);// r calculation
			
			if (r<=rc){
			
			coef1=48*epsilon/(sigma*sigma); // constants for the calculation of LJ forces
			sigma14=pow(sigma,14);																				
			sigma8=pow(sigma,8);
			sigma6=pow(sigma,6);
			sigma12=pow(sigma,12);
			mysqrt6=pow(r,6);// r^6 calculation
			mysqrt8=pow(r,8);// r^8 calculation
			mysqrt12=mysqrt6*mysqrt6;// r^12 calculation
			mysqrt14=mysqrt6*mysqrt8;// r^14 calculation
			coef2=coef1*(sigma14/mysqrt14-0.5*sigma8/mysqrt8);
				
				
			sumV=sumV+sigma12/mysqrt12-sigma6/mysqrt6;
			}
			


		} // end of j
	} // end of i
}