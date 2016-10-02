/*
 *  lj_ts_forces.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 20/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include <math.h>


void lj_ts_forces(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);

void lj_ts_forces(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB)
{
	double coef1,coef2,sigma6,sigma8,sigma12,sigma14,mysqrt6,mysqrt8,mysqrt12,mysqrt14;
	

	
	for(k=0;k<npartA;++k) {sumfxA[k]=0;sumfyA[k]=0;sumfzA[k]=0;}
	for(k=0;k<npartB;++k)   {sumfxB[k]=0;sumfyB[k]=0;sumfzB[k]=0;}

	sumV_ts=0;

	for (i=0;i<npartA;++i)
	{
		for(j=0;j<npartB;++j)
		{
			
			sigma_ts=0.5*(mysigmaA[i]+mysigmaB[j]);
			epsilon_ts=sqrt(myepsilonA[i]*myepsilonB[j]);
			rc_ts=2.5*sigma_ts;   // rc >= 2.5sigma  cutoff distance
			
			xforce=xB[j];yforce=yB[j];zforce=zB[j];

			if (xB[j]-xA[i]>0.5*Lx) {xforce=xB[j]-Lx;}
			if (xB[j]-xA[i]<-0.5*Lx) {xforce=xB[j]+Lx;}
			if (yB[j]-yA[i]>0.5*Ly) {yforce=yB[j]-Ly;}
			if (yB[j]-yA[i]<-0.5*Ly) {yforce=yB[j]+Ly;}
		
			if (zB[j]-zA[i]>0.5*Lz) {zforce=zB[j]-Lz;}
			if (zB[j]-zA[i]<-0.5*Lz) {zforce=zB[j]+Lz;}
	
			r2=(xforce-xA[i])*(xforce-xA[i])+(yforce-yA[i])*(yforce-yA[i])+(zforce-zA[i])*(zforce-zA[i]);		
			r=sqrt(r2);// r calculation	

			if (r<rc_ts){
			
			coef1=48*epsilon_ts/(sigma_ts*sigma_ts); // constants for the calculation of LJ forces
			sigma14=pow(sigma_ts,14);																				
			sigma8=pow(sigma_ts,8);
			sigma6=pow(sigma_ts,6);
			sigma12=pow(sigma_ts,12);
			mysqrt6=pow(r,6);// r^6 calculation
			mysqrt8=pow(r,8);// r^8 calculation
			mysqrt12=mysqrt6*mysqrt6;// r^12 calculation
			mysqrt14=mysqrt6*mysqrt8;// r^14 calculation
			coef2=coef1*(sigma14/mysqrt14-0.5*sigma8/mysqrt8);

				sumfxA[i]=sumfxA[i]+coef2*(xA[i]-xforce);
				sumfyA[i]=sumfyA[i]+coef2*(yA[i]-yforce);
				sumfzA[i]=sumfzA[i]+coef2*(zA[i]-zforce);					
				
				sumfxB[j]=sumfxB[j]-coef2*(xA[i]-xforce);
				sumfyB[j]=sumfyB[j]-coef2*(yA[i]-yforce);
				sumfzB[j]=sumfzB[j]-coef2*(zA[i]-zforce);	
				
				
			sumV_ts = sumV_ts +sigma12/mysqrt12 - sigma6/mysqrt6;
				
			}
			
			
		}
	}

			


}