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
	/////////////forces/////////////////////////////////		  
	double coef1,coef2,sigma6,sigma12,sigma14,sigma8,mysqrt6,mysqrt8,mysqrt12,mysqrt14;
	
	for (k=0;k<npart;++k){sumfx[k]=0;sumfy[k]=0;sumfz[k]=0;}
	
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

				sumfx[i]=sumfx[i]+coef2*(x[i]-xforce);
				sumfy[i]=sumfy[i]+coef2*(y[i]-yforce);
				sumfz[i]=sumfz[i]+coef2*(z[i]-zforce);	
				
				sumfx[j]=sumfx[j]-coef2*(x[i]-xforce);
				sumfy[j]=sumfy[j]-coef2*(y[i]-yforce);
				sumfz[j]=sumfz[j]-coef2*(z[i]-zforce);	
				
				
			sumV=sumV+sigma12/mysqrt12-sigma6/mysqrt6;
			}
			
			//--------------------------------------------start the rdfsum calculation				
			if (n>= N-rdf_avg)
			{
				myindex=(int)floor(r/dr); // see the art of md p.92
				
				if (myindex<bins)
				{
					h[myindex] = h[myindex]+2;
					
				}
				
			}
			
// for matlab				
			if (monitor==1){
				
				myindex=(int)floor(r/dr); // see the art of md p.92
				
				if (myindex<bins)
				{
					hinst[myindex] = hinst[myindex]+2; // +2 when lj
					
				}
				
				
				
			}
			
			//--------------------------------------------end the rdf calculation

		


		} // end of j
	} // end of i
}