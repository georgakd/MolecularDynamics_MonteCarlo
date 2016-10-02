/*
 *  tersoffopt.c
 *  tersoff calculations
 *
 *  Created by Dimitra Georgakaki on 17/01/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */


#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void tersoffopt(double *x,double *y,double *z);
double xi_define(int myZ1, int myZ2);


void tersoffopt(double *x,double *y,double *z){
	double frep,fattr,fcut,theta,gfunc,zeta,fcut_ik,bond,bondvar,Vtersoff,Repterm,Attrterm,paronomastis;
	double S_mean, A_mean, B_mean, R_mean, lamda_mean, mi_mean; 
	
	
	sumenergy=0; // for the calculation of the cohesive energy

	
	for (i=0;i<npart;++i)
		
	{
		
		for (j=0;j<npart;++j)
		{	
			
			
			if (i!=j)  //because we do not have ij=ji like in lennard-jones
			{
				xforce=x[j];
				yforce=y[j];
				zforce=z[j]; // PBC for rij
				if (x[j]-x[i]>0.5*Lx) {xforce=x[j]-Lx;}
				if (x[j]-x[i]<-0.5*Lx) {xforce=x[j]+Lx;}
				if (y[j]-y[i]>0.5*Ly) {yforce=y[j]-Ly;}
				if (y[j]-y[i]<-0.5*Ly) {yforce=y[j]+Ly;}
				
				if (z[j]-z[i]>0.5*Lz) {zforce=z[j]-Lz;}
				if (z[j]-z[i]<-0.5*Lz) {zforce=z[j]+Lz;}
				
				r2=(xforce-x[i])*(xforce-x[i])+(yforce-y[i])*(yforce-y[i])+(zforce-z[i])*(zforce-z[i]);
				r=sqrt(r2); // rij calculation
				
				
				S_mean=sqrt(myS[i]*myS[j]); // Smean
				
	// main cut-off interval r<S
			if (r<=S_mean){
				
				
				A_mean=sqrt(myA[i]*myA[j]);// various means
				B_mean=sqrt(myB[i]*myB[j]);
				R_mean=sqrt(myR[i]*myR[j]);
				lamda_mean=(lamda[i]+lamda[j])/2;
				mi_mean=(mi[i]+mi[j])/2;
				frep=A_mean*exp(-lamda_mean*r);	
				fattr=-B_mean*exp(-mi_mean*r);

					
		// medium cut-off interval r<S	&& r>R
					if (r<=S_mean && r>=R_mean){
						fcut=0.5 +0.5*cos((pi*(r-R_mean))/(S_mean-R_mean));						
					}	//end of r<S && r>R
					
		// short cut-off interval r<R
					if (r<R_mean){
						fcut=1;					
					} //end of r<R
				

		// start the sum over k____________
				zeta=0; // initialization
				
				for(k=0;k<npart;++k)
				{
					if (k!=i && k!=j)
					{
						xikforce=x[k];yikforce=y[k];zikforce=z[k]; // PBC for rik
						if (x[k]-x[i]>0.5*Lx) {xikforce=x[k]-Lx;}
						if (x[k]-x[i]<-0.5*Lx) {xikforce=x[k]+Lx;}
						if (y[k]-y[i]>0.5*Ly) {yikforce=y[k]-Ly;}
						if (y[k]-y[i]<-0.5*Ly) {yikforce=y[k]+Ly;}
						
						if (z[k]-z[i]>0.5*Lz) {zikforce=z[k]-Lz;}
						if (z[k]-z[i]<-0.5*Lz) {zikforce=z[k]+Lz;}
						
						r2ik=(xikforce-x[i])*(xikforce-x[i])+(yikforce-y[i])*(yikforce-y[i])+(zikforce-z[i])*(zikforce-z[i]);
						rik=sqrt(r2ik); // rik calculation
					
						S_mean=sqrt(myS[i]*myS[k]);
						
				// main cut-off interval rik<S		
						if (rik<=S_mean){
							R_mean=sqrt(myR[i]*myR[k]);
							
							if (rik<=S_mean && rik>=R_mean){
								fcut_ik = 0.5 +0.5*cos(pi*(rik-R_mean)/(S_mean-R_mean));
							
							}//end of rik<S && rik >R
						
				// short cut-off interval r<R
						if (rik<R_mean){
								fcut_ik=1;					
							} //end of rik<R	
						
				// calculate cos, g, z and derivatives
							theta = ((xikforce-x[i])*(xforce-x[i])+(yikforce-y[i])*(yforce-y[i])+(zikforce-z[i])*(zforce-z[i]))/(r*rik);
							
							paronomastis=(di[i]*di[i]) + (hi[i]-theta)*(hi[i]-theta);
							gfunc = 1 + (ci[i]*ci[i])/(di[i]*di[i]) - (ci[i]*ci[i])/paronomastis;
												
							zeta = zeta + fcut_ik*gfunc;
						
							
											
						} // end of rik<S
					} // end of if
	
					
				} // end of k loop
				
				bondvar = (1+pow(beta[i],ni[i])*pow(zeta,ni[i]));
				xi = xi_define(Z[i], Z[j]);
				bond = xi*pow(bondvar,-1/(2*ni[i]));
				
				
				Repterm = fcut*frep;
				Attrterm = fcut*bond*fattr;				
				Vtersoff = Repterm + Attrterm;
				
				sumenergy = sumenergy+Vtersoff;
				
				} //end of r<S
				

							
			} // end if i#j
		} // end j
	} //end i
	
		


} //end of function
