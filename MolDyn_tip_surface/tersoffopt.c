/*
 *  tersoffopt.c
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


void tersoffopt(double *x,double *y,double *z);
double xi_define(int myZ1, int myZ2);


void tersoffopt(double *x,double *y,double *z){
	
	
	
	// Tersoff forces and derivatives
	
	double fcut,theta,gfunc,zeta,fcut_ik,frep,fattr,bond,bondvar,paronomastis;	
	double dxr,dyr,dzr,dxrik,dyrik,dzrik,dxfrep,dyfrep,dzfrep,dxfattr,dyfattr,dzfattr;
	double dxtheta,dytheta,dztheta,dxgfunc,dygfunc,dzgfunc,dxfcut,dyfcut,dzfcut,dxfcut_ik,dyfcut_ik,dzfcut_ik;
	double dxzeta,dyzeta,dzzeta,bondvar2,dxbond,dybond,dzbond,d_bij_factor,g_factor,fcut_factor,fcut_ik_factor;	
	double dx_zeta_II,dy_zeta_II,dz_zeta_II,rev,dx_bij_II,dy_bij_II,dz_bij_II; //new terms in bond order!!!	
	
	double S_mean, A_mean, B_mean, R_mean, lamda_mean, mi_mean; 
	
	sumenergy=0; // for the calculation of the cohesive energy
		
	for (k=0;k<npartB;++k)
	{sumfx[k]=0;sumfy[k]=0;sumfz[k]=0;}
	
	for (i=0;i<npartB;++i)		
	{
		for (j=0;j<npartB;++j)
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
				frep = A_mean*exp(-lamda_mean*r);	
				fattr = -B_mean*exp(-mi_mean*r);
				
				//distance derivatives
				dxr = (x[i]-xforce)/r;
				dyr = (y[i]-yforce)/r;
				dzr = (z[i]-zforce)/r;
				// other derivatives
				dxfrep = -lamda_mean*frep*dxr;
				dyfrep = -lamda_mean*frep*dyr;
				dzfrep = -lamda_mean*frep*dzr;
				dxfattr = -mi_mean*fattr*dxr;
				dyfattr = -mi_mean*fattr*dyr;
				dzfattr = -mi_mean*fattr*dzr;
					
		// medium cut-off interval r<S	&& r>R
					if (r<=S_mean && r>=R_mean){
						
						fcut = 0.5 +0.5*cos((pi*(r-R_mean))/(S_mean-R_mean));
						fcut_factor = -0.5* (pi/(S_mean-R_mean)) * sin((pi*(r-R_mean))/(S_mean-R_mean)); //<----------new!
						dxfcut =  fcut_factor*dxr;
						dyfcut =  fcut_factor*dyr;
						dzfcut =  fcut_factor*dzr;
					}	//end of r<S && r>R
					
		// short cut-off interval r<R
					if (r<R_mean){
						fcut=1; dxfcut=0; dyfcut=0; dzfcut=0;					
					} //end of r<R
				

		// start the sum over k____________
				zeta = 0; dxzeta = 0; dyzeta = 0; dzzeta = 0; // initialization
				dx_zeta_II=0;dy_zeta_II=0;dz_zeta_II=0;	//initialize new terms
				
				for(k=0;k<npartB;++k)
				{
					dx_zeta_III[k]=0;	
					dy_zeta_III[k]=0;
					dz_zeta_III[k]=0;
					
					
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
					
				// distance derivatives
						dxrik = (x[i]-xikforce)/rik;
						dyrik = (y[i]-yikforce)/rik;
						dzrik = (z[i]-zikforce)/rik;
						S_mean=sqrt(myS[i]*myS[k]);
						
				// main cut-off interval rik<S		
						if (rik<=S_mean){
							R_mean=sqrt(myR[i]*myR[k]);
							
							if (rik<=S_mean && rik>=R_mean){
								fcut_ik = 0.5 +0.5*cos(pi*(rik-R_mean)/(S_mean-R_mean));
								fcut_ik_factor = -0.5* (pi/(S_mean-R_mean)) * sin((pi*(rik-R_mean))/(S_mean-R_mean)); //<-------- new
								dxfcut_ik = fcut_ik_factor*dxrik;
								dyfcut_ik = fcut_ik_factor*dyrik;
								dzfcut_ik = fcut_ik_factor*dzrik;
							
							}//end of rik<S && rik >R
						
				// short cut-off interval r<R
						if (rik<R_mean){
								fcut_ik=1; dxfcut_ik=0; dyfcut_ik=0; dzfcut_ik=0;					
							} //end of rik<R	
						
				// calculate cos, g, z and derivatives
							
							rev=1.0/(r*rik); //new term	
							theta = rev*((xikforce-x[i])*(xforce-x[i])+(yikforce-y[i])*(yforce-y[i])+(zikforce-z[i])*(zforce-z[i])); //<------new
							dxtheta = (2*x[i]-xforce-xikforce)*rev - theta *((x[i]-xikforce)/r2ik + (x[i]-xforce)/r2);  //<------new
							dytheta = (2*y[i]-yforce-yikforce)*rev - theta *((y[i]-yikforce)/r2ik + (y[i]-yforce)/r2);   //<------new
							dztheta = (2*z[i]-zforce-zikforce)*rev - theta *((z[i]-zikforce)/r2ik + (z[i]-zforce)/r2);   //<------new
							
							paronomastis = (di[i]*di[i]) + (hi[i]-theta)*(hi[i]-theta);
							gfunc = 1 + (ci[i]*ci[i])/(di[i]*di[i]) - (ci[i]*ci[i])/ paronomastis;
							dxgfunc = -dxtheta*((ci[i]*ci[i]*2*(hi[i]- theta))/( paronomastis*paronomastis ));
							dygfunc = -dytheta*((ci[i]*ci[i]*2*(hi[i]- theta))/( paronomastis*paronomastis ));
							dzgfunc = -dztheta*((ci[i]*ci[i]*2*(hi[i]- theta))/( paronomastis*paronomastis ));
							
							g_factor = -((ci[i]*ci[i]*2*(hi[i]- theta))/( paronomastis*paronomastis ));
							
							zeta = zeta + fcut_ik * gfunc;
							dxzeta = dxzeta + (dxfcut_ik * gfunc + fcut_ik *dxgfunc); //<------------ type I
							dyzeta = dyzeta + (dyfcut_ik * gfunc + fcut_ik *dygfunc);
							dzzeta = dzzeta + (dzfcut_ik * gfunc + fcut_ik *dzgfunc);
							
							dx_zeta_II = dx_zeta_II + fcut_ik*g_factor*rev*(xikforce-x[i] + theta*rik*dxr);	//<------------ type II
							dy_zeta_II = dy_zeta_II + fcut_ik*g_factor*rev*(yikforce-y[i] + theta*rik*dyr);
							dz_zeta_II = dz_zeta_II + fcut_ik*g_factor*rev*(zikforce-z[i] + theta*rik*dzr);
							
							dx_zeta_III[k] = -gfunc*dxfcut_ik + fcut_ik*g_factor*rev*(xforce-x[i]+theta*r*dxrik);	//<------------ type III
							dy_zeta_III[k] = -gfunc*dyfcut_ik + fcut_ik*g_factor*rev*(yforce-y[i]+theta*r*dyrik);
							dz_zeta_III[k] = -gfunc*dzfcut_ik + fcut_ik*g_factor*rev*(zforce-z[i]+theta*r*dzrik);
										
						} // end of rik<S
					} // end of if
	
					
				} // end of k loop
				
				bondvar = (1+pow(beta[i],ni[i])*pow(zeta,ni[i]));
				bondvar2 = pow(beta[i],ni[i])*pow(zeta,ni[i]-1); // attention if zeta=0 , we have 0 ^ n-1 thus 0 in a negative power!!!!
				xi = xi_define(Z_B[i], Z_B[j]);
				bond = xi*pow(bondvar,-1/(2*ni[i]));
				
				Repterm = fcut*frep;
				Attrterm = fcut*bond*fattr;				
				Vtersoff = Repterm + Attrterm;				
				sumenergy = sumenergy+Vtersoff;
	
				// force calculations
				if (zeta!=0)
				{
					
					d_bij_factor = -0.5*bond*(bondvar2/bondvar);
					dxbond = d_bij_factor*dxzeta;  	//<---------- type I
					dybond = d_bij_factor*dyzeta;
					dzbond = d_bij_factor*dzzeta;
					dx_bij_II = d_bij_factor*dx_zeta_II;	//<----------  type II
					dy_bij_II = d_bij_factor*dy_zeta_II;
					dz_bij_II = d_bij_factor*dz_zeta_II;
					
					for(k=0;k<npartB;++k)		
					{
						sumfx[k]=sumfx[k]-0.5*fcut*fattr*d_bij_factor*dx_zeta_III[k];
						sumfy[k]=sumfy[k]-0.5*fcut*fattr*d_bij_factor*dy_zeta_III[k];
						sumfz[k]=sumfz[k]-0.5*fcut*fattr*d_bij_factor*dz_zeta_III[k];
					}
					
					
					
				}
				else if (zeta==0)
				{
					bondvar2 = 0; dxbond=0; dybond=0; dzbond=0;
					
					dx_bij_II=0;	
					dy_bij_II=0;
					dz_bij_II=0;
					
				}

				
				sumfx[i]=sumfx[i]-0.5*(dxfcut*(frep + bond*fattr) + fcut*(dxfrep+bond*dxfattr + fattr*dxbond));
				sumfy[i]=sumfy[i]-0.5*(dyfcut*(frep + bond*fattr) + fcut*(dyfrep+bond*dyfattr + fattr*dybond));
				sumfz[i]=sumfz[i]-0.5*(dzfcut*(frep + bond*fattr) + fcut*(dzfrep+bond*dzfattr + fattr*dzbond));
				
				sumfx[j]=sumfx[j]-0.5*(-dxfcut*(frep + bond*fattr) + fcut*(-dxfrep-bond*dxfattr + fattr*dx_bij_II));
				sumfy[j]=sumfy[j]-0.5*(-dyfcut*(frep + bond*fattr) + fcut*(-dyfrep-bond*dyfattr + fattr*dy_bij_II));
				sumfz[j]=sumfz[j]-0.5*(-dzfcut*(frep + bond*fattr) + fcut*(-dzfrep-bond*dzfattr + fattr*dz_bij_II));
				
				
				} //end of r<S
				
				
				
//--------------------------------------------start the rdfsum calculation for pairs i-j---------------------
				
				if (n>= N-rdf_avg)  //for rdf of many time steps, n>=N-rdf_avg
				{
					myindex=(int)floor(r/dr); // see the art of md p.92
					
					if (myindex<bins)
					{
						h[myindex] = h[myindex]+1; // +2 when lj
						
					}
					
				}
				
				// for matlab				
				if (monitor==1){
					
					myindex=(int)floor(r/dr); // see the art of md p.92
					
					if (myindex<bins)
					{
						hinst[myindex] = hinst[myindex]+1; // +2 when lj
						
					}
					
				}
				
//--------------------------------------------end the rdf calculation-----------------------------------------
							
			} // end if i#j
			
		} // end j
	} //end i
	
} //end of function