/*
 *  tersoff_forces_cells.c
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


void tersoff_forces_cells(double *x,double *y,double *z);
double xi_define(int myZ1, int myZ2);

//----------------------------------------- Interactions calculations with Linked Cells------------------------------------------------------------------------

void tersoff_forces_cells(double *x,double *y,double *z)
{
	double fcut,theta,gfunc,zeta,fcut_ik,frep,fattr,bond,bondvar,paronomastis;
	
	double dxr,dyr,dzr,dxrik,dyrik,dzrik,dxfrep,dyfrep,dzfrep,dxfattr,dyfattr,dzfattr;
	double dxtheta,dytheta,dztheta,dxgfunc,dygfunc,dzgfunc,dxfcut,dyfcut,dzfcut,dxfcut_ik,dyfcut_ik,dzfcut_ik;
	double dxzeta,dyzeta,dzzeta,bondvar2,dxbond,dybond,dzbond,d_bij_factor,g_factor,fcut_factor,fcut_ik_factor;
	
	double dx_zeta_II,dy_zeta_II,dz_zeta_II,rev,dx_bij_II,dy_bij_II,dz_bij_II; //new terms in bond order!!!		
	double S_mean, A_mean, B_mean, R_mean, lamda_mean, mi_mean;  
	int icell, cellpart, allcellpart, neighborcellpart;
	int HEADp,LISTp, locali, localj,localk;
	
	
	
	//STEP 1: HEAD - LIST---------------------------------------------------------------------------------------------------------
	
	
	for(cell=0;cell<cellsx*cellsy*cellsz;++cell)
	{
		HEAD[cell] = -1; //because 0 is the index of the first particle!!!!!!!
		
	}
	for(l=0;l<npartB;++l)
	{
		LIST[l] = -1; //because 0 is the index of the first particle!!!!!!!
		
	}
	
	
	for(l=0;l<npartB;++l)
	{
		icell = 
		(int)(floor((x[l] - xmin)/Lcx) + floor((y[l] - ymin)/Lcy)*cellsx + floor((z[l] - zmin)/Lcz)*cellsx*cellsy); // because index paricle starts from 0 we do not add +1 like in mathematica
		
		
		LIST[l] = HEAD[icell];
		HEAD[icell] = l; 
		
	}
//---------------------------Initialization--------------------------------------

		
	for (k=0;k<npartB;++k)
	{sumfx[k]=0;sumfy[k]=0;sumfz[k]=0;}
		
		
	sumenergy = 0; // initialize for calculation of the potential energy when we use cells

	
	
//STEP 2: Populate cellpartarray using the HEAD-LIST reading algorithm (*the cellpart variable counts the number of particles in the central cell*)	
for (cell=0;cell<cellsx*cellsy*cellsz;++cell){	
	   
	   cellpart = 0;
	   HEADp = HEAD[cell];
	if(HEADp != -1){
		  
		  cellpartarray[cellpart] = HEADp;
		  cellpart = cellpart + 1;
		  
		  LISTp = LIST[HEADp];
		while(LISTp != -1){
			if(LISTp != -1){
				
				   cellpartarray[cellpart] = LISTp;
					cellpart = cellpart + 1;
			LISTp = LIST[LISTp];}}}

//------------------------------------------------------------------------------------------------------------------------------
//STEP 3:	Populate the allcellpart array - the first elements are the same with cellpartarray
	allcellpart = 0;
	 HEADp = HEAD[cell];
	if(HEADp != -1){
		allpartarray[allcellpart] = HEADp;
	    allcellpart = allcellpart + 1;
	   
	   LISTp = LIST[HEADp];
	   	while(LISTp != -1){
			if(LISTp != -1){
				
				allpartarray[allcellpart] = LISTp;
				allcellpart = allcellpart + 1;
				
				LISTp = LIST[LISTp];}}}
	
//------------------------------------------------------------------------------------------------------------------------------
//STEP 4: loop over all neighboring cells and store particles indices in allpartarray, continuing from where the previous algorithm stopped
	
	for (next=0;next<26;++next){
	   HEADp = HEAD[neighbours[next][cell]];
	  if(HEADp != -1){
		  
		  allpartarray[allcellpart] = HEADp;
		  allcellpart = allcellpart + 1;
		  
		  LISTp = LIST[HEADp];
		  while(LISTp != -1){
			  if(LISTp != -1){
				
				 allpartarray[allcellpart] = LISTp;
				  
				  allcellpart = allcellpart + 1;
				  
				  LISTp = LIST[LISTp];}}}
				   
	} // end of loop over neighbors
	
	
//------------------------------------------------------------------------------------------------------------------------------
	
//STEP 5: loop over the elements of cellpartarray in order to calculate the interactions in the central cell

	
	//__________________start the calculation of the potential + forces_______________
	
	if (cellpart!=0){
		
		for (i=0;i<cellpart;++i)
		{
	
			
			for (j=0;j<cellpart;++j)
			{	
				
				if (i!=j)  //because we do not have ij=ji like in lennard-jones
				{
					
					locali=cellpartarray[i];
					localj=cellpartarray[j];
					xforce=x[localj];
					yforce=y[localj];
					zforce=z[localj]; // PBC for rij
					if (x[localj]-x[locali]>0.5*Lx) {xforce=x[localj]-Lx;}
					if (x[localj]-x[locali]<-0.5*Lx) {xforce=x[localj]+Lx;}
					if (y[localj]-y[locali]>0.5*Ly) {yforce=y[localj]-Ly;}
					if (y[localj]-y[locali]<-0.5*Ly) {yforce=y[localj]+Ly;}
					
					if (z[localj]-z[locali]>0.5*Lz) {zforce=z[localj]-Lz;}
					if (z[localj]-z[locali]<-0.5*Lz) {zforce=z[localj]+Lz;}
					
					r2=(xforce-x[locali])*(xforce-x[locali])+(yforce-y[locali])*(yforce-y[locali])+(zforce-z[locali])*(zforce-z[locali]);
					r=sqrt(r2); // rij calculation
					S_mean=sqrt(myS[locali]*myS[localj]);
					
					// main cut-off interval r<S
					if (r<=S_mean){
						
						A_mean=sqrt(myA[locali]*myA[localj]);// various means
						B_mean=sqrt(myB[locali]*myB[localj]);
						R_mean=sqrt(myR[locali]*myR[localj]);
						lamda_mean=(lamda[locali]+lamda[localj])/2;
						mi_mean=(mi[locali]+mi[localj])/2;
						
						frep=A_mean*exp(-lamda_mean*r);	
						fattr=-B_mean*exp(-mi_mean*r);
						//distance derivatives
						dxr = (x[locali]-xforce)/r; 
						dyr = (y[locali]-yforce)/r;
						dzr = (z[locali]-zforce)/r;
						// other derivatives
						dxfrep = -lamda_mean*frep*dxr;
						dyfrep = -lamda_mean*frep*dyr;
						dzfrep = -lamda_mean*frep*dzr;
						dxfattr = -mi_mean*fattr*dxr;
						dyfattr = -mi_mean*fattr*dyr;
						dzfattr = -mi_mean*fattr*dzr;
						
						// medium cut-off interval r<S	&& r>R
						if (r<=S_mean && r>=R_mean){
							fcut=0.5 +0.5*cos((pi*(r-R_mean))/(S_mean-R_mean));	
							fcut_factor = -0.5* (pi/(S_mean-R_mean)) * sin((pi*(r-R_mean))/(S_mean-R_mean)); //<----------new!
							dxfcut =  fcut_factor*dxr;
							dyfcut =  fcut_factor*dyr;
							dzfcut =  fcut_factor*dzr;
						}	//end of r<S && r>R
						
						// short cut-off interval r<R
						if (r<R_mean){
							fcut=1;dxfcut=0;dyfcut=0;dzfcut=0;					
						} //end of r<R
						
						
						// start the sum over k____________
						zeta=0; dxzeta=0; dyzeta=0; dzzeta=0; // initialization
						dx_zeta_II=0;dy_zeta_II=0;dz_zeta_II=0;	//initialize new terms
						
						
						
						for(k=0;k<allcellpart;++k)
						{	
							localk=allpartarray[k];
							dx_zeta_III[localk]=0;		
							dy_zeta_III[localk]=0;
							dz_zeta_III[localk]=0;
							
							if (localk!=locali && localk!=localj)
							{	
								
								xikforce=x[localk];yikforce=y[localk];zikforce=z[localk]; // PBC for rik
								if (x[localk]-x[locali]>0.5*Lx) {xikforce=x[localk]-Lx;}
								if (x[localk]-x[locali]<-0.5*Lx) {xikforce=x[localk]+Lx;}
								if (y[localk]-y[locali]>0.5*Ly) {yikforce=y[localk]-Ly;}
								if (y[localk]-y[locali]<-0.5*Ly) {yikforce=y[localk]+Ly;}
								
								if (z[localk]-z[locali]>0.5*Lz) {zikforce=z[localk]-Lz;}
								if (z[localk]-z[locali]<-0.5*Lz) {zikforce=z[localk]+Lz;}
								
								r2ik=(xikforce-x[locali])*(xikforce-x[locali])+(yikforce-y[locali])*(yikforce-y[locali])+(zikforce-z[locali])*(zikforce-z[locali]);
								rik=sqrt(r2ik); // rik calculation
								
								// distance derivatives
								dxrik = (x[locali]-xikforce)/rik;
								dyrik = (y[locali]-yikforce)/rik;
								dzrik = (z[locali]-zikforce)/rik;
								S_mean=sqrt(myS[locali]*myS[localk]);
								
								// main cut-off interval rik<S		
								if (rik<=S_mean){
									
									R_mean=sqrt(myR[locali]*myR[localk]);
									if (rik<=S_mean && rik>=R_mean){
										fcut_ik = 0.5 +0.5*cos(pi*(rik-R_mean)/(S_mean-R_mean));
										fcut_ik_factor = -0.5* (pi/(S_mean-R_mean)) * sin((pi*(rik-R_mean))/(S_mean-R_mean)); //<-------- new
										dxfcut_ik = fcut_ik_factor*dxrik;
										dyfcut_ik = fcut_ik_factor*dyrik;
										dzfcut_ik = fcut_ik_factor*dzrik;
										
									}//end of rik<S && rik >R
									
									// short cut-off interval r<R
									if (rik<R_mean){
										fcut_ik=1;dxfcut_ik=0;dyfcut_ik=0;dzfcut_ik=0;					
									} //end of rik<R	
									
									// calculate cos, g, z and derivatives
									rev=1.0/(r*rik); //new term	
									theta = rev*((xikforce-x[locali])*(xforce-x[locali])+(yikforce-y[locali])*(yforce-y[locali])+(zikforce-z[locali])*(zforce-z[locali])); //<------new
									dxtheta = (2*x[locali]-xforce-xikforce)*rev - theta *((x[locali]-xikforce)/r2ik + (x[locali]-xforce)/r2);  //<------new
									dytheta = (2*y[locali]-yforce-yikforce)*rev - theta *((y[locali]-yikforce)/r2ik + (y[locali]-yforce)/r2);   //<------new
									dztheta = (2*z[locali]-zforce-zikforce)*rev - theta *((z[locali]-zikforce)/r2ik + (z[locali]-zforce)/r2);   //<------new
									
									paronomastis=(di[locali]*di[locali]) + (hi[locali]-theta)*(hi[locali]-theta);
									gfunc = 1 + (ci[locali]*ci[locali])/(di[locali]*di[locali]) - (ci[locali]*ci[locali])/paronomastis;
									dxgfunc = -dxtheta*((ci[locali]*ci[locali]*2*(hi[locali]-theta))/(paronomastis*paronomastis));
									dygfunc = -dytheta*((ci[locali]*ci[locali]*2*(hi[locali]-theta))/(paronomastis*paronomastis));
									dzgfunc = -dztheta*((ci[locali]*ci[locali]*2*(hi[locali]-theta))/(paronomastis*paronomastis));
									
									g_factor = -((ci[locali]*ci[locali]*2*(hi[locali]- theta))/( paronomastis*paronomastis ));
									
									
									zeta = zeta + fcut_ik*gfunc;
									dxzeta = dxzeta + (dxfcut_ik*gfunc + fcut_ik*dxgfunc);
									dyzeta = dyzeta + (dyfcut_ik*gfunc + fcut_ik*dygfunc);
									dzzeta = dzzeta + (dzfcut_ik*gfunc + fcut_ik*dzgfunc);
									
									dx_zeta_II = dx_zeta_II + fcut_ik*g_factor*rev*(xikforce-x[locali] + theta*rik*dxr);	//<------------ type II
									dy_zeta_II = dy_zeta_II + fcut_ik*g_factor*rev*(yikforce-y[locali] + theta*rik*dyr);
									dz_zeta_II = dz_zeta_II + fcut_ik*g_factor*rev*(zikforce-z[locali] + theta*rik*dzr);
									
									dx_zeta_III[localk] = -gfunc*dxfcut_ik + fcut_ik*g_factor*rev*(xforce-x[locali]+theta*r*dxrik);	//<------------ type III
									dy_zeta_III[localk] = -gfunc*dyfcut_ik + fcut_ik*g_factor*rev*(yforce-y[locali]+theta*r*dyrik);
									dz_zeta_III[localk] = -gfunc*dzfcut_ik + fcut_ik*g_factor*rev*(zforce-z[locali]+theta*r*dzrik);
					
									
								} // end of rik<S
							} // end of if
						} // end of k loop
						
						bondvar = (1+pow(beta[locali],ni[locali])*pow(zeta,ni[locali]));
						bondvar2 = pow(beta[i],ni[i])*pow(zeta,ni[i]-1); // attention if zeta=0 , we have 0 ^ n-1 thus 0 in a negative power!!!!
						xi = xi_define(Z_B[locali], Z_B[localj]);
						bond = xi*pow(bondvar,-1/(2*ni[locali]));
						
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
							
							for(k=0;k<allcellpart;++k)				
							{
								localk=allpartarray[k];
								sumfx[localk]=sumfx[localk]-0.5*fcut*fattr*d_bij_factor*dx_zeta_III[localk];
								sumfy[localk]=sumfy[localk]-0.5*fcut*fattr*d_bij_factor*dy_zeta_III[localk];
								sumfz[localk]=sumfz[localk]-0.5*fcut*fattr*d_bij_factor*dz_zeta_III[localk];
							}
							
						}
						else if (zeta==0)
						{
							bondvar2 = 0; dxbond=0; dybond=0; dzbond=0;
							
							dx_bij_II=0;									
							dy_bij_II=0;
							dz_bij_II=0;
							
						}
						
						
						sumfx[locali]=sumfx[locali]-0.5*(dxfcut*(frep + bond*fattr) + fcut*(dxfrep+bond*dxfattr + fattr*dxbond));		
						sumfy[locali]=sumfy[locali]-0.5*(dyfcut*(frep + bond*fattr) + fcut*(dyfrep+bond*dyfattr + fattr*dybond));
						sumfz[locali]=sumfz[locali]-0.5*(dzfcut*(frep + bond*fattr) + fcut*(dzfrep+bond*dzfattr + fattr*dzbond));
						
						sumfx[localj]=sumfx[localj]-0.5*(-dxfcut*(frep + bond*fattr) + fcut*(-dxfrep-bond*dxfattr + fattr*dx_bij_II));
						sumfy[localj]=sumfy[localj]-0.5*(-dyfcut*(frep + bond*fattr) + fcut*(-dyfrep-bond*dyfattr + fattr*dy_bij_II));
						sumfz[localj]=sumfz[localj]-0.5*(-dzfcut*(frep + bond*fattr) + fcut*(-dzfrep-bond*dzfattr + fattr*dz_bij_II));
							
						
					} //end of r<S
					
//--------------------------------------------start the rdfsum calculation for pairs i-j---------------------
					
					if (n>= N-rdf_avg && status[locali]==1 && status[localj]==1)  //for rdf of many time steps, n>=N-rdf_avg
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
	} // end of if cellpart !=0
	
		
// STEP 6: (populate neighborcellpartarray using the HEAD-LIST reading algorithm) 	
	
	for(next=0;next<26;++next){
	   
	   neighborcellpart = 0;//(*use different counter*)
	   
	   HEADp = HEAD[neighbours[next][cell]];
		if (HEADp != -1){
		 
		  neighborcellpartarray[neighborcellpart] = HEADp;
			 neighborcellpart = neighborcellpart + 1;
			
		  LISTp = LIST[HEADp];
			while(LISTp != -1){
				if(LISTp != -1){
				   
				   neighborcellpartarray[neighborcellpart] = LISTp;
					neighborcellpart = neighborcellpart + 1;
				LISTp = LIST[LISTp];}}}
		
// STEP 7: FINAL INTERACTIONS
		
		
		//_________________start the calculation of the potential + forces_______________
		
		
		if (neighborcellpart!=0){
			
			for (i=0;i<cellpart;++i)
			{
							
				for (j=0;j<neighborcellpart;++j)
				{				
					
					//if (i!=j)  //because we do not have ij=ji like in lennard-jones
					//{	
						
						locali=cellpartarray[i];
						localj=neighborcellpartarray[j];
						xforce=x[localj];
						yforce=y[localj];
						zforce=z[localj]; // PBC for rij
						if (x[localj]-x[locali]>0.5*Lx) {xforce=x[localj]-Lx;}
						if (x[localj]-x[locali]<-0.5*Lx) {xforce=x[localj]+Lx;}
						if (y[localj]-y[locali]>0.5*Ly) {yforce=y[localj]-Ly;}
						if (y[localj]-y[locali]<-0.5*Ly) {yforce=y[localj]+Ly;}
						
						if (z[localj]-z[locali]>0.5*Lz) {zforce=z[localj]-Lz;}
						if (z[localj]-z[locali]<-0.5*Lz) {zforce=z[localj]+Lz;}
					r2=(xforce-x[locali])*(xforce-x[locali])+(yforce-y[locali])*(yforce-y[locali])+(zforce-z[locali])*(zforce-z[locali]);
					r=sqrt(r2); // rij calculation
					S_mean=sqrt(myS[locali]*myS[localj]);
					
					// main cut-off interval r<S
					if (r<=S_mean){
						
						A_mean=sqrt(myA[locali]*myA[localj]);// various means
						B_mean=sqrt(myB[locali]*myB[localj]);
						R_mean=sqrt(myR[locali]*myR[localj]);
						lamda_mean=(lamda[locali]+lamda[localj])/2;
						mi_mean=(mi[locali]+mi[localj])/2;
						
						frep=A_mean*exp(-lamda_mean*r);	
						fattr=-B_mean*exp(-mi_mean*r);
						//distance derivatives
						dxr = (x[locali]-xforce)/r; 
						dyr = (y[locali]-yforce)/r;
						dzr = (z[locali]-zforce)/r;
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
							fcut=1;dxfcut=0;dyfcut=0;dzfcut=0;					
						} //end of r<R
						
						
						// start the sum over k____________
						zeta=0; dxzeta=0; dyzeta=0; dzzeta=0; // initialization
						dx_zeta_II=0;dy_zeta_II=0;dz_zeta_II=0;	//initialize new terms
						
						
						for(k=0;k<allcellpart;++k)
						{	
							localk=allpartarray[k];
							dx_zeta_III[localk]=0;					
							dy_zeta_III[localk]=0;
							dz_zeta_III[localk]=0;
							
							if (localk!=locali && localk!=localj)
							{	
								
								xikforce=x[localk];yikforce=y[localk];zikforce=z[localk]; // PBC for rik
								if (x[localk]-x[locali]>0.5*Lx) {xikforce=x[localk]-Lx;}
								if (x[localk]-x[locali]<-0.5*Lx) {xikforce=x[localk]+Lx;}
								if (y[localk]-y[locali]>0.5*Ly) {yikforce=y[localk]-Ly;}
								if (y[localk]-y[locali]<-0.5*Ly) {yikforce=y[localk]+Ly;}
								
								if (z[localk]-z[locali]>0.5*Lz) {zikforce=z[localk]-Lz;}
								if (z[localk]-z[locali]<-0.5*Lz) {zikforce=z[localk]+Lz;}
								
								r2ik=(xikforce-x[locali])*(xikforce-x[locali])+(yikforce-y[locali])*(yikforce-y[locali])+(zikforce-z[locali])*(zikforce-z[locali]);
								rik=sqrt(r2ik); // rik calculation
								
								// distance derivatives
								dxrik = (x[locali]-xikforce)/rik;
								dyrik = (y[locali]-yikforce)/rik;
								dzrik = (z[locali]-zikforce)/rik;
								S_mean=sqrt(myS[locali]*myS[localk]);
								
								// main cut-off interval rik<S		
								if (rik<=S_mean){
									
									R_mean=sqrt(myR[locali]*myR[localk]);
									if (rik<=S_mean && rik>=R_mean){
										fcut_ik = 0.5 +0.5*cos(pi*(rik-R_mean)/(S_mean-R_mean));
										fcut_ik_factor = -0.5* (pi/(S_mean-R_mean)) * sin((pi*(rik-R_mean))/(S_mean-R_mean)); //<-------- new
										dxfcut_ik = fcut_ik_factor*dxrik;
										dyfcut_ik = fcut_ik_factor*dyrik;
										dzfcut_ik = fcut_ik_factor*dzrik;
										
									}//end of rik<S && rik >R
									
									// short cut-off interval r<R
									if (rik<R_mean){
										fcut_ik=1;dxfcut_ik=0;dyfcut_ik=0;dzfcut_ik=0;					
									} //end of rik<R	
									
									// calculate cos, g, z and derivatives
									rev=1.0/(r*rik); //new term	
									theta = rev*((xikforce-x[locali])*(xforce-x[locali])+(yikforce-y[locali])*(yforce-y[locali])+(zikforce-z[locali])*(zforce-z[locali])); //<------new
									dxtheta = (2*x[locali]-xforce-xikforce)*rev - theta *((x[locali]-xikforce)/r2ik + (x[locali]-xforce)/r2);  //<------new
									dytheta = (2*y[locali]-yforce-yikforce)*rev - theta *((y[locali]-yikforce)/r2ik + (y[locali]-yforce)/r2);   //<------new
									dztheta = (2*z[locali]-zforce-zikforce)*rev - theta *((z[locali]-zikforce)/r2ik + (z[locali]-zforce)/r2);   //<------new
									
									paronomastis=(di[locali]*di[locali]) + (hi[locali]-theta)*(hi[locali]-theta);
									gfunc = 1 + (ci[locali]*ci[locali])/(di[locali]*di[locali]) - (ci[locali]*ci[locali])/paronomastis;
									dxgfunc = -dxtheta*((ci[locali]*ci[locali]*2*(hi[locali]-theta))/(paronomastis*paronomastis));
									dygfunc = -dytheta*((ci[locali]*ci[locali]*2*(hi[locali]-theta))/(paronomastis*paronomastis));
									dzgfunc = -dztheta*((ci[locali]*ci[locali]*2*(hi[locali]-theta))/(paronomastis*paronomastis));
									
									g_factor = -((ci[locali]*ci[locali]*2*(hi[locali]- theta))/( paronomastis*paronomastis ));
																		
									zeta = zeta + fcut_ik*gfunc;
									dxzeta = dxzeta + (dxfcut_ik*gfunc + fcut_ik*dxgfunc);
									dyzeta = dyzeta + (dyfcut_ik*gfunc + fcut_ik*dygfunc);
									dzzeta = dzzeta + (dzfcut_ik*gfunc + fcut_ik*dzgfunc);
									
									dx_zeta_II = dx_zeta_II + fcut_ik*g_factor*rev*(xikforce-x[locali] + theta*rik*dxr);	//<------------ type II
									dy_zeta_II = dy_zeta_II + fcut_ik*g_factor*rev*(yikforce-y[locali] + theta*rik*dyr);
									dz_zeta_II = dz_zeta_II + fcut_ik*g_factor*rev*(zikforce-z[locali] + theta*rik*dzr);
									
									dx_zeta_III[localk] = -gfunc*dxfcut_ik + fcut_ik*g_factor*rev*(xforce-x[locali]+theta*r*dxrik);	//<------------ type III
									dy_zeta_III[localk] = -gfunc*dyfcut_ik + fcut_ik*g_factor*rev*(yforce-y[locali]+theta*r*dyrik);
									dz_zeta_III[localk] = -gfunc*dzfcut_ik + fcut_ik*g_factor*rev*(zforce-z[locali]+theta*r*dzrik);
									
									
								} // end of rik<S
							} // end of if
						} // end of k loop
						
						bondvar = (1+pow(beta[locali],ni[locali])*pow(zeta,ni[locali]));
						bondvar2 = pow(beta[i],ni[i])*pow(zeta,ni[i]-1); // attention if zeta=0 , we have 0 ^ n-1 thus 0 in a negative power!!!!
						xi = xi_define(Z_B[locali], Z_B[localj]);
						bond = xi*pow(bondvar,-1/(2*ni[locali]));
						
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
							
							for(k=0;k<allcellpart;++k)	
							{
								localk=allpartarray[k];
								sumfx[localk]=sumfx[localk]-0.5*fcut*fattr*d_bij_factor*dx_zeta_III[localk];
								sumfy[localk]=sumfy[localk]-0.5*fcut*fattr*d_bij_factor*dy_zeta_III[localk];
								sumfz[localk]=sumfz[localk]-0.5*fcut*fattr*d_bij_factor*dz_zeta_III[localk];
							}
							
						}
						else if (zeta==0)
						{
							bondvar2 = 0; dxbond=0; dybond=0; dzbond=0;
							
							dx_bij_II=0;								
							dy_bij_II=0;
							dz_bij_II=0;
							
						}
						
						
						sumfx[locali]=sumfx[locali]-0.5*(dxfcut*(frep + bond*fattr) + fcut*(dxfrep+bond*dxfattr + fattr*dxbond));		
						sumfy[locali]=sumfy[locali]-0.5*(dyfcut*(frep + bond*fattr) + fcut*(dyfrep+bond*dyfattr + fattr*dybond));
						sumfz[locali]=sumfz[locali]-0.5*(dzfcut*(frep + bond*fattr) + fcut*(dzfrep+bond*dzfattr + fattr*dzbond));
						
						sumfx[localj]=sumfx[localj]-0.5*(-dxfcut*(frep + bond*fattr) + fcut*(-dxfrep-bond*dxfattr + fattr*dx_bij_II));	
						sumfy[localj]=sumfy[localj]-0.5*(-dyfcut*(frep + bond*fattr) + fcut*(-dyfrep-bond*dyfattr + fattr*dy_bij_II));
						sumfz[localj]=sumfz[localj]-0.5*(-dzfcut*(frep + bond*fattr) + fcut*(-dzfrep-bond*dzfattr + fattr*dz_bij_II));
						
					
							
						} //end of r<S
			
//--------------------------------------------start the rdfsum calculation for pairs i-j---------------------
					
					if (n>= N-rdf_avg && status[locali]==1 && status[localj]==1)  //for rdf of many time steps, n>=N-rdf_avg
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
					
			
						
					//} // end if i#j
				} // end j
			} //end i
			
		} //end of if neighborcellpart!=0
		


	} // end of steps 6-7

}  // end of loop over all cells

} // end of function