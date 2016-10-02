#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void tersoff_forces_cells(double *x,double *y,double *z);
double xi_define(int myZ1, int myZ2);

//----------------------------------------- Interactions calculations with Linked Cells------------------------------------------------------------------------

void tersoff_forces_cells(double *x,double *y,double *z)
{	
	double frep,fattr,fcut,theta,gfunc,zeta,fcut_ik,bond,bondvar,Vtersoff,Repterm,Attrterm,paronomastis;
	double S_mean, A_mean, B_mean, R_mean, lamda_mean, mi_mean; 

//STEP 1: HEAD - LIST---------------------------------------------------------------------------------------------------------
	
	
	for(cell=0;cell<totalcells;++cell)
	{
		HEAD[cell] = -1; //because 0 is the index of the first particle!!!!!!!
		
	}
	for(l=0;l<npart;++l)
	{
		LIST[l] = -1; //because 0 is the index of the first particle!!!!!!!
		
	}
	
	
	for(l=0;l<npart;++l)
	{
		icell = 
		(int)(floor((x[l] - xmin)/Lcx) + floor((y[l] - ymin)/Lcy)*cellsx + floor((z[l] - zmin)/Lcz)*cellsx*cellsy); // because index paricle starts from 0 we do not add +1 like in mathematica
		
		
		LIST[l] = HEAD[icell];
		HEAD[icell] = l; 
		
	}
//---------------------------Initialization--------------------------------------

		
		
	sumenergy = 0; // initialize for calculation of the potential energy when we use cells

	
//STEP 2: Populate cellpartarray using the HEAD-LIST reading algorithm (*the cellpart variable counts the number of particles in the central cell*)	
for (cell=0;cell<totalcells;++cell){	
	   
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
						
						for(k=0;k<allcellpart;++k)
						{	
							localk=allpartarray[k];
							
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
							
								S_mean=sqrt(myS[locali]*myS[localk]);
								
								// main cut-off interval rik<S		
								if (rik<=S_mean){
									
									R_mean=sqrt(myR[locali]*myR[localk]);
									if (rik<=S_mean && rik>=R_mean){
										fcut_ik = 0.5 +0.5*cos(pi*(rik-R_mean)/(S_mean-R_mean));
									
										
									}//end of rik<S && rik >R
									
									// short cut-off interval r<R
									if (rik<R_mean){
										fcut_ik=1;					
									} //end of rik<R	
									
									// calculate cos, g, z and derivatives
									theta = ((xikforce-x[locali])*(xforce-x[locali])+(yikforce-y[locali])*(yforce-y[locali])+(zikforce-z[locali])*(zforce-z[locali]))/(r*rik);
									
									paronomastis=(di[locali]*di[locali]) + (hi[locali]-theta)*(hi[locali]-theta);
									gfunc = 1 + (ci[locali]*ci[locali])/(di[locali]*di[locali]) - (ci[locali]*ci[locali])/paronomastis;
											
									zeta = zeta + fcut_ik*gfunc;
													
									
								} // end of rik<S
							} // end of if
						} // end of k loop
						
						bondvar = (1+pow(beta[locali],ni[locali])*pow(zeta,ni[locali]));
						xi = xi_define(Z[locali], Z[localj]);
						bond = xi*pow(bondvar,-1/(2*ni[locali]));
						
						Repterm = fcut*frep;
						Attrterm = fcut*bond*fattr;				
						Vtersoff = Repterm + Attrterm;
						
						sumenergy = sumenergy+Vtersoff;
			
						
					} //end of r<S



					
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
						
						// medium cut-off interval r<S	&& r>R
						if (r<=S_mean && r>=R_mean){
							fcut=0.5 +0.5*cos((pi*(r-R_mean))/(S_mean-R_mean));						
					
						}	//end of r<S && r>R
						
						// short cut-off interval r<R
						if (r<R_mean){
							fcut=1;					
						} //end of r<R
						
						
						// start the sum over k____________
						zeta=0;  // initialization
						
						for(k=0;k<allcellpart;++k)
						{	
							localk=allpartarray[k];
							
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
								
							
								S_mean=sqrt(myS[locali]*myS[localk]);
								
								// main cut-off interval rik<S		
								if (rik<=S_mean){
									
									R_mean=sqrt(myR[locali]*myR[localk]);
									if (rik<=S_mean && rik>=R_mean){
										fcut_ik = 0.5 +0.5*cos(pi*(rik-R_mean)/(S_mean-R_mean));
									
										
									}//end of rik<S && rik >R
									
									// short cut-off interval r<R
									if (rik<R_mean){
										fcut_ik=1;				
									} //end of rik<R	
									
									// calculate cos, g, z and derivatives
									theta = ((xikforce-x[locali])*(xforce-x[locali])+(yikforce-y[locali])*(yforce-y[locali])+(zikforce-z[locali])*(zforce-z[locali]))/(r*rik);
							
									
									paronomastis=(di[locali]*di[locali]) + (hi[locali]-theta)*(hi[locali]-theta);
									gfunc = 1 + (ci[locali]*ci[locali])/(di[locali]*di[locali]) - (ci[locali]*ci[locali])/paronomastis;

									
									zeta = zeta + fcut_ik*gfunc;
								
									
									
								} // end of rik<S
							} // end of if
						} // end of k loop
						
						bondvar = (1+pow(beta[locali],ni[locali])*pow(zeta,ni[locali]));
						xi = xi_define(Z[locali], Z[localj]);
						bond = xi*pow(bondvar,-1/(2*ni[locali]));
						
						Repterm = fcut*frep;
						Attrterm = fcut*bond*fattr;				
						Vtersoff = Repterm + Attrterm;
						
						sumenergy = sumenergy+Vtersoff;
							
						} //end of r<S


			
						
					//} // end if i#j
				} // end j
			} //end i
			
		} //end of if neighborcellpart!=0
		


	} // end of steps 6-7

}  // end of loop over all cells

} // end of void function
