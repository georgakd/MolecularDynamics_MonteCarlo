/*
 *  lj_forces_cells.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 31/01/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void lj_forces_cells(double *x,double *y,double *z);

//----------------------------------------- Interactions calculations with Linked Cells------------------------------------------------------------------------

void lj_forces_cells(double *x,double *y,double *z)
{	

	double coef1,coef2,sigma6,sigma12,sigma14,sigma8,mysqrt6,mysqrt8,mysqrt12,mysqrt14;
	
//------------ Build HEAD and LIST-------------------------------------
	
	
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
		(int)(floor((x[l] - xmin)/Lcx) + floor((y[l] - ymin)/Lcy)*cellsx + floor((z[l] - zmin)/Lcz)*cellsx*cellsy);
		
		
		LIST[l] = HEAD[icell];
		HEAD[icell] = l; 
		
	}
	
	
	
//--------------------Initialization--------------------//
	

	
	for (k=0;k<npart;++k){sumfx[k]=0;sumfy[k]=0;sumfz[k]=0;}
	
	sumVcells = 0; // initialize for calculation of the potential energy when we use cells 

//Part 1: Read current cell and define the index of the particles that it contains, do not print empty cells!!!!!

for (cell=0;cell<totalcells;++cell){		  
	
	counter = 0;
	HEADp = HEAD[cell];
	if(HEADp != -1){
		
		
		partindex[counter] = HEADp;
		
		counter = counter + 1;
		LISTp = LIST[HEADp];
				
		while(LISTp != -1){
			if(LISTp != -1){
				
				partindex[counter] = LISTp;

				counter = counter + 1;
				
				LISTp = LIST[LISTp];
				
			} //end of if LISTp
		} //end of while 
	} // end of if HEADp	
	
	
//Part 2: Calculate the particle interactions in the current cell 
	
	
	if(counter != 0){ //The loop must be calculated only if there are atoms in the cell we choose
		
		for(i=0;i<counter;++i){	
			for(j=0;j<counter;++j){
				if (i!=j){
					
					locali = partindex[i];
					localj = partindex[j];
					
					sigma=0.5*(mysigma[locali]+mysigma[localj]);
					epsilon=sqrt(myepsilon[locali]*myepsilon[localj]);
					rc=2.5*sigma;   // rc >= 2.5sigma  cutoff distance	
					
					
					xforce = x[localj];
					yforce = y[localj];
					zforce = z[localj];
					
					if (x[localj]-x[locali]>0.5*Lx) {xforce=x[localj]-Lx;}
					if (x[localj]-x[locali]<-0.5*Lx) {xforce=x[localj]+Lx;}
					if (y[localj]-y[locali]>0.5*Ly) {yforce=y[localj]-Ly;}
					if (y[localj]-y[locali]<-0.5*Ly) {yforce=y[localj]+Ly;}
					
					if (z[localj]-z[locali]>0.5*Lz) {zforce=z[localj]-Lz;}
					if (z[localj]-z[locali]<-0.5*Lz) {zforce=z[localj]+Lz;}	
					
					
						r2 = (xforce - x[locali])*(xforce - x[locali]) + (yforce - y[locali])*(yforce - y[locali]) + (zforce - z[locali])*(zforce - z[locali]);
						r=sqrt(r2);
									
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

			
							sumfx[locali]=sumfx[locali]+coef2*(x[locali]-xforce);
							sumfy[locali]=sumfy[locali]+coef2*(y[locali]-yforce);
							sumfz[locali]=sumfz[locali]+coef2*(z[locali]-zforce);	
		
							
							sumfx[localj]=sumfx[localj]-coef2*(x[locali]-xforce);
							sumfy[localj]=sumfy[localj]-coef2*(y[locali]-yforce);
							sumfz[localj]=sumfz[localj]-coef2*(z[locali]-zforce);									
							
							
						sumVcells=sumVcells+sigma12/mysqrt12-sigma6/mysqrt6;
						
					}
					
//--------------------------------------------start the rdfsum calculation				
					if (n>= N-rdf_avg && status[locali]==1 && status[localj]==1)  //for rdf of many time steps, n>=N-rdf_avg
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
					
					
					
				}//end of if i#j
				
			}}
	} // end of counter
	
		
//Part 3: Read the next cell (neighbour) and define the index of the particles that it contains, do not print empty cells
	
	for(next=0;next<26;++next){
		

		nextcounter = 0;
		HEADp = HEAD[neighbours[next][cell]];
		if(HEADp != -1){
			
			nextpartindex[nextcounter] = HEADp;
			nextcounter = nextcounter + 1;  
			LISTp = LIST[HEADp];
			while(LISTp != -1){
				if(LISTp != -1){
					
					nextpartindex[nextcounter] = LISTp;
					nextcounter = nextcounter + 1;   
					LISTp = LIST[LISTp];
				} //end if
			}//end while	
			
		}
		
		
//	Part 4 of interactions calculations
		
		if(nextcounter != 0){ //The loop must be calculated only if there are atoms in the cell we choose
			
			for(i=0;i<counter;++i){	
				for(j=0;j<nextcounter;++j){
					
					locali = partindex[i];
					localj = nextpartindex[j];
					sigma=0.5*(mysigma[locali]+mysigma[localj]);
					epsilon=sqrt(myepsilon[locali]*myepsilon[localj]);
					rc=2.5*sigma;   // rc >= 2.5sigma  cutoff distance
					
					
					xforce = x[localj];
					yforce = y[localj];
					zforce = z[localj];
					
					if (x[localj]-x[locali]>0.5*Lx) {xforce=x[localj]-Lx;}
					if (x[localj]-x[locali]<-0.5*Lx) {xforce=x[localj]+Lx;}
					if (y[localj]-y[locali]>0.5*Ly) {yforce=y[localj]-Ly;}
					if (y[localj]-y[locali]<-0.5*Ly) {yforce=y[localj]+Ly;}
					
					if (z[localj]-z[locali]>0.5*Lz) {zforce=z[localj]-Lz;}
					if (z[localj]-z[locali]<-0.5*Lz) {zforce=z[localj]+Lz;}	
					
					
					r2 = ((xforce - x[locali])*(xforce - x[locali]) + (yforce - y[locali])*(yforce - y[locali]) + (zforce - z[locali])*(zforce - z[locali]));
					r=sqrt(r2);
					
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
						
						sumfx[locali]=sumfx[locali]+coef2*(x[locali]-xforce);
						sumfy[locali]=sumfy[locali]+coef2*(y[locali]-yforce);
						sumfz[locali]=sumfz[locali]+coef2*(z[locali]-zforce);	
						
						
						sumfx[localj]=sumfx[localj]-coef2*(x[locali]-xforce);
						sumfy[localj]=sumfy[localj]-coef2*(y[locali]-yforce);
						sumfz[localj]=sumfz[localj]-coef2*(z[locali]-zforce);	
						
						
						sumVcells=sumVcells+sigma12/mysqrt12-sigma6/mysqrt6;
					}
					
//--------------------------------------------start the rdfsum calculation				
					if (n>= N-rdf_avg && status[locali]==1 && status[localj]==1)  //for rdf of many time steps, n>=N-rdf_avg
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
					
	
					
				}}
		} // end of nextcounter
		
	
} // end of next
}//end of cell
	
	
} // end of function	