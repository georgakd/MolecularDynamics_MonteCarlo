/*
 *  lj_ts_forces_cells.c
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


void lj_ts_forces_cells(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);

void lj_ts_forces_cells(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB)
{
	double coef1,coef2,sigma6,sigma8,sigma12,sigma14,mysqrt6,mysqrt8,mysqrt12,mysqrt14;
	
	int locali,localj;

	int counterA,counterB,nextcounter;
	int icellA,icellB;
	int HEADp,LISTp;

//----------------------------------------- Build HEAD and LIST for both populations-------------------------------------//
	
		
	for(cell=0;cell<cellsx_supercell*cellsy_supercell*cellsz_supercell;++cell)
	{
		HEAD_A[cell] = -1; //because 0 is the index of the first particle!!!!!!!
		HEAD_B[cell] = -1;	
	}

	for(l=0;l<npartA;++l)
	{
		LIST_A[l] = -1; //because 0 is the index of the first particle!!!!!!!		
	}
	for(l=0;l<npartB;++l)
	{
		LIST_B[l] = -1; //because 0 is the index of the first particle!!!!!!!
	}
	

	for(l=0;l<npartA;++l)
	{
		icellA = 
		(int)(floor((xA[l] - xmin)/Lcx_supercell) + floor((yA[l] - ymin)/Lcy_supercell)*cellsx_supercell + floor((zA[l] - zmin)/Lcz_supercell)*cellsx_supercell*cellsy_supercell);
		
		LIST_A[l] = HEAD_A[icellA];
		HEAD_A[icellA] = l; 

	}


		for(l=0;l<npartB;++l)
	{
		icellB = 
		(int)(floor((xB[l] - xmin)/Lcx_supercell) + floor((yB[l] - ymin)/Lcy_supercell)*cellsx_supercell + floor((zB[l] - zmin)/Lcz_supercell)*cellsx_supercell*cellsy_supercell);
				
		LIST_B[l] = HEAD_B[icellB];
		HEAD_B[icellB] = l; 
		
	}

//-------------------------------------------------------------------------------------------------------------------------------------------//

	for (k=0;k<npartA;++k){sumfxA[k]=0;sumfyA[k]=0;sumfzA[k]=0;}
	for (k=0;k<npartB;++k){sumfxB[k]=0;sumfyB[k]=0;sumfzB[k]=0;}
	
	sumV_ts_cells=0;
	

	// 1. Read each current cell and define particles of population A in each cell
	
for (cell=0;cell<cellsx_supercell*cellsy_supercell*cellsz_supercell;++cell){		  
	
	counterA = 0;
	HEADp = HEAD_A[cell];
	if(HEADp != -1){				
		partindex_A[counterA] = HEADp;
		counterA = counterA + 1;		
		LISTp = LIST_A[HEADp];

		while(LISTp != -1){
			if(LISTp != -1){
				partindex_A[counterA] = LISTp;
				counterA = counterA + 1;				
				LISTp = LIST_A[LISTp];				
			} //end of if LISTp
		} //end of while 

	} // end of if HEADp



	counterB = 0;
	HEADp = HEAD_B[cell];
	if(HEADp != -1){
				
		partindex_B[counterB] = HEADp;
		
		counterB = counterB + 1;
		LISTp = LIST_B[HEADp];
				
		while(LISTp != -1){
			if(LISTp != -1){
				
				partindex_B[counterB] = LISTp;

				counterB = counterB + 1;
	
				LISTp = LIST_B[LISTp];
				
			} //end of if LISTp
		} //end of while 
	} // end of if HEADp

	//printf("%lf",counterB);
	//getchar();



	//2. Calculate the interactions with population B in the current cell
if (counterA!=0 && counterB!=0) // be  careful, this if is very important
{
	for (i=0;i<counterA;++i)
	{
		for(j=0;j<counterB;++j)
		{
			locali = partindex_A[i];
			localj = partindex_B[j];
			sigma_ts=0.5*(mysigmaA[locali]+mysigmaB[localj]);
			epsilon_ts=sqrt(myepsilonA[locali]*myepsilonB[localj]);
			rc_ts=2.5*sigma_ts;   // rc >= 2.5sigma  cutoff distance
			
			xforce=xB[localj];yforce=yB[localj];zforce=zB[localj];

			if (xB[localj]-xA[locali]>0.5*Lx) {xforce=xB[localj]-Lx;}
			if (xB[localj]-xA[locali]<-0.5*Lx) {xforce=xB[localj]+Lx;}
			if (yB[localj]-yA[locali]>0.5*Ly) {yforce=yB[localj]-Ly;}
			if (yB[localj]-yA[locali]<-0.5*Ly) {yforce=yB[localj]+Ly;}
		
			if (zB[localj]-zA[locali]>0.5*Lz) {zforce=zB[localj]-Lz;}
			if (zB[localj]-zA[locali]<-0.5*Lz) {zforce=zB[localj]+Lz;}
	
			r2=(xforce-xA[locali])*(xforce-xA[locali])+(yforce-yA[locali])*(yforce-yA[locali])+(zforce-zA[locali])*(zforce-zA[locali]);		
			
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


				sumfxA[locali]=sumfxA[locali]+coef2*(xA[locali]-xforce);
				sumfyA[locali]=sumfyA[locali]+coef2*(yA[locali]-yforce);
				sumfzA[locali]=sumfzA[locali]+coef2*(zA[locali]-zforce);					
				
				sumfxB[localj]=sumfxB[localj]-coef2*(xA[locali]-xforce);
				sumfyB[localj]=sumfyB[localj]-coef2*(yA[locali]-yforce);
				sumfzB[localj]=sumfzB[localj]-coef2*(zA[locali]-zforce);	
				
				
			sumV_ts_cells = sumV_ts_cells +sigma12/mysqrt12-sigma6/mysqrt6;
			
			
			}
				
		}
	}

} //----------------------- end of part 2-------------------------------

// 3. Read neighbour cells and define particles of population B

for(next=0;next<26;++next){
		
		nextcounter = 0;
		HEADp = HEAD_B[neighbours_supercell[next][cell]];
		if(HEADp != -1){
			
			nextpartindex_B[nextcounter] = HEADp;
			nextcounter = nextcounter + 1;  
			LISTp = LIST_B[HEADp];
			while(LISTp != -1){
				if(LISTp != -1){
					
					nextpartindex_B[nextcounter] = LISTp;
					nextcounter = nextcounter + 1;   
					LISTp = LIST_B[LISTp];
				} //end if
			}//end while	
			
		}

//4. Define interactions of A with B in the neighbour cells

if (nextcounter!=0)
{
	for (i=0;i<counterA;++i)
	{
		for(j=0;j<nextcounter;++j)
		{
			
			locali = partindex_A[i];
			localj = nextpartindex_B[j];
			sigma_ts=0.5*(mysigmaA[locali]+mysigmaB[localj]);
			epsilon_ts=sqrt(myepsilonA[locali]*myepsilonB[localj]);
			rc_ts=2.5*sigma_ts;   // rc >= 2.5sigma  cutoff distance
			
			xforce=xB[localj];yforce=yB[localj];zforce=zB[localj];

			if (xB[localj]-xA[locali]>0.5*Lx) {xforce=xB[localj]-Lx;}
			if (xB[localj]-xA[locali]<-0.5*Lx) {xforce=xB[localj]+Lx;}
			if (yB[localj]-yA[locali]>0.5*Ly) {yforce=yB[localj]-Ly;}
			if (yB[localj]-yA[locali]<-0.5*Ly) {yforce=yB[localj]+Ly;}
		
			if (zB[localj]-zA[locali]>0.5*Lz) {zforce=zB[localj]-Lz;}
			if (zB[localj]-zA[locali]<-0.5*Lz) {zforce=zB[localj]+Lz;}
	
			r2=(xforce-xA[locali])*(xforce-xA[locali])+(yforce-yA[locali])*(yforce-yA[locali])+(zforce-zA[locali])*(zforce-zA[locali]);		
			
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

				sumfxA[locali]=sumfxA[locali]+coef2*(xA[locali]-xforce);
				sumfyA[locali]=sumfyA[locali]+coef2*(yA[locali]-yforce);
				sumfzA[locali]=sumfzA[locali]+coef2*(zA[locali]-zforce);					
				
				sumfxB[localj]=sumfxB[localj]-coef2*(xA[locali]-xforce);
				sumfyB[localj]=sumfyB[localj]-coef2*(yA[locali]-yforce);
				sumfzB[localj]=sumfzB[localj]-coef2*(zA[locali]-zforce);		
							
			sumV_ts_cells = sumV_ts_cells +sigma12/mysqrt12-sigma6/mysqrt6;
					
			}
			
		}
	}

} //----------------------- end of part 4-------------------------------



} // close the reading of neighbour cells


} // close the reading of all cells



}