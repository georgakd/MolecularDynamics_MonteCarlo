/*
 *  read_elements.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 18/01/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */
#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void read_elements(void);
void read_elements(void)
{
	//Tersoff (example Si/C/Ge sample)
	myA=(double*)malloc(npartB*sizeof(double));
	myB=(double*)malloc(npartB*sizeof(double));
	lamda=(double*)malloc(npartB*sizeof(double));
	mi=(double*)malloc(npartB*sizeof(double));
	beta=(double*)malloc(npartB*sizeof(double));
	ni=(double*)malloc(npartB*sizeof(double));
	ci=(double*)malloc(npartB*sizeof(double));
	di=(double*)malloc(npartB*sizeof(double));
	hi=(double*)malloc(npartB*sizeof(double));
	myR=(double*)malloc(npartB*sizeof(double));
	myS=(double*)malloc(npartB*sizeof(double));
	
	dx_zeta_III = (double*)malloc(npartB*sizeof(double));
	dy_zeta_III = (double*)malloc(npartB*sizeof(double));
	dz_zeta_III = (double*)malloc(npartB*sizeof(double));
	
	
	// LJ interactions between 2 populations
	myepsilonA=(double*)malloc(npartA*sizeof(double));
	mysigmaA=(double*)malloc(npartA*sizeof(double));	
	myepsilonB=(double*)malloc(npartB*sizeof(double));
	mysigmaB=(double*)malloc(npartB*sizeof(double));
		
	
	//LJ (example Argon)
	myepsilon=(double*)malloc(npartB*sizeof(double));
	mysigma=(double*)malloc(npartB*sizeof(double));
	
//-----------------------------------Population A-----------------------------------------//
	

	
	for (k=0;k<npartA;++k)
	{	
		
		if (Z_A[k]==0){printf("No element selected........Error!!!!!!!!");}
		
		
		
		if (Z_A[k]==14) // e.g. Silicon tip
		{	
			//LJ
			myepsilonA[k]=0.01182;
			mysigmaA[k]=3.7;			
		}
		
		if (Z_A[k]==6)
		{	
			//LJ
			myepsilonA[k]=0.00284;
			mysigmaA[k]=3.4;	
		}
				
		
	}
	
//-----------------------Population B----------------------------------------------------//
	
	for (k=0;k<npartB;++k)
	{	
		
		if (Z_B[k]==0){printf("No element selected........Error!!!!!!!!");}
		
		
		if (Z_B[k]==18) //Argon
		{		
			myepsilon[k] = 0.0102985; // eV
			mysigma[k] = 3.4; // Angstrom		
		}
		
		if (Z_B[k]==14) //Silicon
		{
			//LJ
			myepsilonB[k]=0.01182;
			mysigmaB[k]=3.7;
			
			// Tersoff
			myA[k]=1830.8;
			myB[k]=471.18;
			lamda[k]=2.4799;
			mi[k]=1.7322;
			beta[k]=1.1e-6;
			ni[k]=0.78734;
			ci[k]=100390.0;
			di[k]=16.217;
			hi[k]=-0.59825;
			myR[k]=2.7;
			myS[k]=3;
		}
		
		if (Z_B[k]==6) //Carbon
		{
			
			//LJ
			myepsilonB[k] = 0.00284;
			mysigmaB[k] = 3.4;			
			
			//(parameters after Tersoff erratum)
			myA[k] = 1393.6; //eV
			myB[k] = 346.74; //eV //corrected in erratum
			lamda[k] = 3.4879; // 1/Angstrom
			mi[k] = 2.2119; // 1/Angstrom
			beta[k] = 1.5724e-7;
			ni[k] = 0.72751;
			ci[k] = 38049;
			di[k] = 4.3484; // corrected in erratum
			hi[k] = -0.57058;
			myR[k] =  1.8; //Angstrom
			myS[k] = 2.1; //Angstrom
			
		}
		
		
		if (Z_B[k]==32) // Germanium
		{
			
			myA[k]=1769;
			myB[k]=419.23;
			lamda[k]=2.4451;
			mi[k]=1.7047;
			beta[k]=9.0166e-7;
			ni[k]=0.75627;
			ci[k]=1.0643e5;
			di[k]=15.652;
			hi[k]=-0.43884;
			myR[k]=2.8;
			myS[k]=3.1;
		}
		

		
		
		
	}
	
	
	
	
	
	
	
	
	
	
	
}




