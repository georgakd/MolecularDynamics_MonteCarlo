/*
 *  read_elements.c
 *  tersoff calculations
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
	
	//Tersoff
	myA=(double*)malloc(npart*sizeof(double));
	myB=(double*)malloc(npart*sizeof(double));
	lamda=(double*)malloc(npart*sizeof(double));
	mi=(double*)malloc(npart*sizeof(double));
	beta=(double*)malloc(npart*sizeof(double));
	ni=(double*)malloc(npart*sizeof(double));
	ci=(double*)malloc(npart*sizeof(double));
	di=(double*)malloc(npart*sizeof(double));
	hi=(double*)malloc(npart*sizeof(double));
	myR=(double*)malloc(npart*sizeof(double));
	myS=(double*)malloc(npart*sizeof(double));
	
	//LJ
	myepsilon=(double*)malloc(npart*sizeof(double));
	mysigma=(double*)malloc(npart*sizeof(double));
	
	
	for (k=0;k<npart;++k)
	{	
		
		if (Z[k]==0){printf("No element selected........Error!!!!!!!!");}
		
		if (Z[k]==18)
		{
			//Argon
			myepsilon[k] = 0.0102985; // eV
			mysigma[k] = 3.4; // Angstrom
		}

		
		if (Z[k]==6)
		{
			// Carbon
			myA[k] = 1.3936e+3; //eV
			myB[k] = 3.4674e+2; //eV
			lamda[k] = 3.4879; // 1/Angstrom
			mi[k] = 2.2119; // 1/Angstrom
			beta[k] = 1.5724e-7;
			ni[k] = 7.2751e-1;
			ci[k] = 3.8049e+4;
			di[k] = 4.3484;
			hi[k] = -5.7058e-1;
			myR[k] =  1.8; //Angstrom
			myS[k] = 2.1; //Angstrom
		}
		
		
		if (Z[k]==14)
		{
			// Silicon
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

			if (Z[k]==32)
		{
			// Germanium
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


		if (Z[k]==29)
		{
			//Copper
			myepsilon[k] = 0.012382; // eV
			mysigma[k] = 3.61; // Angstrom
			mpower=6;
			npower=9;
			c_param=39.432;
						
		}
		
		
		
		
	}
}




