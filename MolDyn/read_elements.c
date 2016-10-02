/*
 *  read_elements.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 27/02/2012.
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
	
	dx_zeta_III=(double*)malloc(npart*sizeof(double));				//<---type III interactions
	dy_zeta_III=(double*)malloc(npart*sizeof(double));
	dz_zeta_III=(double*)malloc(npart*sizeof(double));
	
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
		
//---------------------------------------------------------------------------//
	
		
		
		if (Z[k]==6)
		{
			
			//Carbon (parameters after Tersoff erratum)
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
			

			/*
			// Carbon (parameters Porter)
			myA[k] = 1.5448e+3; //eV
			myB[k] = 389.63; //eV //corrected in erratum
			lamda[k] = 3.4653; // 1/Angstrom
			mi[k] = 2.3064; // 1/Angstrom
			beta[k] = 4.1612e-6;
			ni[k] = 0.99054;
			ci[k] = 19981;
			di[k] = 7.0340; // corrected in erratum
			hi[k] = -0.39953;
			myR[k] =  1.8; //Angstrom
			myS[k] = 2.1; //Angstrom
			*/


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



		
//......more elements.......//
		
		
		
		
	}
}