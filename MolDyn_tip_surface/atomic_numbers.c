/*
 *  atomic_numbers.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 19/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */
#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void atomic_numbers(void);
void atomic_numbers(void)
{
	
	
//-------------------------------------CASE OF TWO DIFFERENT POPULATIONS---------------------------------//	
	
	for (k=0;k<npartA;++k)
	{
		
		if (strcmp(species_A[k],"Ar")==0)
		{
			// Argon A
			Z_A[k]=18;   // atomic number
			m_A[k]=39.948; // atomic mass in amu
			covalent_r_A[k]=1.06; // Covalent Radius in A
			vdW_r_A[k]=1.88; // vand der Waals Radius in A
		}
		
		if (strcmp(species_A[k],"C")==0)
		{
			// Carbon A
			Z_A[k]=6;
			m_A[k]=12.0107;
			covalent_r_A[k]=0.77; 
			vdW_r_A[k]=1.70; 
		}
		
		if (strcmp(species_A[k],"Si")==0)
		{
			// Silicon A
			Z_A[k]=14;
			m_A[k]=28.08553;
			covalent_r_A[k]=1.11; 
			vdW_r_A[k]=2.1;   
		}
		
	}
	

	for (k=0;k<npartB;++k)
	{
		if (strcmp(species_B[k],"Ar")==0)
		{
			// Argon B
			Z_B[k]=18;   // atomic number
			m_B[k]=39.948; // atomic mass in amu
			covalent_r_B[k]=1.06; // Covalent Radius in A
			vdW_r_B[k]=1.88; // vand der Waals Radius in A
		}
		
		if (strcmp(species_B[k],"C")==0)
		{
			// Carbon B
			Z_B[k]=6;
			m_B[k]=12.0107;
			covalent_r_B[k]=0.77;  
			vdW_r_B[k]=1.70; 
		}
				
		if (strcmp(species_B[k],"Si")==0)
		{
			// Silicon B
			Z_B[k]=14;
			m_B[k]=28.08553;
			covalent_r_B[k]=1.11; 
			vdW_r_B[k]=2.1;   
		}	
		
	}
	
}
