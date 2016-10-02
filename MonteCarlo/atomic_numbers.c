/*
 *  atomic numbers.c
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
#include <string.h>

void atomic_numbers(void);
void atomic_numbers(void)
{
	for (k=0;k<npart;++k)
	{
		
		if (strcmp(species[k],"Ar")==0)
		{
			// Argon
			Z[k]=18;   // atomic number
			m[k]=39.948; // atomic mass in amu
			covalent_r[k]=1.06; // Covalent Radius in A
			vdW_r[k]=1.88; // vand der Waals Radius in A
		}	
//---------------------------------------------------------------------------//
		if (strcmp(species[k],"C")==0)
		{
			// Carbon
			Z[k]=6;
			m[k]=12.0107;
			covalent_r[k]=0.77; // Covalent Radius 
			vdW_r[k]=1.70; // vand der Waals Radius 
		}
		
		
		
		if (strcmp(species[k],"Si")==0)
		{
			// Silicon
			Z[k]=14;
			m[k]=28.08553;
			covalent_r[k]=1.11; // Covalent Radius
			vdW_r[k]=2.1;   // vand der Waals Radius
		}
		
		
		if (strcmp(species[k],"Ge")==0)
		{
			// Germanium
			Z[k]=32;
			m[k]=72.59;
			covalent_r[k]=1.22; // Covalent Radius
			vdW_r[k]=2.11; // vand der Waals Radius
		}
		
//---------------------------------------------------------------------------//
		
		if (strcmp(species[k],"Cu")==0)
		{
			// Copper
			Z[k]=29;
			m[k]=63.546;
			covalent_r[k]=1.17; // Covalent Radius
			vdW_r[k]=1.40; // vand der Waals Radius
		}
		
		
		
	}
}
