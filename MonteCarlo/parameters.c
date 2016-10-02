/*
 *  parameters.c
 *  tersoff calculations
 *
 *  Created by Dimitra Georgakaki on 02/09/2011.
 *  Copyright 2011 Home. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

FILE *fp_config;


void parameters(void);

void parameters(void)
{
	
	fp_config=fopen(file_path,"r");
		if(fp_config==NULL){printf("Could not locate %s",file_path);printf("\nAborting Simulation...Press 1 for exit...\n");
		system("PAUSE");
		exit(1);}
		
		fscanf(fp_config,"mypotential: %d\n",&mypotential);
		fscanf(fp_config,"max_cycles: %d\n",&max_cycles);
		fscanf(fp_config,"cellsx: %d\n",&cellsx);
		fscanf(fp_config,"cellsy: %d\n",&cellsy);
		fscanf(fp_config,"cellsz: %d\n",&cellsz);
		fscanf(fp_config,"linkcells: %d\n",&linkcells);
		fscanf(fp_config,"resume_flag: %d\n",&resume_flag);
		fscanf(fp_config,"csvout: %d\n",&csvout);
		fscanf(fp_config,"jmolout: %d\n",&jmolout);
		fscanf(fp_config,"nsample: %d\n",&nsample);
		fscanf(fp_config,"dmax: %lf\n",&dmax);
		fscanf(fp_config,"Tref: %lf\n",&Tref);
		
		fclose(fp_config);


//------------------------fixed parameters and comments--------------------------------------//

//	w = 10745; x = 10923; y = 10578; z = 10334; // initial values for seeding the random number generator for random disp	
// Tref = 300; // simple temperature ref	
// dmax = 0.09; // choose wisely according to case study!!!!!!
	
// the same seed as the other projects
	w = 3456; // for initial seeding for the random number generator
	x = 1357; // for initial seeding for the random number generator
	y = 9876; // for initial seeding for the random number generator
	z = 4792; // for initial seeding for the random number generator

	pi = 3.14159265;
	kbol=8.6173324e-5;  // Boltzmann constant


	
}