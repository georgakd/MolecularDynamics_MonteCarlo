/*
 *  parameters.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 22/03/2012.
 *  Copyright 2012 Home. All rights reserved.
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
		exit(1);}
		
		fscanf(fp_config,"mypotential: %d\n",&mypotential);
		fscanf(fp_config,"dtSI: %lf\n",&dtSI);
		fscanf(fp_config,"N: %d\n",&N);
		fscanf(fp_config,"nstep: %d\n",&nstep);
		fscanf(fp_config,"graph_step: %d\n",&graph_step);
		fscanf(fp_config,"cellsx: %d\n",&cellsx);
		fscanf(fp_config,"cellsy: %d\n",&cellsy);
		fscanf(fp_config,"cellsz: %d\n",&cellsz);
		fscanf(fp_config,"linkcells: %d\n",&linkcells);
		fscanf(fp_config,"resume_flag: %d\n",&resume_flag);
		fscanf(fp_config,"csvout: %d\n",&csvout);
		fscanf(fp_config,"jmolout: %d\n",&jmolout);
		fscanf(fp_config,"monitor: %d\n",&monitor);
		fscanf(fp_config,"aneal: %d\n",&aneal);
		fscanf(fp_config,"choose_thermostat: %d\n",&choose_thermostat);	
		fscanf(fp_config,"Tintervals: %d\n",&Tintervals);
		fscanf(fp_config,"nbath: %d\n",&nbath);
		fscanf(fp_config,"nbath_end: %d\n",&nbath_end);
		fscanf(fp_config,"Tinitial: %lf\n",&Tinitial);
		fscanf(fp_config,"Ttarget: %lf\n",&Ttarget);
		fscanf(fp_config,"Tcool: %lf\n",&Tcool);
		fscanf(fp_config,"nu: %lf\n",&nu);
		fscanf(fp_config,"trise: %lf\n",&trise);
		fscanf(fp_config,"startpoint: %d\n",&startpoint);
		fscanf(fp_config,"rdf_avg: %d\n",&rdf_avg);
		fscanf(fp_config,"dr: %lf\n",&dr);
		fscanf(fp_config,"range: %lf\n",&range);
		
		fclose(fp_config);

//****************************************Fixed Parameters******************************************************//

	w = 3456; // for initial seeding for the random number generator
	x = 1357; // for initial seeding for the random number generator
	y = 9876; // for initial seeding for the random number generator
	z = 4792; // for initial seeding for the random number generator
	dtscale = 10.1805; // see the file in mathematica
	kbol = 8.6173324e-5;  // Boltzmann constant
	pi = 3.14159265;
	



/*
**********************Comments for parameters**************************************************************

	mypotential = 0; // 0 for tersoff (dtSI~0.5),  1 for lj (dtSI~10)

	aneal = 1; // 1 = sim anneal, -1 = sim quench
	choose_thermostat = 2; // 0 = NVE, 1 = ANDERSEN, 2 = Berendsen
	Tintervals = 1; // be very careful, nbath_end/Tintervals must have modulus zero!!!!!!!!
	nbath = 1;   // when to start the thermostat
	nbath_end = N; // when to end the thermostat	
	Tinitial = 0; // start Temperature in sim aneal if we are not supposed to start from T=0K the thermostated intervals
	Ttarget = 300;  // Temp of thermostat in Kelvin that the system should relax to
	Tcool = 10; //  Temp of thermostat in Kelvin that the system should cool to	
	nu = 0.5; //collision parameter with heat bath for Andersen
	trise = 10;  // rise time of Berendsen thermostat ~ (in fs)	
	//nstep = 500; //for printf in console
*************************************************************************************************************
*/	
}
