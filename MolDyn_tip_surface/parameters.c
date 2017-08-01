/*
 *  parameters.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 20/04/2012.
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
		fscanf(fp_config,"nbath_end: %d\n",&nbath_end);
		fscanf(fp_config,"nbath: %d\n",&nbath);
		fscanf(fp_config,"aneal: %d\n",&aneal);
		fscanf(fp_config,"choose_thermostat: %d\n",&choose_thermostat);
		fscanf(fp_config,"resume_surf: %d\n",&resume_surf);
		fscanf(fp_config,"resume_tip: %d\n",&resume_tip);
		fscanf(fp_config,"csvout: %d\n",&csvout);
		fscanf(fp_config,"jmolout: %d\n",&jmolout);
		fscanf(fp_config,"monitor: %d\n",&monitor);
		fscanf(fp_config,"nstep: %d\n",&nstep);
		fscanf(fp_config,"graph_step: %d\n",&graph_step);
		fscanf(fp_config,"tip_apex: %d\n",&tip_apex);
		fscanf(fp_config,"linkcells: %d\n",&linkcells);
		fscanf(fp_config,"linkcells_supercell: %d\n",&linkcells_supercell);
		fscanf(fp_config,"cellsx: %d\n",&cellsx);
		fscanf(fp_config,"cellsy: %d\n",&cellsy);
		fscanf(fp_config,"cellsz: %d\n",&cellsz);
		fscanf(fp_config,"cellsx_supercell: %d\n",&cellsx_supercell);
		fscanf(fp_config,"cellsy_supercell: %d\n",&cellsy_supercell);
		fscanf(fp_config,"cellsz_supercell: %d\n",&cellsz_supercell);
		fscanf(fp_config,"Tintervals: %d\n",&Tintervals);
		fscanf(fp_config,"Ttarget: %lf\n",&Ttarget);
		fscanf(fp_config,"Tcool: %lf\n",&Tcool);
		fscanf(fp_config,"Tinitial: %lf\n",&Tinitial);
		fscanf(fp_config,"nu: %lf\n",&nu);
		fscanf(fp_config,"trise: %lf\n",&trise);
		fscanf(fp_config,"startpoint: %d\n",&startpoint);
		fscanf(fp_config,"rdf_avg: %d\n",&rdf_avg);
		fscanf(fp_config,"dr: %lf\n",&dr);
		fscanf(fp_config,"range: %lf\n",&range);
		

	fclose(fp_config);

//---------------------------------------comments-------------------------------------------------------------------------------//
	
	// mypotential=0; // 0 for tersoff,  1 for lj
	// N :   number of steps for surface relaxation // N/intervals must have modulus zero
	//	nbath = 1;   // when to start the thermostat
	// nbath_end = N; // when to end the thermostat
	//aneal = 1; // 1 = sim anneal, -1 = sim quench
	//choose_thermostat = 2; // 0 = NVE, 1 = ANDERSEN, 2 = Berendsen
	//resume_surf = 0; // if 0 do not resume, if 1 resume
	//resume_tip = 0; // if 0 do not resume, if 1 resume		
	//csvout = 0;   // flag for paraview (0=do not export, 1=export)
	//monitor = 0; // monitor for matlab (0=do not export, 1=export)
	//	tip = 2; // cases for tip movement
	//graph_step  //each ith iteration a jmol file is generated
	//nstep=10; //each ith iteration the results are printed at the screen
	//tip_apex=0; //choose the label of the particle that represents the tip apex
	//linkcells = 1; //flag for the use of linked cells in the sample
	//linkcells_supercell = 1; //flag for the use of linked cells in the supercell
	//	Tintervals = 1; // be very careful, nbath_end/Tintervals must have modulus zero!!!!!!!!

	//The Initial Temp in Kelvin (step n=0) is given by Tinit[0] or Ttarget according to the selected sim type
	//Ttarget = 2;  // Temp of thermostat in Kelvin that the system should relax to
	//Tcool = 20; // Temp of thermostat in Kelvin that the system should cool to
	//Tinitial=0; // if 0 we start sim aneal from the beginning
	//nu=0.5; //collision parameter with heat bath for Andersen
	//trise = 10;  // rise time of Berendsen thermostat ~ (in fs)	



	dtscale = 10.1805; // see the file in mathematica, in fs
	Ftargetscale = 1.60217; //in nN
	kscale = 16.0217; //in N/m
	velocityscale = 9822.67; //in m/s

	pi = 3.14159265;
	kbol = 8.6173324e-5;  // Boltzmann constant
	w=3456;x=1357;y=9876;z=4792; // for initial seeding for the random number generator



	
}
