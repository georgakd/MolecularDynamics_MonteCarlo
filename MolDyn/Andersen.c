/*
 *  Andersen.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 13/03/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double random_gen(void); // include this function to my code

void Andersen(double *vx,double *vy,double *vz);

void Andersen(double *vx,double *vy,double *vz)
{
		
	 uxcm=0;uycm=0;uzcm=0;
	 
	 if (n>=nbath && n<=nbath_end)
	 {
	 number = random_gen();
	 
	 if (number<nu)
	 {
	 for (k=0;k<npart;++k){
	 if (thermostatus[k]==1){   // if we use the configuration fixed-thermo-mobile
		 
	 myrandom1 = random_gen();
	 myrandom2 = random_gen();
	 vx[k]=sqrt((kbol*Ttarget)/m[k])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2); // Box-Muller Method
	 vy[k]=sqrt((kbol*Ttarget)/m[k])*sqrt(-2*log(myrandom1))*sin(2*pi*myrandom2);
	 
	 myrandom1 = random_gen();
	 myrandom2 = random_gen();
	 vz[k]=sqrt((kbol*Ttarget)/m[k])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2);
	

	 uxcm=uxcm+vx[k];uycm=uycm+vy[k];uzcm=uzcm+vz[k];  //velocity center of mass for mobile particles
	 }
	 }
	 
	 uxcm=uxcm/thermo;uycm=uycm/thermo;uzcm=uzcm/thermo;    // if we use the configuration fixed-thermo-mobile
	 
	 for (k=0;k<npart;++k)
	 {	
	 
	 if (thermostatus[k]==1){    // if we use the configuration fixed-thermo-mobile 
	 
	 vx[k]=vx[k]-uxcm;
	 vy[k]=vy[k]-uycm;
	 vz[k]=vz[k]-uzcm;

	
		
	 }
	 }
	 
	 }
	 
	 
	 
	 
	 }
	 

	
}