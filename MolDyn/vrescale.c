/*
 *  vrescale.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 22/03/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void vrescale(double *vx,double *vy,double *vz);

void vrescale(double *vx,double *vy,double *vz)
{
	if (n>=nbath && n<=nbath_end)
	{
		
		sumkin=0;
		for (k=0;k<npart;++k)
		{
			if (thermostatus[k]==1){
				
				sumkin=sumkin+m[k]*vx[k]*vx[k]+m[k]*vy[k]*vy[k]+m[k]*vz[k]*vz[k];
			}
			
		}
		
		K[n]=0.5*sumkin;
		
		T[n]=(2*K[n])/(3*(thermo-1)*kbol); // if we use the configuration fixed-thermo-mobile
		
		// can be imaginary-be careful with the choice of trise!!!!
		constant=sqrt(1+(dtSI/trise)*(thermostat[n]/T[n]-1));
		
		for (k=0;k<npart;++k)
		{
			
			if (thermostatus[k]==1){ // if we use the configuration fixed-thermo-mobile
				
				vx[k]=constant*vx[k];  // Berendsen thermostat
				vy[k]=constant*vy[k];  
				vz[k]=constant*vz[k];  
				
			}
		}
	}
	
	
}