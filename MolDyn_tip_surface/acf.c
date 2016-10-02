/*
 *  acf.c
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

void acf(double *vx,double *vy,double *vz);

void acf(double *vx,double *vy,double *vz)
{
	double sumparticles;
//________________________________start the ACF calculation for Total DOS----------------------------------------------------------------------------------
	
 if (n == startpoint)
 {	
 for (k=0;k<npartB;++k)
 {
 if (status[k]==1)
 {
 vxfixed[k]=vx[k];
 vyfixed[k]=vy[k];
 vzfixed[k]=vz[k];
 acf_t[k]=0;
 }
 }
 sumtime[0]=1;
 
 }
 
 if (n > startpoint)
 {	
 sumparticles=0;
 
 for (k=0;k<npartB;++k)
 {
 if (status[k]==1)
 {
 acf_t[k]=acf_t[k]+(vx[k]*vxfixed[k]+vy[k]*vyfixed[k]+vz[k]*vzfixed[k]);
 sumparticles=sumparticles+acf_t[k]/(vxfixed[k]*vxfixed[k]+vyfixed[k]*vyfixed[k]+vzfixed[k]*vzfixed[k]);
 }
 }
 
 sumtime[n-startpoint] = sumparticles/(mobile*(n-startpoint)); // this is correct, do not change it!!!!!!
 }
 
} // end of function
