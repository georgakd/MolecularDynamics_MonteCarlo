/*
 *  paraview_func.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 30/01/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void paraview_func(double *x,double *y,double *z);

FILE *fp; //file pointer to read and write data


//_________________paraview_____________________________________________________________________

void paraview_func(double *x,double *y,double *z)
{
	sprintf(myexport,"%sout%d.csv",pathcsv,n);
	fp=fopen(myexport,"w+");
	fprintf(fp,"x,y,z,r,color\n");
	for(k=0;k<npart;++k)
	{
		if (status[k]==0){csv_color=0;}
		if (status[k]==1 && thermostatus[k]==1){csv_color=1;}
		if (status[k]==1 && thermostatus[k]==0){csv_color=2;}
		fprintf(fp,"%.16lf,%.16lf,%.16lf,%lf,%d\n",x[k],y[k],z[k],covalent_r[k],csv_color);
	}
	fclose(fp);

}