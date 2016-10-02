/*
 *  paraview_func.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 20/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void paraview_func(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);

FILE *fp; //file pointer to read and write data


//_________________paraview_____________________________________________________________________

void paraview_func(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB)
{
	int csv_color;
	char myexport[100];
	
	sprintf(myexport,"%sout%d.csv",pathcsv,n);
	fp=fopen(myexport,"w+");
	fprintf(fp,"x,y,z,r,color\n");
	//surface
	for(k=0;k<npartB;++k)
	{
		csv_color=2;
		fprintf(fp,"%.16lf,%.16lf,%.16lf,%lf,%d\n",xB[k],yB[k],zB[k],covalent_r_B[k],csv_color);
	}
	//tip
	for(k=0;k<npartA;++k)
	{
		csv_color=1;
		fprintf(fp,"%.16lf,%.16lf,%.16lf,%lf,%d\n",xA[k],yA[k],zA[k],covalent_r_A[k],csv_color);
	}
	
	fclose(fp);
	
}
