/*
 *  jmol_func.c
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

void jmol_func(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);

FILE *fp; //file pointer to read and write data

//_________________jmol________________________________________________________

void jmol_func(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB)
{
	char xyzexport[100]; // for jmol format
	sprintf(title,"tip_sample"); // title for jmol format
	sprintf(xyzexport,"%sxyzout%d.xyz",pathxyz,n);
	fp=fopen(xyzexport,"w+");
	fprintf(fp,"%d\n%s\n",(npartA+npartB),title);
	
	for (k=0;k<npartB;++k) //surface
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species_B[k],xB[k],yB[k],zB[k]);
	}
	
	for (k=0;k<npartA;++k) //tip
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species_A[k],xA[k],yA[k],zA[k]);
	}
	
	fclose(fp);

}