/*
 *  jmol_func.c
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

void jmol_func(double *x,double *y,double *z);

FILE *fp; //file pointer to read and write data

//_________________jmol____________________________________

void jmol_func(double *x,double *y,double *z)
{
	
	sprintf(myexport,"%sxyzout%d.xyz",pathxyz,n);
	fp=fopen(myexport,"w+");
	fprintf(fp,"%d\n%s\n",npart,title);
	
	for (k=0;k<npart;++k)
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species[k],x[k],y[k],z[k]);
	}
	fclose(fp);


}