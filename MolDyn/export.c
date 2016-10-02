/*
 *  export.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 04/06/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

FILE *fp;

void export_info(void);

void export_info(void){

//________________________export last positions and velocities according to structure chosen_______________


	sprintf(path,"%sresume.dat",pathout);
	fp=fopen(path,"w+");
	fprintf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",title,npart,xmin,xmax,ymin,ymax,zmin,zmax);
	fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");	

	for(k=0;k<npart;++k)
	{
		fprintf(fp,"%s\t%.16lf\t%.16lf\t%.16lf\t%d\t%d\t%.16lf\t%.16lf\t%.16lf\n",species[k],x2[k],y2[k],z2[k],status[k],thermostatus[k],vx2[k],vy2[k],vz2[k]);
	}
	fclose(fp);
	
	
//___________write to file energies.dat_______________________________________________	
	
	sprintf(path,"%sthermodynamic_quantities_%s%d.dat",pathout,ensemble,(int)*pointer);
	fp=fopen(path,"w+");

	for(n=0;n<=N;++n)
	{
		fprintf(fp,"%d\t"
				"%.16lf\t"				
				"%.16lf\t"
				"%.16lf\t"
				"%.16lf\n",
				n,
				K[n]/npart,
				U[n]/npart,
				E[n]/npart,
				T[n]
				);
	}
	fclose(fp);
	
	
//___________write to file acf.dat_______________________________________________
	sprintf(path,"%sacf_%s%d.dat",pathout,ensemble,(int)*pointer);
	fp=fopen(path,"w+");
	
	for (n=0;n<=N-startpoint;++n)
	{
		fprintf(fp,"%d\t%.16lf\n",n,sumtime[n]);
	}
	fclose(fp);

	
	
//___________write to file rdf.dat_______________________________________________
	sprintf(path,"%srdf_%s%d.dat",pathout,ensemble,(int)*pointer);
	fp=fopen(path,"w+");
	
	for (p=0;p<bins;++p)
	{
		fprintf(fp,"%.16lf\t%.16lf\n",r1[p],g[p]);
	}
	fclose(fp);


}