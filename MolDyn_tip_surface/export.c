/*
 *  export.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 31/10/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

FILE *fp_export;

void export_info(void);

void export_info(void){

//________________________export last positions and velocities of sample_____________________________//
	
	
	sprintf(path,"%sresume.dat",pathout); 
	fp_export = fopen(path,"w+");
	fprintf(fp_export,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",titleB,npartB,xmin,xmax,ymin,ymax,zmin,zmax);
	fprintf(fp_export,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");	

	for(k=0;k<npartB;++k)
	{
		fprintf(fp_export,"%s\t%.16lf\t%.16lf\t%.16lf\t%d\t%d\t%.16lf\t%.16lf\t%.16lf\n",species_B[k],xfinalB[k],yfinalB[k],zfinalB[k],status[k],thermostatus[k],vxfinalB[k],vyfinalB[k],vzfinalB[k]);
	}
	fprintf(fp_export,"%d\n",n-1+nres);
	fclose(fp_export);

//------------------------------export last positions of tip----------------------------------------//
	
	sprintf(path,"%sresume_tip.dat",pathout);
	fp_export = fopen(path,"w+");
	fprintf(fp_export,"%d\n%s\n",npartA,titleA);
	fprintf(fp_export,"%lf\t%lf\t%lf\n",vx_tip,vy_tip,vz_tip);
	fprintf(fp_export,"%lf\t%lf\t%lf\n",xholder,yholder,zholder);
	
	for (k=0;k<npartA;++k)
	{		
		fprintf(fp_export,"%s\t%lf\t%lf\t%lf\n",species_A[k],xfinalA[k],yfinalA[k],zfinalA[k]);		
	}

	fclose(fp_export);
	

//_______________________________write to file energies.dat____________________________________________//	
	
	sprintf(path,"%sthermodynamic_quantities_%s%d_%d.dat",pathout,ensemble,(int)*pointer,sim_id);	
	fp_export=fopen(path,"w+");
	for(n=0;n<=N;++n)
	{
		fprintf(fp_export,"%d\t"
				"%.16lf\t"				
				"%.16lf\t"
				"%.16lf\t"
				"%.16lf\n",
				n+nres,
				K[n]/npartB,
				U[n]/npartB,
				E[n]/npartB,
				T[n]
				);
	}
	fclose(fp_export);

	
//________________________________write to file acf.dat_________________________________________________//
	
	sprintf(path,"%sacf_%s%d.dat",pathout,ensemble,(int)*pointer);
	fp_export = fopen(path,"w+");	
	for (n=0;n<=N-startpoint;++n)
	{
		fprintf(fp_export,"%d\t%.16lf\n",n,sumtime[n]);
	}
	fclose(fp_export);

	
	
//_______________________________write to file rdf.dat___________________________________________________//
	
	sprintf(path,"%srdf_%s%d.dat",pathout,ensemble,(int)*pointer);
	fp_export = fopen(path,"w+");	
	for (p=0;p<bins;++p)
	{
		fprintf(fp_export,"%.16lf\t%.16lf\n",r1[p],g[p]);
	}
	fclose(fp_export);	 	

} // end EXPORTS
