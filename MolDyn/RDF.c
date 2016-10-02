/*
 *  RDF.c
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 10/02/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */


#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void RDF(void);

void RDF(void)
{
	//---------------------------- Radial Distribution Function-----------------------------------------------------

	
	for (p=0;p<bins;++p) //bin number
	{
		r1[p] = ((p+1)-0.5)*dr;
		//g[p] = (1.0 / rdf_avg )*((Lx*Ly*Lz*h[p])/(4*pi*r1[p]*r1[p]*dr*npart*npart)); // we average over time for smoothness ex: for the last N-rdf_avg = 1000 time steps
		g[p] = (1.0 / rdf_avg )*((Lx*Ly*Lz*h[p])/(4*pi*r1[p]*r1[p]*dr*mobile*mobile)); //rdf only on mobile atoms
	}
	

}