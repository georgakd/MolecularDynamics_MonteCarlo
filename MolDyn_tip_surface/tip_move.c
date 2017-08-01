/*
 *  tip_move.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 25/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

FILE *fp;

void tip_move(void);

void tip_move(void)
{

//............................................ Choose how to move the tip: (1) No feedback (constant-height) ...................................... //

//............................................ Choose how to move the tip: (2) Feedback (constant-force) ...................................... //
		int horiz_iter;
		
		fp=fopen(move_path,"r");
		if(fp==NULL){printf("Could not locate %s",move_path);printf("\nAborting Simulation...Press 1 for exit...\n");
		exit(1);}
		fscanf(fp,"tip: %d\n",&tip);
		fscanf(fp,"tip_version: %d\n",&tip_version); // 1 dynamic surface, 2 static surface
		fscanf(fp,"intervals: %d\n",&intervals);
		fscanf(fp,"horiz_iteration: %d\n",&horiz_iter);
		fscanf(fp,"FtargetxSI: %lf\n",&FtargetxSI);
		fscanf(fp,"FtargetySI: %lf\n",&FtargetySI);
		fscanf(fp,"FtargetzSI: %lf\n",&FtargetzSI);
		fscanf(fp,"vx_tip: %lf\n",&vx_tip);
		fscanf(fp,"vy_tip: %lf\n",&vy_tip);
		fscanf(fp,"vz_tip: %lf\n",&vz_tip);
		fscanf(fp,"kxSI: %lf\n",&kxSI);
		fscanf(fp,"kySI: %lf\n",&kySI);
		fscanf(fp,"kzSI: %lf\n",&kzSI);
		fscanf(fp,"gammaxSI: %lf\n",&gammaxSI);
		fscanf(fp,"gammaySI: %lf\n",&gammaySI);
		fscanf(fp,"gammazSI: %lf\n",&gammazSI);
		fscanf(fp,"f0: %lf\n",&f0);
		fscanf(fp,"Q: %lf\n",&Q);
		fclose(fp);
	
}
