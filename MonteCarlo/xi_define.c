/*
 *  xi_define.c
 *  tersoff calculations
 *
 *  Created by Dimitra Georgakaki on 23/01/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
 
  
double xi_define(int myZ1, int myZ2);

double xi_define(int myZ1, int myZ2)
{
	double myxi_SiC = 0.9776; // SiC parameter
	//double myxi_SiGe = 1.00061; //SiGe parameter
	
	if (myZ1 == 6 && myZ2 == 14)
	{	
		return myxi_SiC;
	}
	
	else if (myZ1 == 14 && myZ2 == 6)
	{
		return myxi_SiC;
	}
	/*
	if (myZ1 == 32 && myZ2 == 14)
	{	
		return myxi_SiGe;
	}
	
	else if (myZ1 == 14 && myZ2 == 32)
	{
		return myxi_SiGe;
	}
	*/

	
	else 
	{
		return 1.0;
	}


}
