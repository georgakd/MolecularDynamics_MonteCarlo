/*
 *  sim_aneal.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 23/11/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void sim_aneal(void);

void sim_aneal(void){
	
	


//*********Here starts the simulated annealing or simulated quenching if chosen******************//		

	if (nbath_end%Tintervals==0){printf("\nThe Temperature intervals are checked and found ok\n");}
	else
	{	printf("Aborting Simulation because the temperature intervals are not ok...Press 1 for exit...\n");
		system("PAUSE");
		exit(1);}
	
if (aneal == 1) //simulated annealing
{
	for (i=0;i<Tintervals;++i)
	{
		Tinit[i]=(i)*TempStep + Tinitial;
		Tfinal[i]=(i+1)*TempStep + Tinitial;
		
		for (n = i*((nbath_end)/Tintervals)+nbath; n <= (i+1)*((nbath_end)/Tintervals)+nbath; ++n){ // we add +1 in order to start the division of segments from n=1 and not from n=0
			if (n>=nbath && n<=nbath_end){	
				
				thermostat[n]=Tfinal[i];
				
				
			}
		}	
	
	}
}

	
	
if (aneal == -1) //simulated quenching
{
	for (i=0;i<Tintervals;++i)
	{
		Tinit[i] = (Tintervals-i)*(-TempStep)+Tcool;

		Tfinal[i] = (Tintervals-(i+1))*(-TempStep)+Tcool;
		
		
		for (n = i*((nbath_end)/Tintervals)+nbath; n <= (i+1)*((nbath_end)/Tintervals)+nbath; ++n){ // we add +1 in order to start the division of segments from n=1 and not from n=0
			if (n>=nbath && n<=nbath_end){	
				
				thermostat[n]=Tfinal[i];				
				
			}
		}		

	}
}
	//for(n=nbath;n<=nbath_end;++n){printf("%lf\n",thermostat[n]);}
	//printf("%lf\t%lf",Tinit[0],Tfinal[Tintervals-1]);
	//getchar();
	

}