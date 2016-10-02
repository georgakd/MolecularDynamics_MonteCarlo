#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"

//--------------Implementation of the Wichmann - Hill II Generator

double random_gen(void);

double random_gen(void){

int intpart; //take the integer part of the random number
double decpart; //take the decimal part of the random number
int m1 = 2147483579;
int	m2 = 2147483543;
int	m3 = 2147483423;
int	m4 = 2147483123;
int a1 = 11600; 
int	a2 = 47003; 
int	a3 = 23000;
int	a4 = 33000;


//--------------------------------------------------------------------------------------------


		w = a1*((int)w%185127)-10379*((int)w/185127);
		x = a2*((int)x%45688)-10479*((int)x/45688);
		y = a3*((int)y%93368)-19423*((int)y/93368);
		z = a4*((int)z%65075)-8123*((int)z/65075);

		if (w<0){w=w+m1;}
		if (x<0){x=x+m2;}
		if (y<0){y=y+m3;}
		if (z<0){z=z+m4;}

		u = w*1.0/m1 + x*1.0/m2 + y*1.0/m3 + z*1.0/m4;
		intpart=(int)u;
		decpart = u - intpart;
		u=decpart;
	
	return u;

}