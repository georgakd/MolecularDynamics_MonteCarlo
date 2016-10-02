/*
 *  global.h
 *  tersoff calculations
 *
 *  Created by Dimitra Georgakaki on 02/09/2011.
 *  Copyright 2011 Home. All rights reserved.
 *
 */
#define cmax_length 1000
//----------------------general parameters for MC---------------------------------------------------------------------

int i,j,k,l,n,nsample,p,npart,mobile,thermo,myrandom; //general
int csvout,jmolout,csv_color,mypotential,resume_flag; //flags
int *status; //choose if an atom is fixed or mobile
int *thermostatus;  //choose if an atom is thermostated or not

double *m,*covalent_r,*vdW_r, rc, rcmax; 
int *Z;
char **species; 
int temp; // random particle
double *x1,*y_1,*z1,*vx,*vy,*vz; // initial positions
double *x2,*y2,*z2; //atomic positions after the update
double xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz; //simulation box
double dmax; // for the maximum allowed displacement
int max_cycles,MC_accept;
double sumU;
double myrandom_choice,myrandom_choicex,myrandom_choicey,myrandom_choicez,random_dispx,random_dispy,random_dispz;
int counter; // for counting the accepted MC steps 
double w,x,y,z,u,Utrial,Uold,dU;
double xtemp,ytemp,ztemp;

//-------------------------------Thermostats--------------------------------------------------------------------------

double Tref,kbol,Wvar;


//------------------------------------------------ For Files-----------------------------------------------------------
char file_path[cmax_length];
char myexport[cmax_length];
char pathcsv[cmax_length]; // for csv output
char pathxyz[cmax_length]; //for jmol output
char title[cmax_length]; // for jmol format

//------------------------------------------ LJ parameters------------------------------------------------------------

double *myepsilon,*mysigma,epsilon,sigma,sumV;

//------------------------------------------Tersoff parameters--------------------------------------------------------

double *beta,*ni,*ci,*di,*hi,xi; // undimensional parameters
double *myA,*myB; // eV
double *lamda,*mi; // 1/Angstrom
double *myR,*myS; // Angstrom
 
double r2,r,r2ik,rik,pi,xforce,yforce,zforce,xikforce,yikforce,zikforce;
double *U,sumenergy; // energies	

//------------------------------------------Sutton-Chen parameters----------------------------------------------------

double npower,mpower;
double c_param,term,coeff_rep,coeff_attr,sum_attr,denfunc,sum_rep;


//---------------------------------Linked Cells------------------------------------------------------------------------
int icell,totalcells,line,myindex,index1,index2,box,myn;
int cellsx,cellsy,cellsz,ghostsx,ghostsy,cell,linkcells; 
double Lcx, Lcy, Lcz,celli,sumVcells,rcells;
int *HEAD,*LIST,*partindex,*nextpartindex; 
int HEADp,LISTp; 
int **cells,**neighbours ;
int counter,nextcounter,allcellcounter,next;
int locali,localj,localk;
int *cellpartarray, *neighborcellpartarray, *allpartarray; 
int cellpart, allcellpart,neighborcellpart; // for tersoff cells
