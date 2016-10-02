/*
 *  global.h
 *  MolDyn
 *
 *  Created by Dimitra Georgakaki on 22/03/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#define cmax_length 1000

//----------------------general parameters for MD---------------------------------------------------------------------

int i,j,k,l,n,p,N,npart,mobile,thermo,graph_step,nstep; //general
int csvout,jmolout,choose_thermostat,csv_color,monitor,resume_flag,mypotential,aneal,tip; //flags
int *status; //choose if an atom is fixed or mobile
int *thermostatus;  //choose if an atom is thermostated or not
double *m,*covalent_r,*vdW_r, rc, rcmax; 
int *Z;
char **species; 

double sumkin,sumkinth;
double *ax1,*ay_1,*az1; // initial accel
double *vx1,*vy_1,*vz1; //initial velocities
double *x1,*y_1,*z1; // initial positions

double *x2,*y2,*z2; //atomic positions
double *vx2,*vy2,*vz2,*ax2,*ay2,*az2; //velocities, accelerations
double dtscale,dtSI,dt; //time scale and time step
double xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz; //simulation box
double *sumfx, *sumfy, *sumfz; // forces

double x,y,z,w,u; // for random number generator

//------------------------------------------------ For Files----------------------------------------------------------

char pathcsv[cmax_length]; // for csv output
char pathxyz[cmax_length]; //for jmol output
char title[cmax_length]; // for jmol format
char myexport[cmax_length]; // for jmol/paraview export
char sim_folder[cmax_length]; // for defining folders
char path[cmax_length]; // for all the exports
char pathout[cmax_length]; // for the \\out folder in the current sim folder
char ensemble[cmax_length]; // for choosing the file extension NVT or NVE
double *pointer; // for naming files

char file_path[cmax_length]; //for config.txt

//-------------------------------Radial Distribution Function plus ACF------------------------------------------------

int bins,myindex,startpoint,rdf_avg;
int *h,*hinst;
double *g,*r1,*ginst,*rinst; 
double dr,range;

double *vxfixed,*vyfixed,*vzfixed,*acf_t,*sumtime; //acf


//----------------------------------------- parameters for Thermostats------------------------------------------------

int nbath,nbath_end,Tintervals;
double number,nu,uxcm,uycm,uzcm,myrandom1,myrandom2; //Andersen
double constant,trise; //Berendsen


double *Tinit,*Tfinal,*thermostat,Ttarget,Tcool,Tinitial,TempStep,kbol; //general

//------------------------------------------ LJ parameters------------------------------------------------------------

double *myepsilon,*mysigma,epsilon,sigma,sumV;

//------------------------------------------Tersoff parameters--------------------------------------------------------

double *beta,*ni,*ci,*di,*hi,xi; // undimensional parameters
double *myA,*myB; // eV
double *lamda,*mi; // 1/Angstrom
double *myR,*myS; // Angstrom
 
double r2,r,r2ik,rik,pi,xforce,yforce,zforce,xikforce,yikforce,zikforce;
double Vtersoff,Repterm,Attrterm,sumenergy;
double *Rep,*Attr,*Bond; // Tersoff terms
double *U,*E,*K,*T; // total energies

double *dx_zeta_III,*dy_zeta_III,*dz_zeta_III;	//new terms, off-diagonal terms calc

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