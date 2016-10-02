/*
 *  global.h
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 19/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 *
 */

#define cmax_length 1000
int sim_id; //simulation id
int nres; // for the multiple resume files
//------------------------------------------------ For Files----------------------------------------------------------//
char sim_folder[cmax_length]; // for defining folders
char pathcsv[cmax_length]; // for csv output
char pathxyz[cmax_length]; //for jmol output
char path[cmax_length],pathout[cmax_length],move_path[cmax_length],file_path[cmax_length]; //paths for exporting
char titleB[cmax_length]; //for surface
char titleA[cmax_length]; //for tip
char title[cmax_length]; // for jmol of tip A - sample B
char ensemble[cmax_length]; // for choosing the file extension NVT or NVE
double *pointer; //for naming files

//----------------------MD parameters for the sample---------------------------------------------------------------------//

int i,j,k,l,n,p,N,nstep; //general vars for indexing
int mobile,thermo; //vars for only sample (pop B)
int csvout,choose_thermostat,csv_color,monitor,jmolout,resume_surf,resume_tip,mypotential,aneal; //general flags

int *status; //choose if an atom is fixed or mobile
int *thermostatus;  //choose if an atom is thermostated or not

double sumkin,sumkinth; //sums for kinetic energy of mobile atoms or thermo atoms
double dtscale,dtSI,dt; //time scale and time step
double xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz; //simulation box

double *xB,*yB,*zB,*xfinalB,*yfinalB,*zfinalB;
double *vxB,*vyB,*vzB,*vxfinalB,*vyfinalB,*vzfinalB;
double *axB,*ayB,*azB,*axfinalB,*ayfinalB,*azfinalB;

double *sumfx, *sumfy, *sumfz; // forces between atoms of sample

double x,y,z,w,u; // for random number generator random_gen.c
double pi; // 3.14

//------------------ MD parameters that combine two different populations, e.g. tip-sample----------------------------------------//

int tip; //flag plugin for tip
int npartB,npartA; //atoms of 2 populations
int *Z_A,*Z_B; // atomic number for each population
char **species_A,**species_B; // atom species for each population
double *m_A,*covalent_r_A,*vdW_r_A,*m_B,*covalent_r_B,*vdW_r_B; //characteristics for different samples

double *xA,*yA,*zA,*xfinalA,*yfinalA,*zfinalA;
double *vxA,*vyA,*vzA,*vxfinalA,*vyfinalA,*vzfinalA;
double *axA,*ayA,*azA,*axfinalA,*ayfinalA,*azfinalA;

double *sumfxA, *sumfyA, *sumfzA; // forces of tip in tip-sample
double *sumfxB, *sumfyB, *sumfzB; // forces of sample in tip-sample
double sumforcesxA,sumforcesyA,sumforceszA,sumforcesxB,sumforcesyB,sumforceszB;
double *Tip_Forcex,*Tip_Forcey,*Tip_Forcez,*Atom_Forcex,*Atom_Forcey,*Atom_Forcez;

//------------------------------TIP MOVES-----------------------//
int tip_version;
int intervals,steps,tip_apex;  	//tip movement
double *xapex,*yapex,*zapex,*Dseper,*xtemp,*ytemp; //tip movement
double Ftargetx,FtargetxSI,Ftargety,FtargetySI,Ftargetz,FtargetzSI,Ftargetscale;
double gammax,gammay,gammaz,gammaxSI,gammaySI,gammazSI;
double kx,ky,kz,kxSI,kySI,kzSI,kscale; //tip movement
int graph_step;
double vx_tip,vy_tip,vz_tip,velocityscale;
double f0, Q, A0; //oscillations of cantilever
double Mtip;  // define tip mass
double xholder,yholder,zholder; //holder of cantilever tip

//-------------------------------Radial Distribution Function plus ACF---------------------------------------------------------//

int bins,myindex,startpoint,rdf_avg;
int *h,*hinst;
double *g,*r1,*ginst,*rinst; 
double dr,range;
double *vxfixed,*vyfixed,*vzfixed,*acf_t,*sumtime; //acf


//----------------------------------------- parameters for Thermostats------------------------------------------------//

int nbath,nbath_end,Tintervals; //general
double *Tinit,*Tfinal,*thermostat,Ttarget,Tcool,TempStep,kbol,Tinitial; //general
double number,nu,uxcm,uycm,uzcm,myrandom1,myrandom2; //Andersen
double constant,trise; //Berendsen

//-----------------------------------Energies----------------------------------------//

double *Vtotal; // LJ energy of total system (tip-sample)
double *U,*E,*K,*T,*Rep,*Attr,*Bond; // energies for sample	only


//------------------------------------------ LJ parameters------------------------------------------------------------//

// example: only Argon
double rc; //LJ cutoff 
double *myepsilon,*mysigma;
double epsilon,sigma;
double sumV;
double sumVcells;

// example: Silicon tip - Silicon surf
double rc_ts; //LJ cutoff
double *myepsilonA,*myepsilonB,*mysigmaA,*mysigmaB;
double epsilon_ts,sigma_ts;
double sumV_ts;
double sumV_ts_cells; 
//------------------------------------------Tersoff parameters--------------------------------------------------------//

double *beta,*ni,*ci,*di,*hi,xi; // undimensional parameters
double *myA,*myB; // eV
double *lamda,*mi; // 1/Angstrom
double *myR,*myS; // Angstrom

double rcmax; //the maximum rcut for linked cells when I have Tersoff multicomponent systems
double *dx_zeta_III,*dy_zeta_III,*dz_zeta_III;	//new terms, off-diagonal terms calc
double Repterm, Attrterm, Vtersoff;
double sumenergy,repsum,attrsum,bondsum;


double xforce,yforce,zforce,r2,r,xikforce,yikforce,zikforce,r2ik,rik;
//---------------------------------Linked Cells------------------------------------------------------------------------//

int linkcells; //if 1, then linked cells method is activated
int linkcells_supercell; //if 1, then linked cells method for supercell is activated
//simulation box for one or two populations
int cellsx,cellsy,cellsz,cellsx_supercell,cellsy_supercell,cellsz_supercell; 
double Lcx, Lcy, Lcz, Lcx_supercell, Lcy_supercell, Lcz_supercell;

int **cells,**neighbours,**cells_supercell,**neighbours_supercell;
int *HEAD,*LIST,*HEAD_A,*HEAD_B, *LIST_A,*LIST_B;

int *cellpartarray, *neighborcellpartarray, *allpartarray; //arrays for Tersoff linked-cells
int	*partindex,*nextpartindex; //arrays for LJ linked-cells
int *partindex_A,*partindex_B,*nextpartindex_B; //arrays for tip-sample LJ linked cells 
int cell,next;
