
//**************************************** MOLDYN_TIP_SURFACE PROJECT ******************************************//
//		RELEASE MODE v.3 16/03/2013																					//
//		Created by Dimitra Georgakaki																			//																					
//																												//
//																												//	
// This project calculates the tip-sample forces at a fixed tip point (x,y,z) or at moving tip points.			//
// The routine tip_move is responsible for the tip movement with or without feedback.							//
//																												//
//**************************************************************************************************************//

#define cmax_length 1000
#ifdef _WIN32
#define OS 1
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#include<direct.h>
#elif _WIN64
#define OS 1
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#include<direct.h>
#elif __linux__
#define OS 2
#include<unistd.h>
#elif __APPLE__
#define OS 3
#include<unistd.h>
#elif __MACH__
#define OS 3
#include<unistd.h>
#else
#define OS 0
#include<unistd.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include <time.h>


FILE *fp,*fp1,*fp2; //file pointers to read and write data

void parameters(void); // declare function to initialize my parameters
void atomic_numbers(void); 
void read_elements(void);
void sim_aneal(void);
double random_gen(void); // include random number generator
void export_info(void); // export info
void tip_move(void);  //choose how to move the tip

void linked_cells(int mycellsx, int mycellsy, int mycellsz); //for sample
void linked_cells_supercell(int mycellsx,int mycellsy, int mycellsz); //for tip-sample
void Andersen(double *vx, double *vy, double *vz); // declare the function that initiates the Andersen thermostat
void vrescale(double *vx, double *vy, double *vz); // declare the function that initiates the Berendsen thermostat
void RDF(void); // Radial Distribution Function
void acf(double *vx,double *vy,double *vz); // ACF

//-------- declare the function that will calculate interaction forces and potential_____________________//
void tersoff_forces_cells(double *x,double *y,double *z);
void tersoffopt(double *x,double *y,double *z);
void lj_forces(double *x,double *y,double *z); 
void lj_forces_cells(double *x,double *y,double *z);
//tip-sample
void lj_ts_forces(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);
void lj_ts_forces_cells(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);

//-------------------------graphics---------------------------------------------------------//

void jmol_func(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);
void paraview_func(double *xA, double *xB,double *yA, double *yB,double *zA, double *zB);



int main(int argc, char *argv[])
{
	
//*************************************************** Define local vars / parameters ********************************************************//
	
	
	double sumvx,sumvy,sumvz; // local sums for velocity center of mass at n=0
	double sumx, sumy, sumz,sumxfinal,sumyfinal,sumzfinal; //local sums for tip center of mass
	char pathf[cmax_length]; // for monitor instant forces in matlab
	char pathrdf[cmax_length]; // for monitor instant rdf in matlab
	char path_animation[cmax_length]; //for jmol video
	char mycommand[cmax_length]; // for mkdir and delete files and folders
	char starttime[cmax_length],stoptime[cmax_length];	// for counting the simulation	
	
	double AvgForcex, AvgForcey, AvgForcez;
	double TotalForcex, TotalForcey, TotalForcez;
	double Full_Fx,Full_Fy,Full_Fz;
	double Dseper;
	
	time_t rawtime;
	time(&rawtime);
	sprintf (starttime, "\nSimulation started at: %s",ctime(&rawtime));


	if (argc==1){printf("Define simulation id...\n");exit(-1);}
	sim_id=atoi(argv[1]);

//--------------------------------------------Read current directory--------------------------------------//
	
if (OS==1){_getcwd(sim_folder,cmax_length);}
	else {getcwd(sim_folder,cmax_length);}


//------------------------------------------------Creating folders, paths--------------------------------------------------------------//

	
	if (OS==1){
		
		sprintf(file_path,"%s\\config.txt",sim_folder);		
		sprintf(move_path,"%s\\move.txt",sim_folder);
		sprintf(pathcsv,"%s\\csv\\",sim_folder);
		sprintf(pathout,"%s\\out\\",sim_folder);
		sprintf(pathf,"%s\\forces\\",sim_folder);
		sprintf(pathrdf,"%s\\radial\\",sim_folder);
		sprintf(pathxyz,"%s\\jmol\\",sim_folder);
		sprintf(path_animation,"%s\\out\\animation%d.xyz",sim_folder,sim_id);
	}
	
	else { 
		
		sprintf(file_path,"%s/config.txt",sim_folder);
		sprintf(move_path,"%s/move.txt",sim_folder);
		sprintf(pathcsv,"%s/csv/",sim_folder);
		sprintf(pathout,"%s/out/",sim_folder);
		sprintf(pathf,"%s/forces/",sim_folder);
		sprintf(pathrdf,"%s/radial/",sim_folder);
		sprintf(pathxyz,"%s/jmol/",sim_folder);
		sprintf(path_animation,"%s/out/animation%d.xyz",sim_folder,sim_id);
	}
	
	
	sprintf(mycommand,"mkdir %s",pathcsv);
	system(mycommand);
	sprintf(mycommand,"mkdir %s",pathout);
	system(mycommand);
	sprintf(mycommand,"mkdir %s",pathf);
	system(mycommand);
	sprintf(mycommand,"mkdir %s",pathrdf);
	system(mycommand);	
	sprintf(mycommand,"mkdir %s",pathxyz);
	system(mycommand);
	
	
	printf("\n\n--------------Define the simulation folder: ");printf("%s",sim_folder);printf("-----\n\n");

	
	parameters();  // call the function parameters.c
	tip_move();  //Import tip parameters
	
	
//----------------------------------------------Delete previous folders csv,jmol,radial,forces if they exist----------------------------------------//
	
		
	if (OS==1){
		sprintf(mycommand,"del /Q %s*.csv",pathcsv);system(mycommand);
		sprintf(mycommand,"del /Q %s*.xyz",pathxyz);system(mycommand);
	}
	else{
		sprintf(mycommand,"find %s -name '*.csv' | xargs rm",pathcsv);system(mycommand);
		sprintf(mycommand,"find %s -name '*.xyz' | xargs rm",pathxyz);system(mycommand);
	}	
	

//--------------------------------------------Read box dimensions, number of particles, init surface positions-----------------------------------//

	if (resume_surf == 0 && resume_tip == 0) {nres=0;} // set flag 0 for multiple files
	
//----------------------Surface information--------------------------//

if (resume_surf == 0){ 	// import for the first time initial positions	


		if (OS==1){sprintf(path,"%s\\positions.dat",sim_folder);}
		else {sprintf(path,"%s/positions.dat",sim_folder);}
		fp=fopen(path,"r");
		fscanf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",titleB,&npartB,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
		fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\n");
		xB=(double*)malloc((npartB)*sizeof(double));
		yB=(double*)malloc((npartB)*sizeof(double));
		zB=(double*)malloc((npartB)*sizeof(double));
		vxB=(double*)malloc((npartB)*sizeof(double));
		vyB=(double*)malloc((npartB)*sizeof(double));
		vzB=(double*)malloc((npartB)*sizeof(double));			
		Z_B=(int*)malloc((npartB)*sizeof(int));
		m_B=(double*)malloc((npartB)*sizeof(double));
		covalent_r_B=(double*)malloc((npartB)*sizeof(double));
		vdW_r_B=(double*)malloc((npartB)*sizeof(double));
		status=(int*)malloc((npartB)*sizeof(int));
		thermostatus=(int*)malloc((npartB)*sizeof(int));
		species_B=(char**)malloc((npartB)*sizeof(char*));		
		for (k = 0; k < npartB; ++k){species_B[k]=(char*)malloc((3)*sizeof(char));}
		for (k=0;k<npartB;++k){
			fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species_B[k],&xB[k],&yB[k],&zB[k],&status[k],&thermostatus[k]);
		}
		fclose(fp);
	}
	
if (resume_surf == 1){ 	// import resume surface file	
		if (OS==1){sprintf(path,"%s\\resume.dat",pathout);}
		else {sprintf(path,"%s/resume.dat",pathout);}
		fp=fopen(path,"r");
		fscanf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",titleB,&npartB,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
		fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");
		xB=(double*)malloc((npartB)*sizeof(double));
		yB=(double*)malloc((npartB)*sizeof(double));
		zB=(double*)malloc((npartB)*sizeof(double));
		vxB=(double*)malloc((npartB)*sizeof(double));
		vyB=(double*)malloc((npartB)*sizeof(double));
		vzB=(double*)malloc((npartB)*sizeof(double));			
		Z_B=(int*)malloc((npartB)*sizeof(int));
		m_B=(double*)malloc((npartB)*sizeof(double));
		covalent_r_B=(double*)malloc((npartB)*sizeof(double));
		vdW_r_B=(double*)malloc((npartB)*sizeof(double));
		status=(int*)malloc((npartB)*sizeof(int));
		thermostatus=(int*)malloc((npartB)*sizeof(int));
		species_B=(char**)malloc((npartB)*sizeof(char*));		
		for (k = 0; k < npartB; ++k){species_B[k]=(char*)malloc((3)*sizeof(char));}
		for (k=0;k<npartB;++k){
			fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",species_B[k],&xB[k],&yB[k],&zB[k],&status[k],&thermostatus[k],&vxB[k],&vyB[k],&vzB[k]);
		}
		fscanf(fp,"%d\n",&nres);
		fclose(fp);
	}	


//------------------------Tip information-----------------------------------//

	//index A represents the initial tip positions
if (resume_tip == 0){
	if (OS==1){sprintf(path,"%s\\tip.dat",sim_folder);}
	else {sprintf(path,"%s/tip.dat",sim_folder);}
	fp=fopen(path,"r");
	fscanf(fp,"%d\n%s\n",&npartA,titleA);	
	xA=(double*)malloc((npartA)*sizeof(double)); //initial positions
	yA=(double*)malloc((npartA)*sizeof(double));
	zA=(double*)malloc((npartA)*sizeof(double));
	Z_A=(int*)malloc((npartA)*sizeof(int));
	m_A=(double*)malloc((npartA)*sizeof(double));
	covalent_r_A=(double*)malloc((npartA)*sizeof(double));
	vdW_r_A=(double*)malloc((npartA)*sizeof(double));
	species_A=(char**)malloc((npartA)*sizeof(char*));		
	for (k = 0; k < npartA; ++k)
	{species_A[k]=(char*)malloc((3)*sizeof(char));}
	
	for (k=0;k<npartA;++k)
	{		
		fscanf(fp,"%s\t%lf\t%lf\t%lf\n",species_A[k],&xA[k],&yA[k],&zA[k]);		
	}
	fclose(fp);
	}
	
if (resume_tip == 1){
		if (OS==1){sprintf(path,"%s\\resume_tip.dat",pathout);}
		else {sprintf(path,"%s/resume_tip.dat",pathout);}
		fp=fopen(path,"r");
		fscanf(fp,"%d\n%s\n",&npartA,titleA);
		fscanf(fp,"%lf\t%lf\t%lf\n",&vx_tip,&vy_tip,&vz_tip);
		fscanf(fp,"%lf\t%lf\t%lf\n",&xholder,&yholder,&zholder);
		xA=(double*)malloc((npartA)*sizeof(double)); //initial positions
		yA=(double*)malloc((npartA)*sizeof(double));
		zA=(double*)malloc((npartA)*sizeof(double));
		Z_A=(int*)malloc((npartA)*sizeof(int));
		m_A=(double*)malloc((npartA)*sizeof(double));
		covalent_r_A=(double*)malloc((npartA)*sizeof(double));
		vdW_r_A=(double*)malloc((npartA)*sizeof(double));
		species_A=(char**)malloc((npartA)*sizeof(char*));		
		for (k = 0; k < npartA; ++k)
		{species_A[k]=(char*)malloc((3)*sizeof(char));}
		
		for (k=0;k<npartA;++k)
		{		
			fscanf(fp,"%s\t%lf\t%lf\t%lf\n",species_A[k],&xA[k],&yA[k],&zA[k]);		
		}
		fclose(fp);
	}


//-------------------- Set tip movement,  import initial apex coords  --------------------------------------------//
	
	
	// if tip apex is not obvious you should choose the apex atom from jmol or mathematica
	
	Dseper = zA[tip_apex] - zB[npartB-1]; // the Dseperation can only be calculated for a static flat surface!
	printf("\nThe initial tip-sample seperation is:  %lf", Dseper);
	printf("\nThe initial tip-height is:  %lf", zA[tip_apex]);

	
//------------------------------------------------------- Allocation of all other arrays -----------------------------------------------------------//	

	xfinalB=(double*)malloc((npartB)*sizeof(double));
	yfinalB=(double*)malloc((npartB)*sizeof(double));
	zfinalB=(double*)malloc((npartB)*sizeof(double));
	xfinalA=(double*)malloc((npartA)*sizeof(double));
	yfinalA=(double*)malloc((npartA)*sizeof(double));
	zfinalA=(double*)malloc((npartA)*sizeof(double));
	
	vxfinalB=(double*)malloc((npartB)*sizeof(double));
	vyfinalB=(double*)malloc((npartB)*sizeof(double));
	vzfinalB=(double*)malloc((npartB)*sizeof(double));
		
	axB=(double*)malloc((npartB)*sizeof(double)); //initial accelerations
	ayB=(double*)malloc((npartB)*sizeof(double));
	azB=(double*)malloc((npartB)*sizeof(double));
	axfinalB=(double*)malloc((npartB)*sizeof(double)); //final accelerations
	ayfinalB=(double*)malloc((npartB)*sizeof(double));
	azfinalB=(double*)malloc((npartB)*sizeof(double));
	
	//sums for LJ interactions
	sumfxA=(double*)malloc((npartA)*sizeof(double)); 
	sumfyA=(double*)malloc((npartA)*sizeof(double));
	sumfzA=(double*)malloc((npartA)*sizeof(double));
	sumfxB=(double*)malloc((npartB)*sizeof(double)); 
	sumfyB=(double*)malloc((npartB)*sizeof(double));
	sumfzB=(double*)malloc((npartB)*sizeof(double));
	
	//Tip-Sample Forces		
	Tip_Forcex=(double*)malloc((N+1)*sizeof(double)); // x-tip-force in time
	Tip_Forcey=(double*)malloc((N+1)*sizeof(double)); // y-tip-force in time
	Tip_Forcez=(double*)malloc((N+1)*sizeof(double)); // z-tip-force in time	
	
	//Thermostat arrays
	thermostat=(double*)malloc((N+1)*sizeof(double));
	Tinit=(double*)malloc(Tintervals*sizeof(double));
	Tfinal=(double*)malloc(Tintervals*sizeof(double));
	
	//Total Energies of the given system
	U=(double*)malloc((N+1)*sizeof(double));
	K=(double*)malloc((N+1)*sizeof(double));
	E=(double*)malloc((N+1)*sizeof(double));
	T=(double*)malloc((N+1)*sizeof(double));		
	Rep=(double*)malloc((N+1)*sizeof(double));
	Attr=(double*)malloc((N+1)*sizeof(double));
	Bond=(double*)malloc((N+1)*sizeof(double));
	
	Vtotal=(double*)malloc((N+1)*sizeof(double)); // potential of the total system in time	
	
	//velocity acf	
	vxfixed=(double*)malloc((npartB)*sizeof(double));
	vyfixed=(double*)malloc((npartB)*sizeof(double));
	vzfixed=(double*)malloc((npartB)*sizeof(double));
	acf_t=(double*)malloc((npartB)*sizeof(double)); // particles sum
	sumtime=(double*)malloc((N+1-startpoint)*sizeof(double)); // time sum
	
	// rdf and rdf monitor
	bins = (int)floor(range/dr+0.5); // we add 0.5 because we want the function round (c has round but too complicated)
	h=(int*)malloc((bins)*sizeof(int));	
	g=(double*)malloc((bins)*sizeof(double));
	r1=(double*)malloc((bins)*sizeof(double));
	for (p=0;p<bins;++p) //initialize
	{
		g[p]=0;r1[p]=0;h[p]=0;
	}	
	hinst=(int*)malloc((bins)*sizeof(int));	
	ginst=(double*)malloc((bins)*sizeof(double));
	rinst=(double*)malloc((bins)*sizeof(double));
	for (p=0;p<bins;++p) //initialize
	{
		ginst[p]=0;rinst[p]=0;hinst[p]=0;
	}
	
	
//Be careful, forces of sample are two-dim arrays but they can be transformed to 1D!
	
	sumfx=(double*)malloc((npartB)*sizeof(double)); //forces 1D (for speed)
	sumfy=(double*)malloc((npartB)*sizeof(double));
	sumfz=(double*)malloc((npartB)*sizeof(double));		
	

//---------------------------------Define secondary parameters and import initial crystal state-------------------------------------------------//
	
	//***********be very careful this is the correct order for the parameters to be read from the program*******************//
	
	
	if (mypotential==0){printf("\n  Tersoff T3 potential...\n");}
	if (mypotential==1){printf("\n  Lennard Jones potential...\n");}
	
	mobile=0;thermo=0;
	for (k=0;k<npartB;++k){mobile=mobile+status[k];} //count the mobile atoms
	for (k=0;k<npartB;++k){thermo=thermo+thermostatus[k];} //count the thermostated atoms
	
	Lx=xmax-xmin; // box simulation
	Ly=ymax-ymin; 
	Lz=zmax-zmin;
	dt=dtSI/dtscale;  // the reduced time step used in the simulation
	Ftargetx = FtargetxSI/Ftargetscale; //the reduced force x used in the simulation
	Ftargety = FtargetySI/Ftargetscale; //the reduced force y used in the simulation
	Ftargetz = FtargetzSI/Ftargetscale; //the reduced force z used in the simulation
	kx = kxSI/kscale;ky = kySI/kscale;kz = kzSI/kscale;	//the reduced spring constant used in the simulation
	gammax = (gammaxSI*1000000)/(kscale*dtscale);
	gammay = (gammaySI*1000000)/(kscale*dtscale);
	gammaz = (gammazSI*1000000)/(kscale*dtscale); //the reduced hydrodynamic term used in the simulation in nN*sec/m
	printf("The reduced Ftargetx is %lf\n",Ftargetx);
	printf("The reduced Ftargety is %lf\n",Ftargety);
	printf("The reduced Ftargetz is %lf\n",Ftargetz);
	printf("The reduced kz is %lf\n",kz);
	printf("The reduced gammaz is %.16lf\n",gammaz);

	TempStep = Ttarget/Tintervals; //for the simulated anealing/quenching
	
	sim_aneal(); // call the function for anealing or quenching
	
	if (choose_thermostat==0 || nbath_end!=N) // CASE 1: NVE ensemble
	{sprintf(ensemble,"%s","NVE");pointer=&Ttarget;}
	
	if (choose_thermostat==1) // CASE 2: NVT with Andersen 
	{sprintf(ensemble,"%s","NVT_And");pointer=&Ttarget;}
	
	if (choose_thermostat==2 && aneal==1 && Tintervals==1) // CASE 3: NVT with Berendsen
	{sprintf(ensemble,"%s","NVT_Ber");pointer=&Ttarget;}
	
	if (choose_thermostat==2 && aneal==1 && Tintervals!=1) //CASE 4: Simulated Annealing
	{sprintf(ensemble,"%s","sim_anealing");pointer = &Tfinal[Tintervals-1];}
	
	if (choose_thermostat==2 && aneal==-1)  //CASE 5: Simulated Cooling
	{sprintf(ensemble,"%s","sim_cooling");pointer=&Tcool;}
	
	atomic_numbers();
	read_elements();

//-------------Tip calculation of mass center before starting the simulation-------------------//
		
	Mtip = npartA*m_A[0]; //total mass of tip in amu
	printf("The Tip Mass is %.lf\n",Mtip);

	if (resume_tip==0){
	sumx=0;sumy=0;sumz=0;
	for (i=0;i<npartA;++i) {
	
		sumx = sumx + m_A[i]*xA[i];
		sumy = sumy + m_A[i]*yA[i];
		sumz = sumz + m_A[i]*zA[i];
		
	}
	xholder = sumx /Mtip; yholder = sumy /Mtip; zholder = sumz/Mtip;
	}

	if (xholder>xmax||xholder<xmin||yholder>ymax||yholder<ymin){printf("The coordinates of the holder are: ");printf("%lf\t%lf\t%lf\n",xholder,yholder,zholder);printf("\nExceed box...Press 1 for exit...\n");exit(1);}
	else {printf("The coordinates of the holder are: ");printf("%lf\t%lf\t%lf\n",xholder,yholder,zholder);}
	
	A0 = (Dseper*kz) / Q; //calculation of initial cantilever amplitude
	
	
	
//--------------------------------------------------- Use or not use LINKED CELLS -----------------------------------------------------------------//

	if (linkcells==1){
		// 1. Divide space for Tersoff sample-sample interactions
		printf("\n\n--------------------------LINKED CELLS----------------------------------- \n");
		printf("Dividing space for sample.......\n");	
		printf("Give me the maximum rcut for the given system: rcut_max = "); // for a multicomponent system i choose the max S
		scanf("%lf",&rcmax);		
		
		Lcx=Lx/cellsx;Lcy=Ly/cellsy;Lcz=Lz/cellsz;
		
		if (Lcx >= rcmax)
			printf("The division of cells in x direction is: %d  %lf %lf\n",cellsx, Lcx, rcmax);
		else 
			printf("Sorry, try again!\n");
		if (Lcy >= rcmax)
			printf("The division of cells in y direction is: %d %lf %lf\n",cellsy, Lcy, rcmax);
		else 
			printf("Sorry, try again!\n"); 
		if (Lcz >= rcmax)
			printf("The division of cells in z direction is: %d %lf %lf\n\n",cellsz, Lcz, rcmax);
		else 
			printf("Sorry, try again!\n\n");
				
		linked_cells(cellsx,cellsy,cellsz); //call the function for building Tersoff neighbors
	}

	if (linkcells_supercell==1){

		// 2. Divide space for LJ supercell interactions

		printf("The tip is on! Dividing supercell space.......\n");	

		Lcx_supercell=Lx/cellsx_supercell;Lcy_supercell=Ly/cellsy_supercell;Lcz_supercell=Lz/cellsz_supercell;
		
		if (Lcx_supercell >= 2.5*mysigmaA[0])
			printf("The division of cells in x direction is: %d  %lf %lf\n",cellsx_supercell, Lcx_supercell, 2.5*mysigmaA[0]);
		else 
			printf("Sorry, try again!\n");
		if (Lcy_supercell >= 2.5*mysigmaA[0])
			printf("The division of cells in y direction is: %d %lf %lf\n",cellsy_supercell, Lcy_supercell, 2.5*mysigmaA[0]);
		else 
			printf("Sorry, try again!\n"); 
		if (Lcz_supercell >= 2.5*mysigmaA[0])
			printf("The division of cells in z direction is: %d %lf %lf\n\n",cellsz_supercell, Lcz_supercell, 2.5*mysigmaA[0]);
		else 
			printf("Sorry, try again!\n\n");
				
		linked_cells_supercell(cellsx_supercell,cellsy_supercell,cellsz_supercell); //call the function for building LJ neighbors

}//end linked-cells		

		printf("------------------------------------------------------------------------- \n");
		
	

//***************************************All initializations done, start the program!!!*********************************************************************//	

			
	printf("All initializations done, start the program!!!\n");
	printf("\nTotal atoms of sample: %d\n", npartB);printf("The mobile atoms are: %d\n", mobile);printf("The thermostated atoms are: %d\n", thermo);
	printf("Total atoms of tip: %d\n", npartA);printf("Total supercell atoms: %d\n\n", npartA+npartB);


// Initial Simulation n=0: Rescaling of velocities drawn from a gaussian distribution
	n=0;	
	TotalForcex=0; //initialize the calculation of Total tip-sample Forces in x y z directions
	TotalForcey=0;
	TotalForcez=0;
	
	if (resume_surf==0){ // the first run is always with random velocities


		sumvx=0;sumvy=0;sumvz=0;
		
		for (l=0;l<npartB;++l)
		{
			if (status[l]==0){vxB[l]=0;vyB[l]=0;vzB[l]=0;}
			
			else if (status[l]==1)
			{
				myrandom1=random_gen();
				myrandom2=random_gen();
					//myrandom1=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);		// Box-Muller method to create gaussian random numbers
					//myrandom2=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				vxB[l]=sqrt((kbol*Ttarget)/m_B[l])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2); // Box-Muller Method
				vyB[l]=sqrt((kbol*Ttarget)/m_B[l])*sqrt(-2*log(myrandom1))*sin(2*pi*myrandom2);
			
				myrandom1=random_gen();
				myrandom2=random_gen();

					//myrandom1=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);		// Box-Muller method to create gaussian random numbers
					//myrandom2=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				vzB[l]=sqrt((kbol*Ttarget)/m_B[l])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2);
			
			
			sumvx=sumvx+vxB[l];sumvy=sumvy+vyB[l];sumvz=sumvz+vzB[l];  //velocity center of mass for mobile particles
			
			}
			
		}	
		
		sumvx=sumvx/mobile; sumvy=sumvy/mobile; sumvz=sumvz/mobile;  //   <v>

		sumkin=0;
		sumkinth=0; //if I want to calculate the sum of only the thermostated atoms
		for  (l=0;l<npartB;++l)
		{				
			if (status[l]==1)
			{
			vxB[l]=(vxB[l]-sumvx); //rescaling velocities
			vyB[l]=(vyB[l]-sumvy);
			vzB[l]=(vzB[l]-sumvz);
			sumkin=sumkin+m_B[l]*vxB[l]*vxB[l]+m_B[l]*vyB[l]*vyB[l]+m_B[l]*vzB[l]*vzB[l];					
			}
			
		}
	} // end of first run
		
	

	if (resume_surf==1){ //resume velocities
		sumkin=0;
		sumkinth=0; //if I want to calculate the sum of only the thermostated atoms
		for  (l=0;l<npartB;++l)
		{					
				sumkin=sumkin+m_B[l]*vxB[l]*vxB[l]+m_B[l]*vyB[l]*vyB[l]+m_B[l]*vzB[l]*vzB[l];
		}
	}


		K[n]=0.5*sumkin;
		
		T[n]=(2*K[n])/(3*kbol*(mobile-1));

	
//-------------------------------------------------We must calculate forces for n=0------------------------------------//
if (tip_version==1){ //dynamic surface
		

	if (mypotential==0){
	
		if (linkcells==1){tersoff_forces_cells(xB,yB,zB);}						
		if (linkcells==0){tersoffopt(xB,yB,zB);}
		U[n]=0.5*sumenergy;
	}
		
	if (mypotential==1){
		
		if (linkcells==1){lj_forces_cells(xB,yB,zB); U[n]=4*0.5*epsilon*sumVcells;} // we must multiply with 0.5 because we use the N*N scheme
		if (linkcells==0){lj_forces(xB,yB,zB); U[n]=4*epsilon*sumV;}   // we do not multiply with 0.5 because we use the N*(N-1)/2 scheme
	}

	
	} 
else if (tip_version==2){ //static surface
	
			U[n]=0;
			for (k=0;k<npartB;++k){sumfx[k]=0;sumfy[k]=0;sumfz[k]=0;}		
	}
		
	
		if (linkcells_supercell==1){lj_ts_forces_cells(xA,xB,yA,yB,zA,zB); Vtotal[n]=4*0.5*epsilon_ts*sumV_ts_cells;}  
		if (linkcells_supercell==0){lj_ts_forces(xA,xB,yA,yB,zA,zB); Vtotal[n]=4*0.5*epsilon_ts*sumV_ts;} 
		
	
		sumforcesxA=0;sumforcesyA=0;sumforceszA=0; //calculation of Ftip
		for (i=0;i<npartA;++i){
			sumforcesxA=sumforcesxA+sumfxA[i];
			sumforcesyA=sumforcesyA+sumfyA[i];
			sumforceszA=sumforceszA+sumfzA[i];
		}		
		Tip_Forcex[n] = sumforcesxA;
		Tip_Forcey[n] = sumforcesyA;
		Tip_Forcez[n] = sumforceszA;
		
		sumforcesxB=0;sumforcesyB=0;sumforceszB=0; //calculation of Fatom
		for (j=0;j<npartB;++j){
			sumforcesxB=sumforcesxB+sumfxB[j];
			sumforcesyB=sumforcesyB+sumfyB[j];
			sumforceszB=sumforceszB+sumfzB[j];
		}		
		
	
		E[n]=K[n]+U[n]; //calculation of initial energy
	
	printf("%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",n,U[n]/npartB,T[n],Vtotal[n]/(npartA+npartB),sumforcesxA,sumforcesyA,sumforceszA);
	TotalForcex = TotalForcex + sumforcesxA; //calculate the average force in x-direction
	TotalForcey = TotalForcey + sumforcesyA; //calculate the average force in y-direction
	TotalForcez = TotalForcez + sumforceszA; //calculate the average force in z-direction
	
//_________________paraview or jmol_____________________________________________________________________//
	
	if (csvout==1){paraview_func(xA,xB,yA,yB,zA,zB);}
	if (jmolout==1){jmol_func(xA,xB,yA,yB,zA,zB);} //export the first .xyz of the file for jmol

	fp=fopen(path_animation,"w+");
	fprintf(fp,"%d\n%s\n",(npartA+npartB),title);	
	for (k=0;k<npartB;++k) //surface
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species_B[k],xB[k],yB[k],zB[k]);
	}	
	for (k=0;k<npartA;++k) //tip
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species_A[k],xA[k],yA[k],zA[k]);
	}
	fclose(fp);

//-----------------------------------forces and topography----------------------------------------------//

	sprintf(path,"%sforces%d.dat",pathout,sim_id);
	fp = fopen(path,"w+");	
	fprintf(fp,"%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",n+nres,Vtotal[n]/(npartA+npartB),sumforcesxA,sumforcesyA,sumforceszA,sumforcesxA,sumforcesyA,sumforceszA);
	fclose(fp);

	sprintf(path,"%stopography%d.dat",pathout,sim_id);
	fp=fopen(path,"w+");
	fprintf(fp,"%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",n+nres,xA[tip_apex],yA[tip_apex],zA[tip_apex],Dseper);
	fclose(fp); 
	
	
//--------calculate the sum of forces from all the particles to the kth particle and accel------------------------//
	

	for (k=0;k<npartB;++k)
	{
		axB[k] = (sumfx[k] + sumfxB[k])/m_B[k];
		ayB[k] = (sumfy[k] + sumfyB[k])/m_B[k];
		azB[k] = (sumfz[k] + sumfzB[k])/m_B[k];
		
		if (status[k]==0)
		{
			axB[k]=0;
			ayB[k]=0;
			azB[k]=0;
		}
	
	}
	
	
//***********************************Loop over time: Velocity Verlet*******************************************************************************************************

	
	for(n=1;n<=N;++n)
	{
		
		// if vx_tip = vy_tip = vz_tip = 0, then we have a tip at a fixed point
			
		switch( tip ){
			
			case 0:   //fixed tip at a certain (x,y,z) point or simple force distance curves (force-spectroscopy)
				for(i=0;i<npartA;++i){
				xfinalA[i] = xA[i] + vx_tip*dt; 
				yfinalA[i] = yA[i] + vy_tip*dt; 
				zfinalA[i] = zA[i] + vz_tip*dt;
				}
				break;
			
			case 1:   //STATIC //constant-height mode
				sumx=0;sumy=0;sumz=0; //calc center of mass
				for (i=0;i<npartA;++i) {
					sumx = sumx + m_A[i]*xA[i];
					sumy = sumy + m_A[i]*yA[i];
					sumz = sumz + m_A[i]*zA[i];}

				for(i=0;i<npartA;++i){
					
						xfinalA[i] = xA[i] + vx_tip*dt + 0.5*dt*dt*(Tip_Forcex[n-1] - Ftargetx)/Mtip;
						yfinalA[i] = yA[i] + vy_tip*dt + 0.5*dt*dt*(Tip_Forcey[n-1] - ky*(sumy/Mtip) + ky*yholder)/Mtip;
						zfinalA[i] = zA[i] + vz_tip*dt + 0.5*dt*dt*(Tip_Forcez[n-1] - kz*(sumz/Mtip) + kz*zholder)/Mtip;

				}
				break;

			case 2:   //STATIC with gamma //constant force mode
				for(i=0;i<npartA;++i){
				xfinalA[i] = xA[i] + vx_tip*dt + 0.5*dt*dt*(Tip_Forcex[n-1] - (gammax*vx_tip))/Mtip;
				yfinalA[i] = yA[i] + vy_tip*dt + 0.5*dt*dt*(Tip_Forcey[n-1] - (gammay*vy_tip))/Mtip;	
				zfinalA[i] = zA[i] + vz_tip*dt + 0.5*dt*dt*(Tip_Forcez[n-1] - Ftargetz - (gammaz*vz_tip))/Mtip; // the tip starts to move in the z direction with a friction term, with an Ftarget
				}
				break;

			case 3:   //STATIC with k in x,y //constant force mode
					sumx=0;sumy=0;sumz=0; //calc center of mass
				for (i=0;i<npartA;++i) {
					sumx = sumx + m_A[i]*xA[i];
					sumy = sumy + m_A[i]*yA[i];
					sumz = sumz + m_A[i]*zA[i];}

				for(i=0;i<npartA;++i){
				xfinalA[i] = xA[i] + vx_tip*dt + 0.5*dt*dt*(Tip_Forcex[n-1] - kx*(sumx/Mtip) + kx*xholder)/Mtip;
				yfinalA[i] = yA[i] + vy_tip*dt + 0.5*dt*dt*(Tip_Forcey[n-1] - ky*(sumy/Mtip) + ky*yholder)/Mtip;
				zfinalA[i] = zA[i] + vz_tip*dt + 0.5*dt*dt*(Tip_Forcez[n-1] - Ftargetz - (gammaz*vz_tip))/Mtip;
				}
				break;

			case 4:   //STATIC with k in x,y and gx,gy //constant force mode (FULL MODEL)
				sumx=0;sumy=0;sumz=0; //calc center of mass
				for (i=0;i<npartA;++i) {
					sumx = sumx + m_A[i]*xA[i];
					sumy = sumy + m_A[i]*yA[i];
					sumz = sumz + m_A[i]*zA[i];}

				// we need to export the full force that the system feels
					Full_Fx = (Tip_Forcex[n-1] - Ftargetx - kx*(sumx/Mtip) + kx*xholder - (gammax*vx_tip));
					Full_Fy = (Tip_Forcey[n-1] - Ftargety - ky*(sumy/Mtip) + ky*yholder - (gammay*vy_tip));
					Full_Fz = (Tip_Forcez[n-1] - Ftargetz - (gammaz*vz_tip));


				for(i=0;i<npartA;++i){
				xfinalA[i] = xA[i] + vx_tip*dt + (0.5*dt*dt*Full_Fx)/Mtip;
				yfinalA[i] = yA[i] + vy_tip*dt + (0.5*dt*dt*Full_Fy)/Mtip;
				zfinalA[i] = zA[i] + vz_tip*dt + (0.5*dt*dt*Full_Fz)/Mtip;
				
				//if (xfinalA[l]>xmax){xfinalA[l]=xfinalA[l]-Lx;}
				//if (xfinalA[l]<xmin){xfinalA[l]=xfinalA[l]+Lx;}
				//if (yfinalA[l]>ymax){yfinalA[l]=yfinalA[l]-Ly;}
				//if (yfinalA[l]<ymin){yfinalA[l]=yfinalA[l]+Ly;}
		
				//if (zfinalA[l]>zmax){zfinalA[l]=zfinalA[l]-Lz;}
				//if (zfinalA[l]<zmin){zfinalA[l]=zfinalA[l]+Lz;}
				}	
				break;

			default: 
				{printf("You have selected a case that does not exist!!!!!!! Try again!!!!!!"); exit(1);}

			}
			

		
		for(l=0;l<npartB;++l) //surface
		{
			xfinalB[l]=xB[l]+vxB[l]*dt+0.5*dt*dt*axB[l];
			yfinalB[l]=yB[l]+vyB[l]*dt+0.5*dt*dt*ayB[l];
			zfinalB[l]=zB[l]+vzB[l]*dt+0.5*dt*dt*azB[l];
			
			if (xfinalB[l]>xmax){xfinalB[l]=xfinalB[l]-Lx;}
			if (xfinalB[l]<xmin){xfinalB[l]=xfinalB[l]+Lx;}
			if (yfinalB[l]>ymax){yfinalB[l]=yfinalB[l]-Ly;}
			if (yfinalB[l]<ymin){yfinalB[l]=yfinalB[l]+Ly;}
		
			if (zfinalB[l]>zmax){zfinalB[l]=zfinalB[l]-Lz;}
			if (zfinalB[l]<zmin){zfinalB[l]=zfinalB[l]+Lz;}
				
		}

 if (tip_version==1){ //dynamic surface

		
	if (mypotential==0){
	
		if (linkcells==1){tersoff_forces_cells(xfinalB,yfinalB,zfinalB);}						
		if (linkcells==0){tersoffopt(xfinalB,yfinalB,zfinalB);}
		U[n]=0.5*sumenergy;
	}
		
	if (mypotential==1){
		
		if (linkcells==1){lj_forces_cells(xfinalB,yfinalB,zfinalB); U[n]=4*0.5*epsilon*sumVcells;} // we must multiply with 0.5 because we use the N*N scheme
		if (linkcells==0){lj_forces(xfinalB,yfinalB,zfinalB); U[n]=4*epsilon*sumV;}   // we do not multiply with 0.5 because we use the N*(N-1)/2 scheme
	}
 }
	else if (tip_version==2){ //static version
	
			U[n]=0;	
	}	



		if (linkcells_supercell==1){lj_ts_forces_cells(xfinalA,xfinalB,yfinalA,yfinalB,zfinalA,zfinalB); Vtotal[n]=4*0.5*epsilon_ts*sumV_ts_cells;}  
		if (linkcells_supercell==0){lj_ts_forces(xfinalA,xfinalB,yfinalA,yfinalB,zfinalA,zfinalB); Vtotal[n]=4*0.5*epsilon_ts*sumV_ts;} 
		
		sumforcesxA=0;sumforcesyA=0;sumforceszA=0; //calculation of Ftip
		for (i=0;i<npartA;++i){
			sumforcesxA=sumforcesxA+sumfxA[i];
			sumforcesyA=sumforcesyA+sumfyA[i];
			sumforceszA=sumforceszA+sumfzA[i];
		}		
		Tip_Forcex[n] = sumforcesxA;
		Tip_Forcey[n] = sumforcesyA;
		Tip_Forcez[n] = sumforceszA;
		
		sumforcesxB=0;sumforcesyB=0;sumforceszB=0; //calculation of Fatom
		for (j=0;j<npartB;++j){
			sumforcesxB=sumforcesxB+sumfxB[j];
			sumforcesyB=sumforcesyB+sumfyB[j];
			sumforceszB=sumforceszB+sumfzB[j];
		}		

		
switch( tip ){	//velocity of the center of mass

			case 0:   //fixed tip at a certain (x,y,z) point or force distance curves (force-spectroscopy)
					break;
			case 1:   //STATIC //constant-height mode
					sumxfinal=0;sumyfinal=0;sumzfinal=0; //calc center of mass
					for (i=0;i<npartA;++i) {
					sumxfinal = sumxfinal + m_A[i]*xfinalA[i];
					sumyfinal = sumyfinal + m_A[i]*yfinalA[i];
					sumzfinal = sumzfinal + m_A[i]*zfinalA[i];}

					
						if (xholder>xmax||xholder<xmin||yholder>ymax||yholder<ymin){ //if tip exceeds box do not move
							vx_tip=0;vy_tip=0;vz_tip=0;
						} 

						else{
						vx_tip = vx_tip + dt*(Tip_Forcex[n]+Tip_Forcex[n-1] - 2*Ftargetx)/(2*Mtip);
						vy_tip = vy_tip + dt*(Tip_Forcey[n]+Tip_Forcey[n-1] - ky*(sumy/Mtip) - ky*(sumyfinal/Mtip) + 2*ky*yholder)/(2*Mtip);
						vz_tip = vz_tip + dt*(Tip_Forcez[n]+Tip_Forcez[n-1] - kz*(sumz/Mtip) - kz*(sumzfinal/Mtip) + 2*kz*zholder)/(2*Mtip);
						}
					break;
			case 2:   //STATIC with gamma //constant force mode
					
					vx_tip = ( (1 - (gammax*dt)/(2*Mtip) )*vx_tip + dt*(Tip_Forcex[n]+Tip_Forcex[n-1])/(2*Mtip)) / ( 1 + (gammax*dt)/(2*Mtip) );
					vy_tip = ( (1 - (gammay*dt)/(2*Mtip) )*vy_tip + dt*(Tip_Forcey[n]+Tip_Forcey[n-1])/(2*Mtip)) / ( 1 + (gammay*dt)/(2*Mtip) );
					vz_tip = ( (1 - (gammaz*dt)/(2*Mtip) )*vz_tip + dt*(Tip_Forcez[n]+Tip_Forcez[n-1]-2*Ftargetz)/(2*Mtip)) / ( 1 + (gammaz*dt)/(2*Mtip) ); // calc velocity with friction term
					
					break;
			case 3:   //STATIC with k in x,y //constant force mode
					sumxfinal=0;sumyfinal=0;sumzfinal=0; //calc center of mass
					for (i=0;i<npartA;++i) {
					sumxfinal = sumxfinal + m_A[i]*xfinalA[i];
					sumyfinal = sumyfinal + m_A[i]*yfinalA[i];
					sumzfinal = sumzfinal + m_A[i]*zfinalA[i];}

					vx_tip = vx_tip + dt*(Tip_Forcex[n]+Tip_Forcex[n-1] - kx*(sumx/Mtip) - kx*(sumxfinal/Mtip) + 2*kx*xholder)/(2*Mtip);
					vy_tip = vy_tip + dt*(Tip_Forcey[n]+Tip_Forcey[n-1] - ky*(sumy/Mtip) - ky*(sumyfinal/Mtip) + 2*ky*yholder)/(2*Mtip);
					vz_tip = ( (1 - (gammaz*dt)/(2*Mtip) )*vz_tip + dt*(Tip_Forcez[n]+Tip_Forcez[n-1]-2*Ftargetz)/(2*Mtip)) / ( 1 + (gammaz*dt)/(2*Mtip) ); // calc velocity with friction term
					
					break;
			case 4:   //STATIC with k in x,y and gx,gy //constant force mode (FULL MODEL)
					sumxfinal=0;sumyfinal=0;sumzfinal=0; //calc center of mass
					for (i=0;i<npartA;++i) {
					sumxfinal = sumxfinal + m_A[i]*xfinalA[i];
					sumyfinal = sumyfinal + m_A[i]*yfinalA[i];
					sumzfinal = sumzfinal + m_A[i]*zfinalA[i];}

					vx_tip = ( (1 - (gammax*dt)/(2*Mtip) )*vx_tip + dt*(Tip_Forcex[n]+Tip_Forcex[n-1]- 2*Ftargetx - kx*(sumx/Mtip) - kx*(sumxfinal/Mtip) + 2*kx*xholder)/(2*Mtip)) / ( 1 + (gammax*dt)/(2*Mtip) );
					vy_tip = ( (1 - (gammay*dt)/(2*Mtip) )*vy_tip + dt*(Tip_Forcey[n]+Tip_Forcey[n-1]- 2*Ftargety - ky*(sumy/Mtip) - ky*(sumyfinal/Mtip) + 2*ky*yholder)/(2*Mtip)) / ( 1 + (gammay*dt)/(2*Mtip) );
					vz_tip = ( (1 - (gammaz*dt)/(2*Mtip) )*vz_tip + dt*(Tip_Forcez[n]+Tip_Forcez[n-1]- 2*Ftargetz)/(2*Mtip)) / ( 1 + (gammaz*dt)/(2*Mtip) ); // calc velocity with friction term
					
					break;

			default: 
					{printf("You have selected a case that does not exist!!!!!!! Try again!!!!!!"); exit(1);}	
	}
		
		
		for (k=0;k<npartB;++k)
		{
			axfinalB[k] = (sumfx[k]+sumfxB[k])/m_B[k];
			ayfinalB[k] = (sumfy[k]+sumfyB[k])/m_B[k];
			azfinalB[k] = (sumfz[k]+sumfzB[k])/m_B[k];
			
			if (status[k]==0)
			{
				axfinalB[k]=0;
				ayfinalB[k]=0;
				azfinalB[k]=0;
			}
					
			vxfinalB[k]=vxB[k]+0.5*axB[k]*dt+0.5*axfinalB[k]*dt;
			vyfinalB[k]=vyB[k]+0.5*ayB[k]*dt+0.5*ayfinalB[k]*dt;
			vzfinalB[k]=vzB[k]+0.5*azB[k]*dt+0.5*azfinalB[k]*dt;			
			
			
		} // end of loop over k 
		
//************************************* Untill here we have NVE**************************************************************
		
//-----------------------------------------ANDERSEN or Berendsen THERMOSTAT-------------------------------------------------------------------------------//	
	if (tip_version==1){ //dynamic surface
		if 	(choose_thermostat==1){
			
		Andersen(vxfinalB,vyfinalB,vzfinalB);}

		else if (choose_thermostat==2){
		
		vrescale(vxfinalB,vyfinalB,vzfinalB);}
		
		 sumkin=0;
		 for (k=0;k<npartB;++k)
		 {sumkin=sumkin+m_B[k]*vxfinalB[k]*vxfinalB[k]+m_B[k]*vyfinalB[k]*vyfinalB[k]+m_B[k]*vzfinalB[k]*vzfinalB[k];}
		 
		 
		K[n]=0.5*sumkin;
		 
		
		}

	else if (tip_version==2){ //static surface
			K[n]=0;
			for (k=0;k<npartB;++k){vxfinalB[k]=0;vyfinalB[k]=0;vzfinalB[k]=0;axfinalB[k]=0;ayfinalB[k]=0;azfinalB[k]=0;}
		}

		T[n]=(2*K[n])/(3*(mobile-1)*kbol);	
		E[n]=K[n]+U[n];
		
		if (n%nstep==0){printf("%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",n,U[n]/npartB,T[n],Vtotal[n]/(npartA+npartB),sumforcesxA,sumforcesyA,sumforceszA);}
		TotalForcex = TotalForcex + sumforcesxA; //calculate the average force in x-direction
		TotalForcey = TotalForcey + sumforcesyA; //calculate the average force in y-direction
		TotalForcez = TotalForcez + sumforceszA; //calculate the average force in z-direction
		
//----------------------------------------------end of THERMOSTATS-----------------------------------------------------------------------------------------		
		
		
//---------------------------------start the ACF calculation for Total DOS----------------------------------------------------------------------------------//

		acf(vxfinalB,vyfinalB,vzfinalB);
	
				
//****************************  update **************************************************
		
		// update tip positions and velocity of the center of mass 
		for(i=0;i<npartA;++i)
		{
			xA[i]=xfinalA[i];
			yA[i]=yfinalA[i];
			zA[i]=zfinalA[i];
		}
		
		// update surface positions, velocities and accelerations
		for(l=0;l<npartB;++l)
		{
			xB[l]=xfinalB[l];
			yB[l]=yfinalB[l];
			zB[l]=zfinalB[l];
			
			vxB[l]=vxfinalB[l];
			vyB[l]=vyfinalB[l];
			vzB[l]=vzfinalB[l];
			
			axB[l]=axfinalB[l];
			ayB[l]=ayfinalB[l];
			azB[l]=azfinalB[l];
			
		}
		
		
		
		
//_________________paraview or jmol________________________________________________________//

	if (n%graph_step==0) {	

		if (csvout==1){paraview_func(xA,xB,yA,yB,zA,zB);}
		
		//-----------jmol animation------------------------//
		fp=fopen(path_animation,"a");
		fprintf(fp,"%d\n%s\n",(npartA+npartB),title);
		
		for (k=0;k<npartB;++k) //surface
		{
			fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species_B[k],xB[k],yB[k],zB[k]);
		}
		
		for (k=0;k<npartA;++k) //tip
		{
			fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species_A[k],xA[k],yA[k],zA[k]);
		}
		
		fclose(fp);
		//------------------------------------------------//
	}

//-------------------------export forces-------------------------------------//
	sprintf(path,"%sforces%d.dat",pathout,sim_id);
	fp = fopen(path,"a");		
	fprintf(fp,"%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",n+nres,Vtotal[n]/(npartA+npartB),sumforcesxA,sumforcesyA,sumforceszA,Full_Fx,Full_Fy,Full_Fz);
	fclose(fp);

		
//--------------------- export the tip positions for topography-----------------------//		
		
			Dseper = zA[tip_apex] - zB[npartB-1]; // tip - surface seperation if surface flat static
			sprintf(path,"%stopography%d.dat",pathout,sim_id);
			fp=fopen(path,"a");
			fprintf(fp,"%d\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",n+nres,xA[tip_apex],yA[tip_apex],zA[tip_apex],Dseper);
			fclose(fp);
			
//------export last positions of surface for certain intervals----------------------------------//

if (n%intervals==0){	
	sprintf(path,"%sresume.dat",pathout); 
	fp = fopen(path,"w+");
	fprintf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",titleB,npartB,xmin,xmax,ymin,ymax,zmin,zmax);
	fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");	

	for(k=0;k<npartB;++k)
	{
		fprintf(fp,"%s\t%.16lf\t%.16lf\t%.16lf\t%d\t%d\t%.16lf\t%.16lf\t%.16lf\n",species_B[k],xfinalB[k],yfinalB[k],zfinalB[k],status[k],thermostatus[k],vxfinalB[k],vyfinalB[k],vzfinalB[k]);
	}
	fprintf(fp,"%d\n",n-1+nres);
	fclose(fp);

//------------------------------export last positions of tip for certain intervals----------------//
	
	sprintf(path,"%sresume_tip.dat",pathout);
	fp = fopen(path,"w+");
	fprintf(fp,"%d\n%s\n",npartA,titleA);
	fprintf(fp,"%lf\t%lf\t%lf\n",vx_tip,vy_tip,vz_tip);
	fprintf(fp,"%lf\t%lf\t%lf\n",xholder,yholder,zholder);
	
	for (k=0;k<npartA;++k)
	{		
		fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species_A[k],xfinalA[k],yfinalA[k],zfinalA[k]);		
	}

	fclose(fp);			
}			 
		
//------------------------monitor various quantities-----------------------------------------------------------------------------------------
		
		
		if (monitor==1){
			
			sprintf(path,"%sfi%d.dat",pathf,n);
			
			sprintf(path,"%srdf%d.dat",pathrdf,n);
			fp=fopen(path,"w+");
			
			for (p=0;p<bins;++p)			
			{
				rinst[p] = (p+1-0.5)*dr;
				ginst[p] = ((Lx*Ly*Lz*hinst[p])/(4*pi*rinst[p]*rinst[p]*dr*npartB*npartB)); // we calc rdf for each time step
				
				fprintf(fp,"%.16lf\t%.16lf\n",rinst[p],ginst[p]);
				hinst[p]=0;				
			}
			
			fclose(fp);
			
		} //end monitor		
		
	
		
		
}  //*********************************end time of simulation*****************************************************************************

	if (jmolout==1){jmol_func(xA,xB,yA,yB,zA,zB);} //export the last .xyz of the file for jmol


		//---------------Print averages------------------------------------//
		AvgForcex = (TotalForcex/N);
		AvgForcey = (TotalForcey/N);
		AvgForcez = (TotalForcez/N);

		printf("\n\nThe Force averages are:\t  ");
		printf("%lf\t%lf\t%lf\n",AvgForcex,AvgForcey,AvgForcez);
	
	time(&rawtime);
	sprintf (stoptime, "\nSimulation ended at: %s",ctime(&rawtime));
	
	printf ("%s\n", starttime);
	printf ("%s\n", stoptime);
	
// call RDF
	
	RDF();
//------------------------EXPORT VARIOUS DATA---------------------------------------------------------------------------------------------------------------	
	
	export_info();

//_______________________free vectors from memory_________________________________________________________________________________________//
	
	free(m_A);free(Z_A);free(covalent_r_A);free(vdW_r_A);
	free(m_B);free(Z_B);free(covalent_r_B);free(vdW_r_B);		
	free(status); free(thermostatus);
	for (i=0; i<npartA; ++i){free(species_A[i]);}	
	free(species_A);
	for (i=0; i<npartB; ++i){free(species_B[i]);}	
	free(species_B);

	free(Tinit);free(Tfinal);free(thermostat);
	
	free(myA);free(myB);free(lamda);free(mi);free(beta);
	free(ni);free(ci);free(di);free(hi);free(myR);free(myS);
	free(myepsilon);free(mysigma);
	
	free(xA); free(yA); free(zA); //fixed tip
	free(xfinalA); free(yfinalA); free(zfinalA); //moving tip
	free(xB); free(yB); free(zB);
	free(xfinalB); free(yfinalB); free(zfinalB);
	
	free(vxB); free(vyB); free(vzB);
	free(vxfinalB); free(vyfinalB); free(vzfinalB);
	
	free(axB); free(ayB); free(azB);
	free(axfinalB); free(ayfinalB); free(azfinalB);
		
	free(U); free(K); free(E); free(T);
	free(Vtotal);
	free(Rep); free(Attr); free(Bond);
	free(dx_zeta_III);
	free(dy_zeta_III);
	free(dz_zeta_III);
	
	free(vxfixed); free(vyfixed); free(vzfixed); free(acf_t); free(sumtime);
	free(h); free(g); free(r1); free(hinst); free(ginst); free(rinst);
	
	free(sumfx);free(sumfy);free(sumfz);	
	
	free(sumfxA);free(sumfyA);free(sumfzA);
	free(sumfxB);free(sumfyB);free(sumfzB);
	free(myepsilonA);free(myepsilonB);free(mysigmaA);free(mysigmaB);
	
	free(Tip_Forcex);free(Tip_Forcey);free(Tip_Forcez);
	

	if (linkcells==1)
	{	
		if (mypotential==0){
		free(neighborcellpartarray);
		free(allpartarray);
		free(cellpartarray);}
		
		if (mypotential==1){		
		free(partindex);
		free(nextpartindex);}	
	
		for (i = 0; i < cellsy+2; ++i){free(cells[i]);}
		free(cells);
		
		for (i = 0; i < 26; ++i){free(neighbours[i]);}
		free(neighbours);
		
		free(HEAD);free(LIST);
	}
	
	if (linkcells_supercell==1)
	{
		free(HEAD_A);free(HEAD_B);free(LIST_A);free(LIST_B);
		free(partindex_A);free(partindex_B);free(nextpartindex_B);
		
		for (i = 0; i < cellsy_supercell+2; ++i){free(cells_supercell[i]);}
		free(cells_supercell);
		
		for (i = 0; i < 26; ++i){free(neighbours_supercell[i]);}
		free(neighbours_supercell);
	}

	if(OS==1){_CrtDumpMemoryLeaks();}
	
	getchar();getchar();
	return 0;
	
	
}	// end main
