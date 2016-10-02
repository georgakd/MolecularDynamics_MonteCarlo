

//********************************************* MOLDYN PROJECT *************************************************//
//		RELEASE MODE v.2 17/12/2012																				//
//		Created by Dimitra Georgakaki																			//																					
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

void linked_cells(int cellsx, int cellsy, int cellsz); 
void Andersen(double *vx, double *vy, double *vz); // declare the function that initiates the Andersen thermostat
void vrescale(double *vx, double *vy, double *vz); // declare the function that initiates the Berendsen thermostat
void RDF(void); // Radial Distribution Function
void acf(double *vx,double *vy,double *vz); // ACF

//-------- declare the function that will calculate interaction forces and potential
 
void tersoff_forces_cells(double *x,double *y,double *z);
void tersoffopt(double *x,double *y,double *z);
void lj_forces(double *x,double *y,double *z); 
void lj_forces_cells(double *x,double *y,double *z);
//----------------------------------------------------------------------------------

void jmol_func(double *x,double *y, double *z);
void paraview_func(double *x,double *y, double *z);
void export_info(void);



int main(void)
{
	
//*************************************************** Define local vars ******************************************************************************//
	

	double sumvx,sumvy,sumvz; // local sums
	double sumU; //initialize for avg[U]
	double sumKin; //initialize for avg[Kin]
	
	char pathf[cmax_length]; // for monitor instant forces in matlab
	char pathrdf[cmax_length]; // for monitor instant rdf in matlab	
	char mycommand[cmax_length]; // for mkdir and delete files and folders
	char path_animation[cmax_length]; //for jmol animation
	
	char starttime[cmax_length],stoptime[cmax_length];	// for counting the simulation

	
	time_t rawtime;
	time(&rawtime);
	sprintf (starttime, "\nSimulation started at: %s",ctime(&rawtime));
	

//--------------------------------------------Read external parameters--------------------------------------//
	if (OS==1)
		{_getcwd(sim_folder,cmax_length);}
	else
		{getcwd(sim_folder,cmax_length);}

//--------------------------------------------------------------Creating folders, paths--------------------------------------------------------------//
	
	if (OS==1)
	{
	
	sprintf(file_path,"%s\\config.txt",sim_folder);
	sprintf(pathcsv,"%s\\csv\\",sim_folder);
	sprintf(pathout,"%s\\out\\",sim_folder);
	sprintf(pathf,"%s\\forces\\",sim_folder);
	sprintf(pathrdf,"%s\\radial\\",sim_folder);
	sprintf(pathxyz,"%s\\jmol\\",sim_folder);
	sprintf(path_animation,"%s\\out\\animation.xyz",sim_folder);
	}
	else 
	{	
	sprintf(file_path,"%s/config.txt",sim_folder);
	sprintf(pathcsv,"%s/csv/",sim_folder);
	sprintf(pathout,"%s/out/",sim_folder);
	sprintf(pathf,"%s/forces/",sim_folder);
	sprintf(pathrdf,"%s/radial/",sim_folder);
	sprintf(pathxyz,"%s/jmol/",sim_folder);
	sprintf(path_animation,"%s/out/animation.xyz",sim_folder);
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
	
	printf("START by defining the simulation folder: ");printf("%s",sim_folder);printf("\n");

	parameters();  // call the function parameters.c
	
//----------------------------------------------Delete previous folders csv,jmol,radial,forces if they exist----------------------------------------//
	
	printf("\nDeleting all previous existing files in the selected folders!! \n ");
	
	if (OS==1){
		sprintf(mycommand,"del /Q %s*.csv",pathcsv);system(mycommand);
		sprintf(mycommand,"del /Q %s*.xyz",pathxyz);system(mycommand);
	}
	else{
		sprintf(mycommand,"find %s -name '*.csv' | xargs rm",pathcsv);system(mycommand);
		sprintf(mycommand,"find %s -name '*.xyz' | xargs rm",pathxyz);system(mycommand);
	}



//----------------------Structure information--------------------------//

if (resume_flag==0){ 	// import for the first time initial positions
	if (OS==1){sprintf(path,"%s\\positions.dat",sim_folder);}
	else {sprintf(path,"%s/positions.dat",sim_folder);} 
	
	fp=fopen(path,"r");
	if(fp==NULL){printf("could not locate %s\n",path);exit(-1);}
	fscanf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",title,&npart,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
	fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\n");	

	//preallocations for surface		
	species=(char**)malloc((npart)*sizeof(char*));		
	for (k = 0; k < npart; ++k)
	{species[k]=(char*)malloc((3)*sizeof(char));}
	x1=(double*)malloc((npart)*sizeof(double));
	y_1=(double*)malloc((npart)*sizeof(double));
	z1=(double*)malloc((npart)*sizeof(double));
	status=(int*)malloc((npart)*sizeof(int));
	thermostatus=(int*)malloc((npart)*sizeof(int));
	vx1=(double*)malloc((npart)*sizeof(double));
	vy_1=(double*)malloc((npart)*sizeof(double));
	vz1=(double*)malloc((npart)*sizeof(double));
		
	for (k=0;k<npart;++k){		
		fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[k],&x1[k],&y_1[k],&z1[k],&status[k],&thermostatus[k]);}
	fclose(fp);
}
	
	
else if (resume_flag==1){ // resume simulation 

	if (OS==1){sprintf(path,"%s\\resume.dat",pathout);}
	else {sprintf(path,"%s/resume.dat",pathout);}  
	fp=fopen(path,"r");
	if(fp==NULL){printf("could not locate %s\n",path);exit(-1);}
	fscanf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",title,&npart,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
	fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");

	//preallocations for surface		
	species=(char**)malloc((npart)*sizeof(char*));		
	for (k = 0; k < npart; ++k)
	{species[k]=(char*)malloc((3)*sizeof(char));}
	x1=(double*)malloc((npart)*sizeof(double));
	y_1=(double*)malloc((npart)*sizeof(double));
	z1=(double*)malloc((npart)*sizeof(double));
	status=(int*)malloc((npart)*sizeof(int));
	thermostatus=(int*)malloc((npart)*sizeof(int));
	vx1=(double*)malloc((npart)*sizeof(double));
	vy_1=(double*)malloc((npart)*sizeof(double));
	vz1=(double*)malloc((npart)*sizeof(double));


	for(k=0;k<npart;++k){
		fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",species[k],&x1[k],&y_1[k],&z1[k],&status[k],&thermostatus[k],&vx1[k],&vy_1[k],&vz1[k]);}
	fclose(fp);
}
	


// ------------------------------------------------------- Allocation of all other arrays -----------------------------------------------------------
// Species information		
	Z=(int*)malloc((npart)*sizeof(int));
	m=(double*)malloc((npart)*sizeof(double));
	covalent_r=(double*)malloc((npart)*sizeof(double));
	vdW_r=(double*)malloc((npart)*sizeof(double));	

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


//positions, velocities, accelerations after Verlet
	x2=(double*)malloc((npart)*sizeof(double));
	y2=(double*)malloc((npart)*sizeof(double));
	z2=(double*)malloc((npart)*sizeof(double));
	
	vx2=(double*)malloc((npart)*sizeof(double));
	vy2=(double*)malloc((npart)*sizeof(double));
	vz2=(double*)malloc((npart)*sizeof(double));
	
	ax1=(double*)malloc((npart)*sizeof(double));
	ay_1=(double*)malloc((npart)*sizeof(double));
	az1=(double*)malloc((npart)*sizeof(double));
	ax2=(double*)malloc((npart)*sizeof(double));
	ay2=(double*)malloc((npart)*sizeof(double));
	az2=(double*)malloc((npart)*sizeof(double));	
	
//velocity acf
	
	vxfixed=(double*)malloc((npart)*sizeof(double));
	vyfixed=(double*)malloc((npart)*sizeof(double));
	vzfixed=(double*)malloc((npart)*sizeof(double));
	acf_t=(double*)malloc((npart)*sizeof(double)); // particles sum
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
	
	
//Be careful, forces are two-dim arrays but they can be transformed to 1D!
	
	sumfx=(double*)malloc((npart)*sizeof(double)); //forces 1D (for speed)
	sumfy=(double*)malloc((npart)*sizeof(double));
	sumfz=(double*)malloc((npart)*sizeof(double));	


//---------------------------------------------Define secondary parameters-------------------------------------------------//
	
//***********be very careful this is the correct order for the parameters to be read from the program*******************//
		
	
	if (mypotential==0){printf("\n  Tersoff T3 potential...\n");}
	if (mypotential==1){printf("\n  Lennard Jones potential...\n");}

	mobile=0;thermo=0;
	for (k=0;k<npart;++k){mobile=mobile+status[k];} //count the mobile atoms
	for (k=0;k<npart;++k){thermo=thermo+thermostatus[k];} //count the thermostated atoms

	Lx=xmax-xmin; // box simulation
	Ly=ymax-ymin; 
	Lz=zmax-zmin;
	
	
	dt=dtSI/dtscale;  // the reduced time step used in the simulation

	TempStep = (Ttarget-Tcool)/Tintervals; //for the simulated anealing/quenching (if anealing do not forget that Tcool=0)

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
	
//--------------------------------------------------- Use or not use LINKED CELLS ------------------------------------------------------------------
	
	if (linkcells==1){
		// 1. Divide space for Tersoff interactions
		printf("\n\n--------------------------LINKED CELLS----------------------------------- \n");
		printf("Dividing space for sample.......\n");
		printf("\n\nGive me the maximum rcut for the given system: rcut_max = "); // for a multicomponent system i choose the max S
		scanf("%lf",&rcmax);
		
		Lcx=Lx/cellsx;Lcy=Ly/cellsy;Lcz=Lz/cellsz;
		
		if (Lcx>=rcmax)
			printf("The division of cells in x direction is: %d  %lf %lf\n",cellsx, Lcx, rcmax);
		else 
			printf("Sorry, try again!\n");
		if (Lcy>=rcmax)
			printf("The division of cells in y direction is: %d %lf %lf\n",cellsy, Lcy, rcmax);
		else 
			printf("Sorry, try again!\n"); 
		if (Lcz>=rcmax)
			printf("The division of cells in z direction is: %d %lf %lf\n\n",cellsz, Lcz, rcmax);
		else 
			printf("Sorry, try again!\n\n");
		
		
	linked_cells(cellsx,cellsy,cellsz);
		
	printf("------------------------------------------------------------------------- \n");	
	}

//***************************************All initializations done, start the program!!!*********************************************************************	

			
	printf("All initializations done, start the program!!!\n");
	printf("\nTotal atoms: %d\n", npart);printf("The mobile atoms are: %d\n", mobile);printf("The thermostated atoms are: %d\n\n", thermo);

// Initial Simulation n=0: Rescaling of velocities drawn from a gaussian distribution
	n=0; sumU=0; sumKin=0;
	
	if (resume_flag==0){ // the first run is always with random velocities


		sumvx=0;sumvy=0;sumvz=0;		
		for (l=0;l<npart;++l)
		{
			if (status[l]==0){vx1[l]=0;vy_1[l]=0;vz1[l]=0;}
			
			else if (status[l]==1)
			{
			
			myrandom1=random_gen(); //my own Witchmann-Hill generator
			myrandom2=random_gen();
				
			vx1[l]=sqrt((kbol*Ttarget)/m[l])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2); // Box-Muller Method
			vy_1[l]=sqrt((kbol*Ttarget)/m[l])*sqrt(-2*log(myrandom1))*sin(2*pi*myrandom2);
			
			myrandom1=random_gen();
			myrandom2=random_gen();
			
			vz1[l]=sqrt((kbol*Ttarget)/m[l])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2);		
			
			sumvx=sumvx+vx1[l];sumvy=sumvy+vy_1[l];sumvz=sumvz+vz1[l];  //velocity center of mass for mobile particles
			
			}
			
		}	

		
		sumvx=sumvx/mobile; sumvy=sumvy/mobile; sumvz=sumvz/mobile;  //   <v>

		sumkin=0;sumkinth=0;
		for  (l=0;l<npart;++l)
		{	
			
			if (status[l]==1)
			{
			vx1[l]=(vx1[l]-sumvx); //rescaling velocities
			vy_1[l]=(vy_1[l]-sumvy);
			vz1[l]=(vz1[l]-sumvz);
			sumkin=sumkin+m[l]*vx1[l]*vx1[l]+m[l]*vy_1[l]*vy_1[l]+m[l]*vz1[l]*vz1[l];
					
			}
			
		}
	} // end of first run
		
	

	if (resume_flag==1){ //resume velocities
		sumkin=0;sumkinth=0;
		for  (l=0;l<npart;++l)
		{	
				
			sumkin=sumkin+m[l]*vx1[l]*vx1[l]+m[l]*vy_1[l]*vy_1[l]+m[l]*vz1[l]*vz1[l];
			
		}
	}


		K[n]=0.5*sumkin;
		
		T[n]=(2*K[n])/(3*kbol*(mobile-1));

	
//-----------------------------We must calculate forces for n=0------------------------------//
	if (mypotential==0){
	
		if (linkcells==1){tersoff_forces_cells(x1,y_1,z1);}						
		if (linkcells==0){tersoffopt(x1,y_1,z1);}
		U[n]=0.5*sumenergy;
	}
		
	if (mypotential==1){
		
		if (linkcells==1){lj_forces_cells(x1,y_1,z1); U[n]=4*0.5*epsilon*sumVcells;} // we must multiply with 0.5 because we use the N*N scheme
		if (linkcells==0){lj_forces(x1,y_1,z1); U[n]=4*epsilon*sumV;}   // we do not multiply with 0.5 because we use the N*(N-1)/2 scheme
	}


		E[n]=K[n]+U[n];
	
	printf("%d\t%.16lf\t%.16lf\n",n,U[0]/npart,T[0]);
	
	
	sumU = sumU + U[n]; //calculate the average
	sumKin = sumKin + K[n]; //calculate the average	
	
//_________________paraview or jmol for n=0;_____________________________________________________________________
	
	if (csvout==1){paraview_func(x1,y_1,z1);}
	if (jmolout==1){jmol_func(x1,y_1,z1);} //export the first .xyz file for jmol

	fp=fopen(path_animation,"w+");
	fprintf(fp,"%d\n%s\n",npart,title);
	for (k=0;k<npart;++k) // n=0
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species[k],x1[k],y_1[k],z1[k]);
	}
	fclose(fp);
	
	
//--------calculate the sum of forces from all the particles to the ith particle and accel------------------------
	

	for (k=0;k<npart;++k)
	{
		ax1[k]=sumfx[k]/m[k];
		ay_1[k]=sumfy[k]/m[k];
		az1[k]=sumfz[k]/m[k];
		if (status[k]==0)
		{
			ax1[k]=0;
			ay_1[k]=0;
			az1[k]=0;
		}
	
	}
	
	
//***********************************Loop over time*******************************************************************************************************
	
	
	for(n=1;n<=N;++n)
	{
		for(l=0;l<npart;++l)
		{
			x2[l]=x1[l]+vx1[l]*dt+0.5*dt*dt*ax1[l];
			y2[l]=y_1[l]+vy_1[l]*dt+0.5*dt*dt*ay_1[l];
			z2[l]=z1[l]+vz1[l]*dt+0.5*dt*dt*az1[l];
			
			if (x2[l]>xmax){x2[l]=x2[l]-Lx;}
			if (x2[l]<xmin){x2[l]=x2[l]+Lx;}
			if (y2[l]>ymax){y2[l]=y2[l]-Ly;}
			if (y2[l]<ymin){y2[l]=y2[l]+Ly;}
		
			if (z2[l]>zmax){z2[l]=z2[l]-Lz;}
			if (z2[l]<zmin){z2[l]=z2[l]+Lz;}
				
		}
		
	
	if (mypotential==0){
	
		if (linkcells==1){tersoff_forces_cells(x2,y2,z2);}						
		if (linkcells==0){tersoffopt(x2,y2,z2);}
		U[n]=0.5*sumenergy;
	}
		
	if (mypotential==1){
		
		if (linkcells==1){lj_forces_cells(x2,y2,z2); U[n]=4*0.5*epsilon*sumVcells;} // we must multiply with 0.5 because we use the N*N scheme
		if (linkcells==0){lj_forces(x2,y2,z2); U[n]=4*epsilon*sumV;}   // we do not multiply with 0.5 because we use the N*(N-1)/2 scheme
	}

		
		for (k=0;k<npart;++k)
		{
			ax2[k]=sumfx[k]/m[k];
			ay2[k]=sumfy[k]/m[k];
			az2[k]=sumfz[k]/m[k];
			if (status[k]==0)
			{
				ax2[k]=0;
				ay2[k]=0;
				az2[k]=0;
			}
		
			
			vx2[k]=vx1[k]+0.5*ax1[k]*dt+0.5*ax2[k]*dt;
			vy2[k]=vy_1[k]+0.5*ay_1[k]*dt+0.5*ay2[k]*dt;
			vz2[k]=vz1[k]+0.5*az1[k]*dt+0.5*az2[k]*dt;			
			
			
		} // end of loop over k 
		
//************************************* Untill here we have NVE**************************************************************
		
//-----------------------------------------ANDERSEN or Berendsen THERMOSTAT------------------------------------------------------------------------------------	
		if 	(choose_thermostat==1){
			
		Andersen(vx2,vy2,vz2);}

		else if (choose_thermostat==2){
		
		vrescale(vx2,vy2,vz2);}
		
		 sumkin=0;
		 for (k=0;k<npart;++k)
		 {sumkin=sumkin+m[k]*vx2[k]*vx2[k]+m[k]*vy2[k]*vy2[k]+m[k]*vz2[k]*vz2[k];}
		 
		 
		K[n]=0.5*sumkin;
		 
		T[n]=(2*K[n])/(3*(mobile-1)*kbol);
		
		E[n]=K[n]+U[n];

		if (n%nstep==0){printf("%d\t%.16lf\t%.16lf\n",n,U[n]/npart,T[n]);}

		sumU = sumU + U[n]; //calculate the average
		sumKin = sumKin + K[n]; //calculate the average
		
//----------------------------------------------end of THERMOSTATS-----------------------------------------------------------------------------------------		
		
		
//---------------------------------start the ACF calculation for Total DOS----------------------------------------------------------------------------------

		acf(vx2,vy2,vz2);
		
//****************************  update **************************************************
		
		for(l=0;l<npart;++l)
		{
			x1[l]=x2[l];
			y_1[l]=y2[l];
			z1[l]=z2[l];
			
			vx1[l]=vx2[l];
			vy_1[l]=vy2[l];
			vz1[l]=vz2[l];
			
			ax1[l]=ax2[l];
			ay_1[l]=ay2[l];
			az1[l]=az2[l];
			
		}
		
//_________________paraview or jmol________________________________________________________

	if (n%graph_step==0) {

		if (csvout==1){paraview_func(x1,y_1,z1);}

				
		//-----------jmol animation------------------------//
		fp=fopen(path_animation,"a");
		fprintf(fp,"%d\n%s\n",npart,title);
		for (k=0;k<npart;++k) 
		{
			fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species[k],x1[k],y_1[k],z1[k]);
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
				ginst[p] = ((Lx*Ly*Lz*hinst[p])/(4*pi*rinst[p]*rinst[p]*dr*npart*npart)); // we calc rdf for each time step
				
				fprintf(fp,"%.16lf\t%.16lf\n",rinst[p],ginst[p]);
				hinst[p]=0;				
			}
			
			fclose(fp);

		} //end monitor		
		
		
		
		}  //*********************************end time of simulation*****************************************************************************//

		if (jmolout==1){jmol_func(x1,y_1,z1);} //export the last .xyz file for jmol

//---------------Print averages------------------------------------//
		sumU = (sumU/N)/npart;
		sumKin = (sumKin/N)/npart;

		printf("\n\nThe averages are:\t  ");
		printf("%lf\t%lf\n\n",sumU,sumKin);

//------------------------------------------------------------------//
	
	time(&rawtime);
	sprintf (stoptime, "\nSimulation ended at: %s",ctime(&rawtime));
	
	printf ("%s\n", starttime);
	printf ("%s\n", stoptime);
	
// call RDF
	
	RDF();
	
	
//------------------------EXPORT VARIOUS DATA---------------------------------------------------------------------------------------------------------------	
	
	export_info();

	
//_______________________free vectors from memory________________________________________________________________________________________________________________
	
	free(m);free(Z);free(covalent_r);free(vdW_r);
	free(Tinit);free(Tfinal);free(thermostat);
		
	free(myA);free(myB);free(lamda);free(mi);free(beta);
	free(ni);free(ci);free(di);free(hi);free(myR);free(myS);
	free(dx_zeta_III);free(dy_zeta_III);free(dz_zeta_III);

	free(myepsilon);free(mysigma);
	
	free(x1); free(y_1); free(z1); 
	free(status); free(thermostatus); 
	free(vx1); free(vy_1); free(vz1);
	free(ax1); free(ay_1); free(az1);	
	free(U); free(K); free(E); free(T);
	free(Rep); free(Attr); free(Bond);
	free(x2); free(y2); free(z2);
	free(vx2); free(vy2); free(vz2);  
	free(ax2); free(ay2); free(az2); 
	
	free(vxfixed); free(vyfixed); free(vzfixed); free(acf_t); free(sumtime);
	free(h); free(g); free(r1); free(hinst); free(ginst); free(rinst);
	
	free(sumfx);free(sumfy);free(sumfz);
	
	for (i=0; i<npart; ++i)
	{
		free(species[i]);
	}
	
	free(species);


	if (linkcells==1)
	{
		
		
		
		if (mypotential==0){

		

		free(neighborcellpartarray);
		free(allpartarray);

		free(cellpartarray);
		

		}
		if (mypotential==1){
		
		
		free(partindex);
		free(nextpartindex);	
		
		}

		
		for (i = 0; i < ghostsy; ++i)
		{
		free(cells[i]);
		}
		free(cells);
		
		for (i = 0; i < 26; ++i)
		{
		free(neighbours[i]);
		}
		free(neighbours);
		
		free(HEAD);
		free(LIST);

	
	}

	if (OS==1)	{_CrtDumpMemoryLeaks();}

	getchar();getchar();
	return 0;
	
	
}	// end main