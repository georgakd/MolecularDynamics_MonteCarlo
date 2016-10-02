
//********************************************* MONTE CARLO PROJECT *************************************************//
//		RELEASE MODE v.2 17/12/2012																					    //
//		Created by Dimitra Georgakaki																			   //																					
//****************************************************************************************************************//


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

double random_gen(void);  //call the Wichmann-Hill random number generator
void parameters(void); // declare function to initialize my parameters
void atomic_numbers(void); 
void read_elements(void);
void linked_cells(int cellsx, int cellsy, int cellsz);

//-------- declare the function that will calculate interaction forces and potential

void tersoff_forces_cells(double *x,double *y,double *z);
void tersoffopt(double *x,double *y,double *z);
void lj_forces(double *x,double *y,double *z);
void lj_forces_cells(double *x,double *y,double *z);

//-------- declare functions for imaging

void jmol_func(double *x,double *y, double *z);
void paraview_func(double *x,double *y, double *z);



int main(void)
{
	
//*************************************************** Define local vars / parameters ************************************************************************//
	
	char mycommand[cmax_length]; // for mkdir and delete files and folders
	char starttime[cmax_length],stoptime[cmax_length];	// for counting the simulation
	char sim_folder[cmax_length]; // for defining folders
	char path[cmax_length],pathout[cmax_length],path_animation[cmax_length];
	double sumvx,sumvy,sumvz;
	double myrandom1,myrandom2;
	
	time_t rawtime;
	time(&rawtime);
	sprintf (starttime, "\nSimulation started at: %s",ctime(&rawtime));

//--------------------------------------------Read current directory--------------------------------------//
	
if (OS==1)
	{_getcwd(sim_folder,cmax_length);}
else 
	{getcwd(sim_folder,cmax_length);}
	

//-------------------------------------------------- Creating folders and paths ------------------------------------------------------------------------------//	
	
if (OS==1)
	{
	//sprintf(sim_folder,"%s","C:\\MD_files\\Silicon_MC_bulk"); // choose the simulation folder
	sprintf(file_path,"%s\\config.txt",sim_folder);
	sprintf(pathcsv,"%s\\csv\\",sim_folder);
	sprintf(pathout,"%s\\out\\",sim_folder);
	sprintf(pathxyz,"%s\\jmol\\",sim_folder);
	sprintf(path_animation,"%s\\out\\animation.xyz",sim_folder);
	}
	else 
	{
	//sprintf(sim_folder,"%s","/Users/dimitrageorgakaki/Documents/C/MD_files/Silicon_MC_bulk"); // choose the simulation folder
	sprintf(file_path,"%s/config.txt",sim_folder);
	sprintf(pathcsv,"%s/csv/",sim_folder);
	sprintf(pathout,"%s/out/",sim_folder);
	sprintf(pathxyz,"%s/jmol/",sim_folder);
	sprintf(path_animation,"%s/out/animation.xyz",sim_folder);
	}
	
	sprintf(mycommand,"mkdir %s",pathcsv);
	system(mycommand);
	sprintf(mycommand,"mkdir %s",pathout);
	system(mycommand);
	sprintf(mycommand,"mkdir %s",pathxyz);
	system(mycommand);
	
	
	printf("\n\n--------------Define the simulation folder: ");printf("%s",sim_folder);printf("-----\n\n");

	
	parameters();  // call the function parameters.c
	
//----------------------------------------------------Delete previous folders csv,jmol,radial,forces if they exist-------------------------------------------//
	
	printf("\nDeleting all previous existing files in the selected folders!! \n ");
	
	if (OS==1){
		sprintf(mycommand,"del /Q %s*.csv",pathcsv);system(mycommand);
		sprintf(mycommand,"del /Q %s*.xyz",pathxyz);system(mycommand);
	}
	else{
		sprintf(mycommand,"find %s -name '*.csv' | xargs rm",pathcsv);system(mycommand);
		sprintf(mycommand,"find %s -name '*.xyz' | xargs rm",pathxyz);system(mycommand);
	}

//-----------------------------------------------Read box dimensions, number of particles, init surface positions---------------------------------------------//
	
	
//----------------------Sample information--------------------------//
	
if (resume_flag==0){ 	// import for the first time initial positions
		if (OS==1){sprintf(path,"%s\\positions.dat",sim_folder);}
		else {sprintf(path,"%s/positions.dat",sim_folder);} 
		fp=fopen(path,"r");
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
		vx=(double*)malloc((npart)*sizeof(double));
		vy=(double*)malloc((npart)*sizeof(double));
		vz=(double*)malloc((npart)*sizeof(double));
			
		for (k=0;k<npart;++k){		
		fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[k],&x1[k],&y_1[k],&z1[k],&status[k],&thermostatus[k]);}
		fclose(fp);
	}
	
	
if (resume_flag==1){ // resume simulation (from MD or from MC) for better surface relaxation
		
	if (OS==1){sprintf(path,"%s\\resume.dat",pathout);}
	else {sprintf(path,"%s/resume.dat",pathout);}
	fp=fopen(path,"r");
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
	vx=(double*)malloc((npart)*sizeof(double));
	vy=(double*)malloc((npart)*sizeof(double));
	vz=(double*)malloc((npart)*sizeof(double));


	for(k=0;k<npart;++k){
		fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",species[k],&x1[k],&y_1[k],&z1[k],&status[k],&thermostatus[k],&vx[k],&vy[k],&vz[k]);}
	fclose(fp);
	}

//-------------------------------------------------------- Allocation of all other arrays ------------------------------------------------------------------//	
	

	x2=(double*)malloc((npart)*sizeof(double));
	y2=(double*)malloc((npart)*sizeof(double));
	z2=(double*)malloc((npart)*sizeof(double));
	
	U=(double*)malloc((max_cycles+1)*sizeof(double)); //new configuration

	Z=(int*)malloc((npart)*sizeof(int));
	m=(double*)malloc((npart)*sizeof(double));
	covalent_r=(double*)malloc((npart)*sizeof(double));
	vdW_r=(double*)malloc((npart)*sizeof(double));
	
	

//------------------------------------------------Define secondary parameters and import initial crystal state-----------------------------------------------//
	
	//***********be very careful this is the correct order for the parameters to be read from the program*******************//
	
	
	if (mypotential==0){printf("\n  Tersoff T3 potential...\n");}
	if (mypotential==1){printf("\n  Lennard Jones potential...\n");}
	if (mypotential==2){printf("\n  Sutton Chen potential...\n");}
	
	mobile=0;thermo=0;
	for (k=0;k<npart;++k){mobile=mobile+status[k];} //count the mobile atoms
	for (k=0;k<npart;++k){thermo=thermo+thermostatus[k];} //count the thermostated atoms
	
	
	Lx=xmax-xmin; // box simulation
	Ly=ymax-ymin; 
	Lz=zmax-zmin;

	atomic_numbers();
	read_elements();	

	
//--------------------------------------------------------- Use or not use LINKED CELLS --------------------------------------------------------------------//
	
	
		// Divide space for interactions
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
	
	
	
//***************************************All initializations done, start the program!!!*********************************************************************//	
	
	printf("All initializations done, start the Monte Carlo program!!!\n");
	printf("\nTotal atoms: %d\n", npart);printf("The mobile atoms are: %d\n", mobile);printf("The thermostated atoms are: %d\n\n", thermo);

	
//***********************************Loop through all particles *******************************************************************************************//

	srand((unsigned)time(NULL)); // seed for the random particle (or srand(number)?????)


// 1. calculate the potential energy for the initial-crystal configuration	
	
	n=0;
	if (mypotential==0){
	
		if (linkcells==1){tersoff_forces_cells(x1,y_1,z1);}						
		if (linkcells==0){tersoffopt(x1,y_1,z1);}
		U[n]=0.5*sumenergy;
	}
		
	if (mypotential==1){
		
		if (linkcells==1){lj_forces_cells(x1,y_1,z1); U[n]=4*0.5*epsilon*sumVcells;} // we must multiply with 0.5 because we use the N*N scheme
		if (linkcells==0){lj_forces(x1,y_1,z1); U[n]=4*epsilon*sumV;}   // we do not multiply with 0.5 because we use the N*(N-1)/2 scheme
	}

	if (csvout==1){paraview_func(x1,y_1,z1);}

	fp=fopen(path_animation,"w+");
	fprintf(fp,"%d\n%s\n",npart,title);
	for (k=0;k<npart;++k) // n=0
	{
		fprintf(fp,"%s %.16lf %.16lf %.16lf\n",species[k],x1[k],y_1[k],z1[k]);
	}
	fclose(fp);

	printf("%d\t%.16lf\n",n,U[n]/npart);
	getchar();
/*	
//------------------------Sutton-Chen-------------------------------------------//	
	
	sumV=0;sum_rep=0;
	for (i=0;i<npart;++i)
	{	
		sum_attr=0;denfunc=0;
		for (j=0;j<npart;++j)
		{
			if (i!=j)
			{
			sigma=0.5*(mysigma[i]+mysigma[j]);
			epsilon=sqrt(myepsilon[i]*myepsilon[j]);
			rc = sigma;   // cutoff distance
			
			xforce=x1[j];yforce=y_1[j];zforce=z1[j];
			if (x1[j]-x1[i]>0.5*Lx) {xforce=x1[j]-Lx;}
			if (x1[j]-x1[i]<-0.5*Lx) {xforce=x1[j]+Lx;}
			if (y_1[j]-y_1[i]>0.5*Ly) {yforce=y_1[j]-Ly;}
			if (y_1[j]-y_1[i]<-0.5*Ly) {yforce=y_1[j]+Ly;}
			
			if (z1[j]-z1[i]>0.5*Lz) {zforce=z1[j]-Lz;}
			if (z1[j]-z1[i]<-0.5*Lz) {zforce=z1[j]+Lz;}
			
			r2=(xforce-x1[i])*(xforce-x1[i])+(yforce-y_1[i])*(yforce-y_1[i])+(zforce-z1[i])*(zforce-z1[i]);
			
			
			r=sqrt(r2);// r calculation
				
				//printf("%lf\t%lf\t%lf\t%lf\t%lf\n",r2,r,xforce,yforce,zforce);	
				
				//getchar();
	
				if (r<=rc){
					term = sigma/r;
					coeff_rep = pow(term,npower);
					sum_rep=sum_rep+coeff_rep;
					coeff_attr = pow(term,mpower);
					sum_attr = sum_attr + coeff_attr;
					//sumV = sumV + 0.5*epsilon*coeff_rep - epsilon*c_param*denfunc;

				}	
				
	
			}	// i not equal to j
	
		} //end of j
		denfunc = denfunc+sqrt(sum_attr);
	
		sumV = sumV + 0.5*epsilon*coeff_rep - epsilon*c_param*denfunc;

		
	} //end of i
	
*/	

	//U[0]=(sum_rep*epsilon*(2*npower-mpower))/(2*mpower);
	//printf("%lf\t",U[0]/npart);
	//getchar();
//------------------------Sutton-Chen-------------------------------------------//	

sumU=0; // average potential energy	
MC_accept=0; // count the accepted MC moves

for (n=1;n<=max_cycles;++n){     //****************** define the Monte Carlo cycles ******************//
	
		temp = rand()%npart; // 2. pick a random atom, at each MC cycle one atom is randomly picked to be moved

	
//***************** 3a. move one atom at a time by applying PBC ***********************************//
	
			if (status[temp]==1)  // if the atom is mobile, then it is allowed to make random displacements
			{ 
			myrandom_choicex=random_gen();
			random_dispx = (myrandom_choicex - 0.5)*dmax;
			xtemp = random_dispx+x1[temp];
				
				if (xtemp>xmax) {x2[temp]=xtemp-Lx;}	//PBC in x
				else if (xtemp<xmin) {x2[temp]=xtemp+Lx;}	//PBC in x 
				else x2[temp]=xtemp; //for every other case
				
			myrandom_choicey=random_gen();
			random_dispy = (myrandom_choicey - 0.5)*dmax;
			ytemp = random_dispy+y_1[temp];
				
				if (ytemp>ymax) {y2[temp]=ytemp-Ly;} //PBC in y
				else if (ytemp<ymin) {y2[temp]=ytemp+Ly;} //PBC in y
				else y2[temp]=ytemp; //for every other case
				
			myrandom_choicez=random_gen();
			random_dispz = (myrandom_choicez - 0.5)*dmax;
			ztemp = random_dispz+z1[temp];
				
				if (ztemp>zmax) {z2[temp]=ztemp-Lz;} //PBC in z
				else if (ztemp<zmin) {z2[temp]=ztemp+Lz;} //PBC in z
				else z2[temp]=ztemp; //for every other case
				
			}

			else if (status[temp]==0)  // if the atom is fixed, keep the configuration
			{ 
				x2[temp]=x1[temp];
				y2[temp]=y_1[temp];
				z2[temp]=z1[temp];
			}
	
		
//******************* 3b. for every atom keep positions, for the random atom store the displaced positions *************************//
	
	for (l=0;l<npart;++l) { 
			
		if (l!=temp) 
			{
				x2[l]=x1[l];
				y2[l]=y_1[l];
				z2[l]=z1[l];
				
			} 
	}	
	

	 
//************ end of particle loop (if I want, I can move all atoms at the same time, not recommended!!!!!!!!) ******************//
	
	if (mypotential==0){
	
		if (linkcells==1){tersoff_forces_cells(x2,y2,z2);}						
		if (linkcells==0){tersoffopt(x2,y2,z2);}
		Utrial=0.5*sumenergy;
	}
		
	if (mypotential==1){
		
		if (linkcells==1){lj_forces_cells(x2,y2,z2); Utrial=4*0.5*epsilon*sumVcells;} // we must multiply with 0.5 because we use the N*N scheme
		if (linkcells==0){lj_forces(x2,y2,z2); Utrial=4*epsilon*sumV;}   // we do not multiply with 0.5 because we use the N*(N-1)/2 scheme
	}

	
			dU = Utrial - U[n-1]; // calculate energy difference

//********************************** 5. Perform Metropolis Criterion *************************************************************************//
	
if (Utrial<= U[n-1]) { counter = 1; }
	else {
			Wvar = exp(-(Utrial-U[n-1])/(kbol*Tref));
			myrandom_choice=random_gen();
		if (Wvar>=myrandom_choice) {counter = 1;}
		else {counter = 0;}
	}
		
		
	if (counter==1) {
		//keep new positions and update positions
		for (l=0;l<npart;++l){
			x1[l]=x2[l];
			y_1[l]=y2[l];
		z1[l]=z2[l];}
		
		MC_accept=MC_accept+1;
		U[n]=Utrial; //new potential energy
	}	
	else {
		//restore old positions
		for (l=0;l<npart;++l){
		x2[l]=x1[l];
		y2[l]=y_1[l];
		z2[l]=z1[l];}
		
		U[n]=U[n-1]; //old potential energy
	
	
	}
			
		
		
// sample every nsample cycle
		
if (n%nsample==0){		
		printf("%d\t%lf\t%lf\t%.16lf\t%.16lf\t%.16lf\n",n,myrandom_choice,Wvar,U[n]/npart,Utrial/npart,dU/npart);
							
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


		//sumU = sumU + U[n]; //calculate the average potential energy ??? where to put it???
	
} //end of MC cycles

//assign random velocities
	if (resume_flag==0){ // if we do not resume from MD we have to assign random velocities
				
		sumvx=0;sumvy=0;sumvz=0;
		
		for (l=0;l<npart;++l)
		{
			if (status[l]==0){vx[l]=0;vy[l]=0;vz[l]=0;}
			
			else if (status[l]==1)
			{
				myrandom1=random_gen();
				myrandom2=random_gen();
				//myrandom1=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);		// Box-Muller method to create gaussian random numbers
				//myrandom2=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				vx[l]=sqrt((kbol*Tref)/m[l])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2); // Box-Muller Method
				vy[l]=sqrt((kbol*Tref)/m[l])*sqrt(-2*log(myrandom1))*sin(2*pi*myrandom2);
				
				myrandom1=random_gen();
				myrandom2=random_gen();
				
				//myrandom1=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);		// Box-Muller method to create gaussian random numbers
				//myrandom2=(double)((int)rand()+1)/((unsigned int)RAND_MAX+1);
				vz[l]=sqrt((kbol*Tref)/m[l])*sqrt(-2*log(myrandom1))*cos(2*pi*myrandom2);
				
				
				sumvx=sumvx+vx[l];sumvy=sumvy+vy[l];sumvz=sumvz+vz[l];  //velocity center of mass for mobile particles
				
			}
			
		}	
		
		sumvx=sumvx/mobile; sumvy=sumvy/mobile; sumvz=sumvz/mobile;  //   <v>
		
		for  (l=0;l<npart;++l)
		{				
			if (status[l]==1)
			{
				vx[l]=(vx[l]-sumvx); //rescaling velocities
				vy[l]=(vy[l]-sumvy);
				vz[l]=(vz[l]-sumvz);
			}
			
		}
	} // end of first run
	
	
if (jmolout==1){jmol_func(x1,y_1,z1);} //export the last .xyz file in jmol
	
	printf("\nThe accepted MC moves are:\t");
	printf("%d",MC_accept);
	
	
	time(&rawtime);
	sprintf (stoptime, "\nSimulation ended at: %s",ctime(&rawtime));
	
	printf ("%s\n", starttime);
	printf ("%s\n", stoptime);
	
	
	
//------------------------EXPORT VARIOUS DATA-------------------------------------------------------------------------------------------//

	sprintf(path,"%sresume.dat",pathout);
	fp=fopen(path,"w+");
	fprintf(fp,"title\t%s\nparticles\t%d\nsupercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",title,npart,xmin,xmax,ymin,ymax,zmin,zmax);
	fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");	
	
	for(k=0;k<npart;++k)
	{
		fprintf(fp,"%s\t%.16lf\t%.16lf\t%.16lf\t%d\t%d\t%lf\t%lf\t%lf\n",species[k],x2[k],y2[k],z2[k],status[k],thermostatus[k],vx[k],vy[k],vz[k]);
	}
	fclose(fp);
	
	
	
	sprintf(path,"%sthermodynamic_quantities.dat",pathout);
	fp=fopen(path,"w+");

	for(n=0;n<max_cycles;++n)
	{
		fprintf(fp,"%d\t%.16lf\n",n,U[n]/npart);
	}
	fclose(fp);

	
	
	
//_______________________free vectors from memory__________________________________________
	

	
	free(m);free(Z);free(covalent_r);free(vdW_r);
	free(myA);free(myB);free(lamda);free(mi);free(beta);
	free(ni);free(ci);free(di);free(hi);free(myR);free(myS);
	free(myepsilon);free(mysigma);
	free(vx);free(vy);free(vz);
	free(x1); free(y_1); free(z1); free(x2); free(y2); free(z2);
	free(status); free(thermostatus);	
	free(U); 	
	
	for (i=0; i<npart; ++i)
	{
		free(species[i]);
	}
	
	free(species);

		free(HEAD);
		free(LIST);
	
		if (mypotential==1){
		free(partindex);
		free(nextpartindex);
		}
		if (mypotential==0){
		free(cellpartarray);free(neighborcellpartarray);free(allpartarray);
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


				
		if(OS==1){_CrtDumpMemoryLeaks();}
		getchar();getchar();


	return 0;
	
	
}	// end main