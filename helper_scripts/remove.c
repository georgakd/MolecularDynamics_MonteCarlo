#define equal_diff 10e-4
#define cmax_length 1000
#ifdef _WIN32
#define OS 1
#include<direct.h>
#elif _WIN64
#define OS 1
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

FILE *fp;
int equal(double A,double B);
int main(int argc, char *argv[])
{
	int i,j,flag,mode;
	int *positions,length;
	char current_folder[cmax_length],file_path[cmax_length];
	int particles;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double *x,*y,*z;//*vx,*vy,*vz;
	int *mobileFlag,*thermoFlag;
	char **species,title[cmax_length], buffer[cmax_length];
	
	double *x_,*y_,*z_;
	char **species_;
	int particles_;
	
	printf("\n"
		"run 'remove' to generate the original XYZ file from which you will select the atoms to be removed\n"
		"run 'remove -1' to apply the removal\n"
		"\n\n");


//	char buffer[cmax_length];

	getcwd(current_folder,cmax_length);
	
	length=argc-1;
	
	if(length==1&&atoi(argv[1])==-1){mode=1;}else{mode=0;}
	
	if(mode==1)
	{	
		sprintf(file_path,"C:\\jmol\\out\\selection.xyz"); //the folder where the selected atoms from jmol are stored
		fp=fopen(file_path,"r");
		
		if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}

		fscanf(fp,"%d\n\n",&particles_);
		
		species_=(char**)malloc(particles_*sizeof(char*));
		for (i=0;i<particles_;++i)
		{
			species_[i]=(char*)malloc(4*sizeof(char));
		}
		x_=(double*)malloc(particles_*sizeof(double));
		y_=(double*)malloc(particles_*sizeof(double));
		z_=(double*)malloc(particles_*sizeof(double));
		
		for(i=0;i<particles_;++i){fscanf(fp,"%s\t%lf\t%lf\t%lf\n",species_[i],&x_[i],&y_[i],&z_[i]);}
		
		fclose(fp);
	}	
	
	if(mode==0)
	{
		positions=(int*)malloc(length*sizeof(int));

		for(i=0;i<length;++i)
		{
			positions[i]=atoi(argv[i+1]);
		}
	}
	
	sprintf(file_path,"%s/positions.dat",current_folder);
	fp=fopen(file_path,"r");	

	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}

	fscanf(fp,"title\t%s\nparticles\t%d\n",title,&particles);
	
	species=(char**)malloc(particles*sizeof(char*));
	for (i=0;i<particles;++i)
	{
		species[i]=(char*)malloc(4*sizeof(char));
	}
	x=(double*)malloc(particles*sizeof(double));
	y=(double*)malloc(particles*sizeof(double));
	z=(double*)malloc(particles*sizeof(double));
	mobileFlag=(int*)malloc(particles*sizeof(int));
	thermoFlag=(int*)malloc(particles*sizeof(int));
//	vx=(double*)malloc(particles*sizeof(double));
//	vy=(double*)malloc(particles*sizeof(double));
//	vz=(double*)malloc(particles*sizeof(double));
	
	fscanf(fp,"supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
	//fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");
	fgets(buffer,cmax_length,fp);
	//fgets(buffer,cmax_length,fp);
	for(i=0;i<particles;++i)
	{
		fgets(buffer,cmax_length,fp);
		sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%d\t%d",species[i],&x[i],&y[i],&z[i],&mobileFlag[i],&thermoFlag[i]);
		//fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",species[i],&x[i],&y[i],&z[i],&mobileFlag[i],&thermoFlag[i],&vx[i],&vy[i],&vz[i]);
	}
	fclose(fp);
	
	sprintf(file_path,"%s/positions_removed.dat",current_folder);
	fp=fopen(file_path,"w+");
	
	if(mode==1){fprintf(fp,"title\t%s\nparticles\t%d\n",title,particles-particles_);}
	if(mode==0){fprintf(fp,"title\t%s\nparticles\t%d\n",title,particles-length);}
	fprintf(fp,"supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
	//fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\tvx\tvy\tvz\n");
	fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\n");
	if(mode==1)
	{
		for(i=0;i<particles;++i)
		{
			flag=0;
			for(j=0;j<particles_;++j){if(equal(x[i],x_[j])==1&&equal(y[i],y_[j])==1&&equal(z[i],z_[j])==1){flag=1;}}
			if(flag==0)
			{
				//fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i],mobileFlag[i],thermoFlag[i],vx[i],vy[i],vz[i]);
				fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i],x[i],y[i],z[i],mobileFlag[i],thermoFlag[i]);
			}
		}
	}
	
	if(mode==0)
	{
		for(i=1;i<particles+1;++i)
		{
			flag=0;
			for (j=0;j<length;++j)
			{
				if(i==positions[j]){flag=1;}
			}
			if (flag==0)
			{
					//fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",species[i-1],x[i-1],y[i-1],z[i-1],mobileFlag[i-1],thermoFlag[i-1],vx[i-1],vy[i-1],vz[i-1]);	
				fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i-1],x[i-1],y[i-1],z[i-1],mobileFlag[i-1],thermoFlag[i-1]);	
			}
		}
	}
	
	fclose(fp);
	
	sprintf(file_path,"%s/XYZ.dat",current_folder);
	fp=fopen(file_path,"w+");
	
	if(mode==1){fprintf(fp,"%d\n%s\n",particles-particles_,title);}
	if(mode==0){fprintf(fp,"%d\n%s\n",particles-length,title);}
	
	if(mode==1)
	{
		for(i=0;i<particles;++i)
		{
			flag=0;
			for(j=0;j<particles_;++j){if(equal(x[i],x_[j])==1&&equal(y[i],y_[j])==1&&equal(z[i],z_[j])==1){flag=1;}}
			if(flag==0){fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);}
		}
	}
	
	if(mode==0)
	{
		for(i=1;i<particles+1;++i)
		{
			flag=0;
			for (j=0;j<length;++j)
			{
				if(i==positions[j]){flag=1;}
			}
			if (flag==0)
			{
					fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species[i-1],x[i-1],y[i-1],z[i-1]);
			}
		}
	}
	
	fclose(fp);
	
	for (i=0;i<particles;++i){free(species[i]);}
	free(species);free(x);free(y);free(z);free(mobileFlag);free(thermoFlag);
	//free(vx);free(vy);free(vz);
	if(mode==0){free(positions);}
	
	if(mode==1){
	for (i=0;i<particles_;++i){free(species_[i]);}
	free(species_);free(x_);free(y_);free(z_);
	}
	
	return 0;
}
int equal(double A,double B);
int equal(double A,double B)
{
	int value;
	if(fabs(A-B)<equal_diff){value=1;}else{value=0;}
	return value;
}
