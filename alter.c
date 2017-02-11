#define equal_diff 10e-4
#define cmax_length 1000
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<string.h>
FILE *fp;
int equal(double A,double B);
int main(int argc, char *argv[])
{
	int i,j,flag,mode;
	int *positions,length;
	char current_folder[cmax_length],file_path[cmax_length];
	int particles;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double *x,*y,*z;
	int *mobileFlag,*thermoFlag;
	char **species,title[cmax_length];
	
	double *x_,*y_,*z_;
	char **species_;
	int particles_;
	
	char buffer[cmax_length];
	
	int type;
	
	char label[10];
	
	//------------------------------------------------------------------
	
	// print exe info
	if(argc==1)
	{
		printf("\n./alter <mode> <filename> <atoms-list>\n\n"
		"mode=1: remove\n"
		"mode=2: fix\n"
		"Options:\n"
		"- no atoms list defined:\n  strips from the input file the velocities info (if any) and generates xyz preview.\n"
		"- first element of the atoms list equal to -1:\n  reads atoms.xyz from the home directory and removes or fixes the atoms defined in the file.\n"
		"- normal atoms list:\n  removes or fixes the listed atoms.\n\n"
		"Generates <filename>.out for MD input and preview.dat for Jmol rendering.\n\n"
		);
		exit(-1);
	}
	
	// type: 1 remove, 2 fix
	type=atoi(argv[1]);
	
	// read current working directory
	getcwd(current_folder,cmax_length);
	
	// number of arguments
	length=argc-3;
	
	// select mode:
	// mode=1 --> use atoms.xyz exported by Jmol after manual selection
	// mode=0 --> define atoms from arguments
	if(length==1&&atoi(argv[3])==-1){mode=1;}else{mode=0;}
	
	if(mode==1)
	{
		// the file must be in the home directory	
		fp=popen("echo ~","r");												// home path
		fgets(buffer,cmax_length,fp);
		pclose(fp);
		buffer[strcspn(buffer,"\n")]='\0';
		sprintf(file_path,"%s/atoms.xyz",buffer);	
		fp=fopen(file_path,"r");
		// file check
		if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
		// read data
		fscanf(fp,"%d\n\n",&particles_);
		// preallocations
		species_=(char**)malloc(particles_*sizeof(char*));
		for (i=0;i<particles_;++i)
		{
			species_[i]=(char*)malloc(4*sizeof(char));
		}
		x_=(double*)malloc(particles_*sizeof(double));
		y_=(double*)malloc(particles_*sizeof(double));
		z_=(double*)malloc(particles_*sizeof(double));
		// read data
		for(i=0;i<particles_;++i){fscanf(fp,"%s\t%lf\t%lf\t%lf\n",species_[i],&x_[i],&y_[i],&z_[i]);}
		
		fclose(fp);
	}	
	
	if(mode==0)
	{
		// preallocate
		positions=(int*)malloc(length*sizeof(int));
		// save arguments
		for(i=0;i<length;++i)
		{
			positions[i]=atoi(argv[i+3]);
		}
	}
	
	// open init file
	sprintf(file_path,"%s/%s",current_folder,argv[2]);
	fp=fopen(file_path,"r");
	// file check
	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	// read data
	fscanf(fp,"title\t%s\nparticles\t%d\n",title,&particles);
	// preallocations
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
	// read data
	fscanf(fp,"supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
	//fscanf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\n");
	fgets(buffer,cmax_length,fp);
	for(i=0;i<particles;++i)
	{
		//fscanf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i],&x[i],&y[i],&z[i],&mobileFlag[i],&thermoFlag[i]);
		fgets(buffer,cmax_length,fp);
		sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%d\t%d",species[i],&x[i],&y[i],&z[i],&mobileFlag[i],&thermoFlag[i]);
	}
	fclose(fp);
	
	// write output to file
	sprintf(file_path,"%s/%s.out",current_folder,argv[2]);
	fp=fopen(file_path,"w+");
	
	if(type==1)
	{
		if(mode==1){fprintf(fp,"title\t%s\nparticles\t%d\n",title,particles-particles_);}
		if(mode==0){fprintf(fp,"title\t%s\nparticles\t%d\n",title,particles-length);}
	}
	else if(type==2)
	{
		fprintf(fp,"title\t%s\nparticles\t%d\n",title,particles);
	}
	fprintf(fp,"supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
	fprintf(fp,"species\tx\ty\tz\tmobileFlag\tthermoFlag\n");
	
	if(type==1)
	{
	
		if(mode==1)
		{
			for(i=0;i<particles;++i)
			{
				flag=0;
				for(j=0;j<particles_;++j){if(equal(x[i],x_[j])==1&&equal(y[i],y_[j])==1&&equal(z[i],z_[j])==1){flag=1;}}
				if(flag==0){fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i],x[i],y[i],z[i],mobileFlag[i],thermoFlag[i]);}
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
						fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i-1],x[i-1],y[i-1],z[i-1],mobileFlag[i-1],thermoFlag[i-1]);	
				}
			}
		}
	
	}
	
	else if (type==2)
	{
		
		if(mode==1)
		{
			for(i=0;i<particles;++i)
			{
				flag=0;
				for(j=0;j<particles_;++j){if(equal(x[i],x_[j])==1&&equal(y[i],y_[j])==1&&equal(z[i],z_[j])==1){flag=1;}}
				if(flag==0){fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i],x[i],y[i],z[i],mobileFlag[i],thermoFlag[i]);}
				else
				{fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i],x[i],y[i],z[i],0,0);}
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
						fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i-1],x[i-1],y[i-1],z[i-1],mobileFlag[i-1],thermoFlag[i-1]);	
				}
				else
				{
						fprintf(fp,"%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i-1],x[i-1],y[i-1],z[i-1],0,0);	
				}
			}
		}
		
	}
		
	fclose(fp);
	
	// create Jmol xyz file
	sprintf(file_path,"%s/preview.dat",current_folder);
	fp=fopen(file_path,"w+");
	
	if (type==1)
	{
	
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
	
	}
	
	else if (type==2)
	{
		
		fprintf(fp,"%d\n%s\n",particles,title);
	
		for(i=0;i<particles;++i)
		{
			fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);
		}
	
	}
	fclose(fp);
	
	//------------------------------------------------------------------
	
	for (i=0;i<particles;++i){free(species[i]);}
	free(species);free(x);free(y);free(z);free(mobileFlag);free(thermoFlag);
	
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
