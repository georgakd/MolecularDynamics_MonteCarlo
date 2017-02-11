#define cmax_length 1000
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<string.h>

FILE *fp;


int main(int argc, char *argv[])
{	
	
	char path[cmax_length],current_folder[cmax_length],buffer[cmax_length],title[cmax_length];
	int i,particles,*mobileFlag,*thermoFlag;
	char **species;
	double *x,*y,*z,xmin,xmax,ymin,ymax,zmin,zmax,Lx,a,exx;
	char dummy[cmax_length];
	
	// read current working directory
	getcwd(current_folder,cmax_length);
	
	
	// read arguments
	if (argc==1) {printf (" This script takes the following arguments:  \n 1) file_name \n 2) a (strain = a-1) \n");exit(-1);}
		
	// read arguments
	a=atof(argv[2]);
	
	
	// read init.dat
	sprintf(path,"%s/%s",current_folder,argv[1]);
	fp = fopen(path,"r");
	
	if(fp == NULL) {
      printf("Error opening file");
      return(-1);
   }
   
 
   fgets(buffer,cmax_length,fp); // title
   //printf("%s",buffer);
    sscanf(buffer,"%s\t%s",dummy,title);
   
   fgets(buffer,cmax_length,fp); // particles
   //printf("%s",buffer);   
   sscanf(buffer,"%s\t%d",dummy,&particles);
   
   fgets(buffer,cmax_length,fp); // supercell
   //printf("%s",buffer);
   
   fgets(buffer,cmax_length,fp); // xmin
  // printf("%s",buffer);   
   sscanf(buffer,"%s\t%lf",dummy,&xmin);
      fgets(buffer,cmax_length,fp); // xmax
  // printf("%s",buffer);   
   sscanf(buffer,"%s\t%lf",dummy,&xmax);
      fgets(buffer,cmax_length,fp); // ymin
  // printf("%s",buffer);   
   sscanf(buffer,"%s\t%lf",dummy,&ymin);
      fgets(buffer,cmax_length,fp); // ymax
 //  printf("%s",buffer);   
   sscanf(buffer,"%s\t%lf",dummy,&ymax);
      fgets(buffer,cmax_length,fp); // zmin
  // printf("%s",buffer);   
   sscanf(buffer,"%s\t%lf",dummy,&zmin);
      fgets(buffer,cmax_length,fp); // zmax
 //  printf("%s",buffer);   
   sscanf(buffer,"%s\t%lf",dummy,&zmax);
   
    fgets(buffer,cmax_length,fp); // coords header
 //  printf("%s",buffer);
   
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
   
   
   for (i=0;i<particles;++i){
	   
	    fgets(buffer,cmax_length,fp); 
		sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%d\t%d",species[i],&x[i],&y[i],&z[i],&mobileFlag[i],&thermoFlag[i]);
	      
	   }
   
   fclose(fp);
   
   sprintf(path,"%s/preview.dat",current_folder);
   fp=fopen(path,"w+");
   fprintf(fp,"%d\n\n",particles);
   for(i=0;i<particles;++i){fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);}
   
   
  
   //choose elongation percentage, if a=1 there is no elongation
    
   xmin = a*xmin;
   xmax = a*xmax;
   Lx = a*Lx;
   
   for (i=0;i<particles;++i) {x[i]=a*x[i];}
   fprintf(fp,"%d\n\n",particles);
   for(i=0;i<particles;++i){fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);} 
   fclose(fp);
   
	//write final file to console
	
	printf("title\t%s\nparticles\t%d\n",title,particles);
	
	printf("supercell\nxmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
	printf("species\tx\ty\tz\tmobileFlag\tthermoFlag\n");
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\t%d\t%d\n",species[i],x[i],y[i],z[i],mobileFlag[i],thermoFlag[i]);}
	
	//write output strain
	printf("The strain is: ");
	exx=a-1;
	printf ("%lf\n",exx);
	
	
   
   
   
   for(i=0;i<particles;++i){free(species[i]);}
   free(species);
   free(x);free(y);free(z);free(mobileFlag);free(thermoFlag);


}
