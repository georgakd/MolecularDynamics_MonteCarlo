#define cmax_length 1000
#define species_length 4
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
int main(int argc,char *argv[])
{
	FILE *fp;
	char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length];
	double xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz;
	int i,particles;
	char **species;
	double *x,*y,*z;
	
	if(argc==1)
	{
		printf(
		"\n*Applies a periodic expansion on the xy plane.\n"
		"*Reads the supercell parameters from the init_file (standard MD format).\n"
		"*Prints data on screen - use pipeline to file.\n"
		"\n./xy_expand <xyz_file> <init_file>\n\n"
		);
		exit(-1);
	}
	
	getcwd(current_folder,cmax_length);
	
	sprintf(file_path,"%s/%s",current_folder,argv[2]);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"xmin\t%lf",&xmin);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"xmax\t%lf",&xmax);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"ymin\t%lf",&ymin);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"ymax\t%lf",&ymax);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"zmin\t%lf",&zmin);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"zmax\t%lf",&zmax);
	Lx=xmax-xmin;
	Ly=ymax-ymin;
	Lz=zmax-zmin;
	fclose(fp);
	
	sprintf(file_path,"%s/%s",current_folder,argv[1]);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&particles);
	fgets(buffer,cmax_length,fp);
	species=(char**)malloc(particles*sizeof(char*));
	for(i=0;i<particles;++i){species[i]=(char*)malloc(species_length*sizeof(char));}
	x=(double*)malloc(particles*sizeof(double));
	y=(double*)malloc(particles*sizeof(double));
	z=(double*)malloc(particles*sizeof(double));
	for(i=0;i<particles;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf\t%lf",species[i],&x[i],&y[i],&z[i]);
	}
	fclose(fp);

	printf("%d\n\n",particles*9);
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i]+Lx,y[i],z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i]-Lx,y[i],z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i]+Ly,z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i]-Ly,z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i]+Lx,y[i]+Ly,z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i]-Lx,y[i]-Ly,z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i]+Lx,y[i]-Ly,z[i]);}
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i]-Lx,y[i]+Ly,z[i]);}

	free(x);free(y);free(z);
	for(i=0;i<particles;++i){free(species[i]);}free(species);
	return 0;

}
