#define cmax_length 1000
#define species_length 4
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

int main(int argc,char **argv)
{
	
	FILE *fp;
	
	char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length];
	
	int m,n,i,j,k,c,t;
	double top,h,d,x0,y0,S;
	
	int particles;
	double *x,*y,*z;
	char **species;
	
	double *xmesh,*ymesh,*zmesh;
	
	double r2;
	int b;
	
	int *flag;
	int s_particles;
	
	int *id,col_particles;
	
	if(argc==1)
	{
		printf("\n./z_projection <filename> <top> <dz> <x0> <y0> <d> <m> <n> <S> <N>\n\n");
		printf("top:\tinitial z position of the mesh above the substrate\n");
		printf("dz:\tlowering step\n");
		printf("x0,y0:\tmesh origin\n");
		printf("d:\tmesh spacing\n");
		printf("m,n:\tmesh points along the x and y directions\n");
		printf("S:\tatomic radius\n");
		printf("N:\tnumber of iterations\n\n");
		exit(-1);
	}
	
	top=atof(argv[2]);
	h=atof(argv[3]);
	x0=atof(argv[4]);
	y0=atof(argv[5]);
	d=atof(argv[6]);
	m=atof(argv[7]);
	n=atof(argv[8]);
	S=atof(argv[9]);
	c=atoi(argv[10]);

	getcwd(current_folder,cmax_length);
	sprintf(file_path,"%s/%s",current_folder,argv[1]);
	/*
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}

	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&particles);
	x=(double*)malloc(particles*sizeof(double));
	y=(double*)malloc(particles*sizeof(double));
	z=(double*)malloc(particles*sizeof(double));
	flag=(int*)malloc(particles*sizeof(int));
	id=(int*)malloc(particles*sizeof(int));
	species=(char**)malloc(particles*sizeof(char*));
	for(i=0;i<particles;++i)
	{
		species[i]=(char*)malloc(species_length*sizeof(char));
	}
	fgets(buffer,cmax_length,fp);
	for(i=0;i<particles;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf\t%lf",species[i],&x[i],&y[i],&z[i]);
	}

	fclose(fp);
	*/ 
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}

	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"particles\t%d",&particles);
	for(i=0;i<8;++i){fgets(buffer,cmax_length,fp);}
	x=(double*)malloc(particles*sizeof(double));
	y=(double*)malloc(particles*sizeof(double));
	z=(double*)malloc(particles*sizeof(double));
	flag=(int*)malloc(particles*sizeof(int));
	id=(int*)malloc(particles*sizeof(int));
	species=(char**)malloc(particles*sizeof(char*));
	for(i=0;i<particles;++i)
	{
		species[i]=(char*)malloc(species_length*sizeof(char));
	}
	
	for(i=0;i<particles;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf\t%lf",species[i],&x[i],&y[i],&z[i]);
	}

	fclose(fp);
	
	xmesh=(double*)malloc(m*n*sizeof(double));
	ymesh=(double*)malloc(m*n*sizeof(double));
	zmesh=(double*)malloc(m*n*sizeof(double));
	
	k=-1;
	for(i=0;i<m;++i)
	{
		for(j=0;j<n;++j)
		{
			k=k+1;
			xmesh[k]=x0+d*i;
			ymesh[k]=y0+d*j;
			zmesh[k]=top;
		}
	}
	
	for(i=0;i<particles;++i){flag[i]=0;}
	
	//
	for(k=0;k<m*n;++k)
	{
		col_particles=-1;
		for(i=0;i<particles;++i)
		{
			if((x[i]>=xmesh[k]-S && x[i]<=xmesh[k]+S) && (y[i]>=ymesh[k]-S && y[i]<=ymesh[k]+S))
			{
				col_particles=col_particles+1;
				id[col_particles]=i;
			}
		}
		
		
		b=0;
		for(t=0;t<c;++t)
		{
			zmesh[k]=top-h*t;
			
			//for (i=0;i<particles;++i)
			for(j=0;j<col_particles;++j)
			{
				i=id[j];
				r2=(x[i]-xmesh[k])*(x[i]-xmesh[k])+(y[i]-ymesh[k])*(y[i]-ymesh[k])+(z[i]-zmesh[k])*(z[i]-zmesh[k]);
				if(r2<=S*S){b=1;flag[i]=1;break;}
			}
			if(b==1){break;}
			
		}
	}
	
	printf("%d\n\n",particles+m*n);
	for(i=0;i<particles;++i){printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);}
	for(i=0;i<m*n;++i){printf("point\t%lf\t%lf\t%lf\n",xmesh[i],ymesh[i],zmesh[i]);}
	
	s_particles=0;
	for(i=0;i<particles;++i)
	{
		if(flag[i]==1){s_particles=s_particles+1;}
	}
	
	printf("%d\n\n",s_particles);
	for(i=0;i<particles;++i)
	{
		if(flag[i]==1){
		printf("%s\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i]);}
	}
	
	sprintf(file_path,"%s/surface.dat",current_folder);
	fp=fopen(file_path,"w+");
	for(i=0;i<m*n;++i){fprintf(fp,"%lf\t%lf\t%lf\n",xmesh[i],ymesh[i],zmesh[i]);}
	fclose(fp);
	
		
	free(x);
	free(y);
	free(z);
	for(i=0;i<particles;++i){free(species[i]);}free(species);
	
	free(xmesh);
	free(ymesh);
	free(zmesh);
	
	free(flag);
	
	free(id);
	
	return 0;
}
