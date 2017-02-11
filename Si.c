#define pi 3.14159265359
#define equal_diff 10e-4
#define cmax_length 1000
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

int equal(double A,double B);
int main(int argc,char *argv[])
{
	FILE *fp;
	
	char current_folder[cmax_length],file_path[cmax_length];
	
	double d=5.432;
	int xpart,ypart,zpart;
	double x0=0,y0=0,z0=0;

	int particlesFCC;
	int xypart,A;
	
	double *xFCC,*yFCC,*zFCC;
	double *x,*y,*z;

	double ux,uy,uz;
	double x1,y1,z1,x2,y2,z2;
	double x_0,y_0,z_0;
	
	double x_temp,y_temp,z_temp;
	double norm,rads,theta;
	
	double vacuum;

	int surf_rec;
	
	int i,j,k,l,particles;

	double xmin,xmax,ymin,ymax,zmin,zmax;
	double Lx,Ly,Lz;
	
	int fixed;
		
	//------------------------------------------------------------------
	// read current working directory
	getcwd(current_folder,cmax_length);
	// print exe info
	if(argc==1)
	{
		printf(
		"\n*Generates a diamond lattice of Si atoms with its [001] direction parallel to the z axis.\n"
		"*The final number of atoms is 8*xpart*ypart*zpart.\n"
		"*Lattice info: d=5.432A, origin: (0,0,0).\n"
		"*All atoms are initially considered as mobile and thermostated.\n"
		"*Generates Si.dat for Jmol and prints on screen the MD init file (use pipeline to file).\n"
		"*Additional info:\n"
		"\t- Bulk crystall: vacuum=0, surf_rec=0, fixed=0\n"
		"\t- Slab geometry: define the vacuum and the number of fixed layers\n"
		"\t- Set surf_rec=1 to create a 2x1 reconstruction on the upper side of the slab\n"
		"\t- Due to symmetry reasons, for the reconstructed surface, xpart and ypart must be even\n"
		"\n./Si <xpart> <ypart> <zpart> <vacuum> <surf_rec> <fixed>\n\n"
		);
		exit(-1);
	}
	// read arguments
	xpart=atoi(argv[1]);
	ypart=atoi(argv[2]);
	zpart=atoi(argv[3]);
	vacuum=atof(argv[4]);
	surf_rec=atoi(argv[5]);
	fixed=atoi(argv[6]);
	
	//------------------------------------------------------------------
	// build a bulk diamond crystal with z || [001]
	ypart=ypart*2;
	zpart=zpart*2;
	particlesFCC=xpart*ypart*zpart;
	particles=particlesFCC*2;

	x=(double*)malloc(particles*sizeof(double));
	y=(double*)malloc(particles*sizeof(double));
	z=(double*)malloc(particles*sizeof(double));
	xFCC=(double*)malloc(particlesFCC*sizeof(double));
	yFCC=(double*)malloc(particlesFCC*sizeof(double));
	zFCC=(double*)malloc(particlesFCC*sizeof(double));

	l=-1;
	for(k=0;k<zpart;++k)
	{
	  if(k%2==1)
	{
	  for(j=0;j<ypart;++j)
		{
		  if(j%2==1)
		{
		  for(i=0;i<xpart;++i)
			{
			  l=l+1;
			  xFCC[l]=x0+i*d;
			  yFCC[l]=y0+j*d/2;
			  zFCC[l]=z0+k*d/2;
			}
		}
		  else
		{
		  for(i=0;i<xpart;++i)
					{
					  l=l+1;
					  xFCC[l]=x0+i*d+d/2;
					  yFCC[l]=y0+j*d/2;
					  zFCC[l]=z0+k*d/2;
					}
		}

		}
	}
	  else
	{
	  for(j=0;j<ypart;++j)
			{
			  if(j%2==0)
		{
				  for(i=0;i<xpart;++i)
					{
					  l=l+1;
					  xFCC[l]=x0+i*d;
					  yFCC[l]=y0+j*d/2;
					  zFCC[l]=z0+k*d/2;
					}
				}
			  else
		{
				  for(i=0;i<xpart;++i)
					{
					  l=l+1;
					  xFCC[l]=x0+i*d+d/2;
					  yFCC[l]=y0+j*d/2;
					  zFCC[l]=z0+k*d/2;
					}
				}

			}
	}
	}
	xypart=xpart*ypart;
	A=2*xypart;
	k=-1;
	for(i=0;i<zpart;++i)
	{
	  for(j=0;j<xypart;++j)
	{
	  k=k+1;
	  x[A*i+j]=xFCC[k];
	  y[A*i+j]=yFCC[k];
	  z[A*i+j]=zFCC[k];
	  x[A*i+j+xypart]=xFCC[k]+d/4;
		  y[A*i+j+xypart]=yFCC[k]+d/4;
		  z[A*i+j+xypart]=zFCC[k]+d/4;

	}
	}
	
	//------------------------------------------------------------------
	// apply the 2x1 reconstruction
	if(surf_rec==1)
	{
		for(i=0;i<xpart;++i)
		{
			if(i%2==0){theta=33.3536;}else{theta=-33.3536;}
			//x0+d/2+i*d,y0,z0+(zpart/2)*d-d/4,x0+d/2+i*d+d/4,y0+d/4,z0+(zpart/2)*d-d/4
			x1=x0+d/2+i*d;
			y1=y0;
			z1=z0+(zpart/2)*d-d/4;
			x2=x1+d/4;
			y2=y1+d/4;
			z2=z1;
			x_0=x0+d/2+i*d;
			y_0=y0;
			z_0=z0+(zpart/2)*d-d/4-d/4;
			ux=1;
			uy=1;
			uz=0;
			norm=sqrt(ux*ux+uy*uy+uz*uz);
			ux=ux/norm;
			uy=uy/norm;
			uz=uz/norm;
			rads=(pi*theta)/180.0;
			//printf("Reference atom for rotation:\n[%lf\t%lf\t%lf]\nAtoms to define the line:\n[%lf\t%lf\t%lf]\t[%lf\t%lf\t%lf]\n",x_0,y_0,z_0,x1,y1,z1,x2,y2,z2);
			for(j=0;j<particles;++j)
			{
				if(equal((y2-y1)*(z[j]-z1),(z2-z1)*(y[j]-y1))==1 && equal((x2-x1)*(z[j]-z1),(z2-z1)*(x[j]-x1))==1 && equal((x2-x1)*(y[j]-y1),(y2-y1)*(x[j]-x1))==1)
				{

					x_temp=(ux*ux*(1.0-cos(rads))+cos(rads))*(x[j]-x_0)+(ux*uy*(1.0-cos(rads))-uz*sin(rads))*(y[j]-y_0)+(ux*uz*(1.0-cos(rads))+uy*sin(rads))*(z[j]-z_0)+x_0;
					y_temp=(ux*uy*(1.0-cos(rads))+uz*sin(rads))*(x[j]-x_0)+(uy*uy*(1.0-cos(rads))+cos(rads))*(y[j]-y_0)+(uy*uz*(1.0-cos(rads))-ux*sin(rads))*(z[j]-z_0)+y_0;
					z_temp=(ux*uz*(1.0-cos(rads))-uy*sin(rads))*(x[j]-x_0)+(uy*uz*(1.0-cos(rads))+ux*sin(rads))*(y[j]-y_0)+(uz*uz*(1.0-cos(rads))+cos(rads))*(z[j]-z_0)+z_0;
					
					//printf("**\t%lf\t%lf\t%lf\n\n",x[j],y[j],z[j]);
					x[j]=x_temp;
					y[j]=y_temp;
					z[j]=z_temp;
				}
			}//printf("----------\n");
		}
		
		for(i=0;i<ypart/2;++i)
		{
			if(i%2==1){theta=33.3536;}else{theta=-33.3536;}
			//x0,y0+d/2+i*d,z0+(zpart/2)*d-d/4,x0,y0+d/4+d/2+i*d+d/4,z0+(zpart/2)*d-d/4
			x1=x0;
			y1=y0+d/2+i*d;
			z1=z0+(zpart/2)*d-d/4;
			x2=x1+d/4;
			y2=y1+d/4;
			z2=z1;
			x_0=x0;
			y_0=y0+d/2+i*d;
			z_0=z0+(zpart/2)*d-d/4-d/4;
			ux=1;
			uy=1;
			uz=0;
			norm=sqrt(ux*ux+uy*uy+uz*uz);
			ux=ux/norm;
			uy=uy/norm;
			uz=uz/norm;
			rads=(pi*theta)/180.0;
			//printf("Reference atom for rotation:\n[%lf\t%lf\t%lf]\nAtoms to define the line:\n[%lf\t%lf\t%lf]\t[%lf\t%lf\t%lf]\n",x_0,y_0,z_0,x1,y1,z1,x2,y2,z2);
			for(j=0;j<particles;++j)
			{
				if(equal((y2-y1)*(z[j]-z1),(z2-z1)*(y[j]-y1))==1 && equal((x2-x1)*(z[j]-z1),(z2-z1)*(x[j]-x1))==1 && equal((x2-x1)*(y[j]-y1),(y2-y1)*(x[j]-x1))==1)
				{
					x_temp=(ux*ux*(1.0-cos(rads))+cos(rads))*(x[j]-x_0)+(ux*uy*(1.0-cos(rads))-uz*sin(rads))*(y[j]-y_0)+(ux*uz*(1.0-cos(rads))+uy*sin(rads))*(z[j]-z_0)+x_0;
					y_temp=(ux*uy*(1.0-cos(rads))+uz*sin(rads))*(x[j]-x_0)+(uy*uy*(1.0-cos(rads))+cos(rads))*(y[j]-y_0)+(uy*uz*(1.0-cos(rads))-ux*sin(rads))*(z[j]-z_0)+y_0;
					z_temp=(ux*uz*(1.0-cos(rads))-uy*sin(rads))*(x[j]-x_0)+(uy*uz*(1.0-cos(rads))+ux*sin(rads))*(y[j]-y_0)+(uz*uz*(1.0-cos(rads))+cos(rads))*(z[j]-z_0)+z_0;
					
					//printf("**\t%lf\t%lf\t%lf\n\n",x[j],y[j],z[j]);
					x[j]=x_temp;
					y[j]=y_temp;
					z[j]=z_temp;
				}
			}//printf("----------\n");
		}
	}
	//------------------------------------------------------------------
	// supercell
	
	xmin = x0 - d/8;
	xmax = x0 + (xpart - 1)*d + d/2 + d/4 + d/8;
	ymin = y0 - d/8;
	ymax = y0 + (ypart/2 - 1)*d + d/2 + d/4 + d/8;
	zmin = z0 - d/8 - vacuum;
	zmax = z0 + (zpart/2 - 1)*d + d/2 + d/4 + d/8 + vacuum;
	Lx=xmax-xmin;
	Ly=ymax-ymin;
	Lz=zmax-zmin;
	
	//------------------------------------------------------------------
	// periodic wrapping
	
	for(i=0;i<particles;++i)
	{
		if(x[i]<xmin){x[i]=x[i]+Lx;}
		if(x[i]>xmax){x[i]=x[i]-Lx;}
		if(y[i]<ymin){y[i]=y[i]+Ly;}
		if(y[i]>ymax){y[i]=y[i]-Ly;}
		
		
	}
	
	//------------------------------------------------------------------
	// exports
	
	sprintf(file_path,"%s/Si.dat",current_folder);
	fp=fopen(file_path,"w+");

	fprintf(fp,"%d\nSi_diamond\n",particles);
	for(i=0;i<particles;++i)
	{
		//if(i>=particles-xpart*ypart){fprintf(fp,"F\t%lf\t%lf\t%lf\n",x[i],y[i],z[i]);}else{fprintf(fp,"Si\t%lf\t%lf\t%lf\n",x[i],y[i],z[i]);}
		fprintf(fp,"Si\t%lf\t%lf\t%lf\n",x[i],y[i],z[i]);
	}
	
	
	fclose(fp);//printf("Generated %s\n",file_path);
	
	printf("title\tSi_diamond\n");
	printf("particles\t%d\n",particles);
	printf("supercell\n");
	printf("xmin\t%lf\nxmax\t%lf\nymin\t%lf\nymax\t%lf\nzmin\t%lf\nzmax\t%lf\n",xmin,xmax,ymin,ymax,zmin,zmax);
	printf("species\tx\ty\tz\tmobileFlag\tthermoFlag\n");
	for(i=0;i<particles;++i)
	{
		if(i<fixed*xpart*ypart){
		printf("Si\t%lf\t%lf\t%lf\t%d\t%d\n",x[i],y[i],z[i],0,0);}
		else
		{printf("Si\t%lf\t%lf\t%lf\t%d\t%d\n",x[i],y[i],z[i],1,1);}
	}
	
	//------------------------------------------------------------------
	
	free(xFCC);free(yFCC);free(zFCC);
	free(x);free(y);free(z);
	return 0;
}
int equal(double A,double B);
int equal(double A,double B)
{
	int value;
	if(fabs(A-B)<equal_diff){value=1;}else{value=0;}
	return value;
}



		
