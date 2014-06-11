#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "field.h"

//#ifndef NOFFTWLIB
#include <fftw3.h>
#include "fft_diff.h"
//#endif

#ifndef NOFFTWLIB
#define USE_FFTW 1
#else
#define USE_FFTW 0
#endif

#define MOD(A,B) ((A>=0)?(A%B):((A*(1-B)))%B)

double IndGen_3D(int i,int j,int k,int n)
{
    double a,b,c;

    a=(i<=n/2)?i:n-i;
    b=(j<=n/2)?j:n-j;
    c=(k<=n/2)?k:n-k;

    return (sqrt(a*a+b*b+c*c)/n);
    
}

double IndGen_3D_bis(int i, int j, int k, int n)
{
    double a, b, c, kk;
    

    a=(i<=n/2)?i:n-i;
    b=(j<=n/2)?j:n-j;
    c=(k<=n/2)?k:n-k;

    kk = sqrt(sin(0.5*PI*a/n)*sin(0.5*PI*a/n) +
	      sin(0.5*PI*b/n)*sin(0.5*PI*b/n) +
	      sin(0.5*PI*c/n)*sin(0.5*PI*c/n));
    return kk;
}

//compute the potential from the density
int ComputePotential(density_grid *grid)
{

    
    int i,j,k,l,m;

    size_t Nalloc = 2*(grid->Nx/2+1)*grid->Ny*grid->Nz;

    FLOAT *dummy=NULL;
    FLOAT *mygrid;
    

    FFTW_COMPLEX *cpl_grid;
    FFTW_PLAN p1,p3;



    //Resort data for inplace fft ...
    fprintf(stdout, "Nalloc %zd \n", Nalloc);
    fflush(stdout);
    grid->grid = realloc(grid->grid,Nalloc*sizeof(FLOAT));

    dummy = calloc(Nalloc, sizeof(FLOAT));

    mygrid = grid->grid;

    for (i=0;i<Nalloc;i++){
      dummy[i]=mygrid[i];
    }

    puts("E");
    //Padds ...
    for (k=0,l=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,l++,m++)
	    {
		mygrid[m]=dummy[l];
	    }
    

    free(dummy);


#ifdef FFTW_NTHREADS
    printf ("Using fftw3 with %d threads.\n",FFTW_NTHREADS);
    FFTWNAME(init_threads)();
    FFTWNAME(plan_with_nthreads)(FFTW_NTHREADS);
#endif

    cpl_grid = (FFTW_COMPLEX *) mygrid;



    p1=FFTWNAME(plan_dft_r2c_3d)(grid->Nz,grid->Ny,grid->Nx,mygrid,(FFTW_COMPLEX *)mygrid,FFTW_ESTIMATE);
    FFTWNAME(execute)(p1);
    puts("computing potential, finished first FFT");


    

    for (k=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++)
	    for (i=0;i<grid->Nx/2+1;i++,m++)
	    {
	      double k2;
	      k2 = pow(IndGen_3D(i,j,k,grid->Nx),2);
	      cpl_grid[m] *= (-1.0/k2)/grid->NNodes;
	    }

    cpl_grid[0] = 0.0;
  
    p3=FFTWNAME(plan_dft_c2r_3d)(grid->Nz,grid->Ny,grid->Nx,(FFTW_COMPLEX *)mygrid,mygrid,FFTW_ESTIMATE);
    FFTWNAME(execute)(p3);



#ifdef FFTW_NTHREADS
    FFTWNAME(cleanup_threads)();
#endif


    for (k=0,l=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,l++,m++)
	    {
		mygrid[l]=mygrid[m]/(grid->NNodes);
	    }
    


    mygrid = realloc(mygrid,grid->NNodes*sizeof(FLOAT));


  
    return 0;
}

//#endif

//compute the gradiant of the density
int ComputeDensityGradiant(density_grid *density,int mode)
{
    unsigned int i,j,k,n;
    int dx,dy,dz;
    long long n_to_alloc;
    FLOAT *grid;
    FLOAT *tab;
    //grid = density->grid;

    
    if (density->HasGradiant) return -1;
    dx=1;
    dy=density->Nx;
    dz=density->Nx*density->Ny;

    grid=density->grid;
    density->grad= (FLOAT *) realloc(density->grad,3*sizeof(FLOAT)*density->NNodes);

//#if (USE_FFTW==1)
    if (mode == USE_FFT)
      {

	n_to_alloc = sizeof(FLOAT) * density->NNodes;
	tab = (FLOAT *) malloc (n_to_alloc);
	
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_DX,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) density->grad[3*i]=tab[i];
	
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_DY,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) density->grad[3*i+1]=tab[i];
	
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_DZ,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) density->grad[3*i+2]=tab[i];
	
	free(tab);
	density->HasGradiant=1;
	
      }

//#endif
    else
    {
	for (k=0,n=0;k<density->Nz;k++)
	    for (j=0;j<density->Ny;j++)
		for (i=0;i<density->Nx;i++,n++)
		{
		    //if (grid[density->maxindex]<grid[n]) density->maxindex=n;
		    if (i==0) {density->grad[3*n]=(grid[n+dx]-grid[n+(density->Nx-2)*dx])/2.;}
		    else if (i==density->Nx-1) {density->grad[3*n]=(grid[n-i*dx+dx]-grid[n-dx])/2.;}
		    else density->grad[3*n]=(grid[n+dx]-grid[n-dx])/2.;
		    
		    if (j==0) {density->grad[3*n+1]=(grid[n+dy]-grid[n+(density->Ny-2)*dy])/2.;}
		    else if (j==density->Ny-1) {density->grad[3*n+1]=(grid[n-j*dy+dy]-grid[n-dy])/2.;}
		    else density->grad[3*n+1]=(grid[n+dy]-grid[n-dy])/2.;
		    
		    if (k==0) {density->grad[3*n+2]=(grid[n+dz]-grid[n+(density->Nz-2)*dz])/2.;}
		    else if (k==density->Nz-1) {density->grad[3*n+2]=(grid[n-k*dz+dz]-grid[n-dz])/2.;}
		    else density->grad[3*n+2]=(grid[n+dz]-grid[n-dz])/2.;
		    //printf ("_%lg_ _%lg_ _%lg_\n",density->grad[3*n],density->grad[3*n+1],density->grad[3*n+2]);
		}
    }
    density->HasGradiant=1;
    return 1;
}


//compute the hassian of the density
int ComputeDensityHessian(density_grid *density,int mode)
{
    unsigned int i,j,k,n;
    int dx,dy,dz;

    int dxp,dyp,dzp;
    int dxm,dym,dzm;

    long long n_to_alloc;
    FLOAT *hess;

    FLOAT *grid;
    FLOAT *tab;

    grid=density->grid;
    if (density->HasHessian) return -1;
//    if (!density->HasGradiant) ComputeDensityGradiant(density,mode);
    
//    grad = density->grad;

    dx=1; //3 for second  method, 1 for first ...
    dy=dx*density->Nx;
    dz=dy*density->Ny;
    
    density->hessian= (FLOAT *) realloc(density->hessian,6*sizeof(FLOAT)*density->NNodes);
    hess = density->hessian;
    //density->maxindex=0;
  
    //#if (USE_FFTW==1)
    if (mode == USE_FFT)
    {
      n_to_alloc = sizeof(FLOAT)*density->NNodes;
	tab = (FLOAT *) malloc (n_to_alloc);
	
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	printf("starts first FFT\n");
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_D2X,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) hess[6*i]=tab[i];
	printf("finishes first FFT\n");

	printf("starts second FFT\n");
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_D2Y,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) hess[6*i+3]=tab[i];
	
	printf("starts 3 FFT\n");
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_D2Z,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) hess[6*i+5]=tab[i];

	printf("starts 4 FFT\n");
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_DX|FLAG_DY,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) hess[6*i+1]=tab[i];

	printf("starts 5 FFT\n");
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_DX|FLAG_DZ,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) hess[6*i+2]=tab[i];

	printf("starts 6 FFT\n");	
	//for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) tab[i]=grid[i];
	FFT_diff(grid,density->Nx,density->Ny,density->Nz,1.,1.,1.,FLAG_DY|FLAG_DZ,&tab);
	for (i=0,n=density->NNodes-1;i<density->NNodes;i++,n--) hess[6*i+4]=tab[i];
	
	free(tab);
	density->HasHessian=1;
	//return 1;
    }
    //#endif
    else
    {
	//Compute the hessian
	for (k=0,n=0;k<density->Nz;k++)
	    for (j=0;j<density->Ny;j++)
		for (i=0;i<density->Nx;i++,n++)
		{
		    dxp=dx;dxm=-dx;
		    dyp=dy;dym=-dy;
		    dzp=dz;dzm=-dz;
		    
		    if (i==0) dxm=dy-dx;
		    if (i==density->Nx-1) dxp=-dy+dx;
		    if (j==0) dym=dz-dy;
		    if (j==density->Ny-1) dyp=-dz+dy;
		    if (k==0) dzm=density->NNodes-dz;
		    if (k==density->Nz-1) dzp=-density->NNodes+dz;
		    
		    //xx, yy and zz
		    hess[6*n] = grid[n+dxp]+grid[n+dxm]-2*grid[n];
		    hess[6*n+3] = grid[n+dyp]+grid[n+dym]-2*grid[n];
		    hess[6*n+5] = grid[n+dzp]+grid[n+dzm]-2*grid[n];
		    
		    //xy, xz and yz
		    hess[6*n+1]= 0.5*(grid[n+dxp+dyp]+grid[n+dxm+dym]+2*grid[n]
				      -grid[n+dxp]-grid[n+dxm]-grid[n+dyp]-grid[n+dym]);
		    hess[6*n+2]= 0.5*(grid[n+dxp+dzp]+grid[n+dxm+dzm]+2*grid[n]
				      -grid[n+dxp]-grid[n+dxm]-grid[n+dzp]-grid[n+dzm]);
		    hess[6*n+4]= 0.5*(grid[n+dyp+dzp]+grid[n+dym+dzm]+2*grid[n]
				      -grid[n+dyp]-grid[n+dym]-grid[n+dzp]-grid[n+dzm]);
		}
    }
    /*	
    for (i=0;i<6*density->NNodes;i+=6)
    {
	hess[i]/=(density->dx*density->dx);
	hess[i+1]/=(density->dx*density->dy);
	hess[i+2]/=(density->dx*density->dz);
	hess[i+3]/=(density->dy*density->dy);
	hess[i+4]/=(density->dy*density->dz);
	hess[i+5]/=(density->dz*density->dz);
    }
    */
    //printf("max=%f\n",grid[density->maxindex]);
    density->HasHessian=1;
    return 1;
}

int ComputeVelocityTorsion(density_grid *density)
{
  //    unsigned int i,j,k,n;
    long long i,j,k,n;
    int dx,dy,dz;
    
    int dxp,dyp,dzp;
    int dxm,dym,dzm;
    int ll;

    FLOAT *torsion;
    FLOAT *grid;
    FLOAT *vel_x;
    FLOAT *vel_y;
    FLOAT *vel_z;
    FLOAT *dens;
    FLOAT *tab;
    FLOAT * my_grid;
    
    grid  = density->grid;
    vel_x = density->grid_vx;
    vel_y = density->grid_vy;
    vel_z = density->grid_vz;

    dx=1; //3 for second  method, 1 for first ...
    dy=dx*density->Nx;
    dz=dy*density->Ny;

    fprintf(stdout, "Starting to compute torsion on %d grid size\n", density->Nx);

    density->torsion= (FLOAT *) realloc(density->torsion,6*sizeof(FLOAT)*density->NNodes);
    torsion = density->torsion;

    fprintf(stdout, "Starting to compute torsion on %d grid size\n", density->Nx);

    /*    
    for (k=0,n=0;k<density->Nz;k++){
      for (j=0;j<density->Ny;j++){
	for (i=0;i<density->Nx;i++,n++){
	  if(grid[n]>0.001){	    
	    vel_x[n] = vel_x[n]/grid[n];
	    vel_y[n] = vel_y[n]/grid[n];
	    vel_z[n] = vel_z[n]/grid[n];
	  }
	}
      }
    }
    */

    fprintf(stdout, "Finished the initialization - Starting to compute torsion on %d grid size\n", density->Nx);    
    //Compute the torsion
    for (k=0,n=0;k<density->Nz;k++){
      for (j=0;j<density->Ny;j++){
	for (i=0;i<density->Nx;i++,n++){
	  //	  fprintf(stdout, "%d\n",n);
	    dxp=dx;dxm=-dx;
	    dyp=dy;dym=-dy;
	    dzp=dz;dzm=-dz;
	    
	    if (i==0) dxm=dy-dx;
	    if (i==density->Nx-1) dxp=-dy+dx;
	    if (j==0) dym=dz-dy;
	    if (j==density->Ny-1) dyp=-dz+dy;
	    if (k==0) dzm=density->NNodes-dz;
	    if (k==density->Nz-1) dzp=-density->NNodes+dz;
	    
	    //xx, yy and zz
	    torsion[6*n] = 0.5*(vel_x[n+dxp] - vel_x[n+dxm]);
	    torsion[6*n+3] = 0.5*(vel_y[n+dyp] - vel_y[n+dym]);
	    torsion[6*n+5] = 0.5*(vel_z[n+dzp] - vel_z[n+dzm]);
	    
	    //xy, xz and yz
	    torsion[6*n+1]= 0.25*(vel_x[n+dyp] - vel_x[n+dym]
				  + vel_y[n+dxp] - vel_y[n+dxm]);

	    torsion[6*n+2]= 0.25*(vel_x[n+dzp] - vel_x[n+dzm]
				  + vel_z[n+dxp] - vel_z[n+dxm]);

	    torsion[6*n+4]= 0.25*(vel_y[n+dzp] - vel_y[n+dzm]
				  + vel_z[n+dyp] - vel_z[n+dym]);
	    /*
	    for(ll=0;ll<6;ll++){
	      if(torsion[6*n+ll]>1.0e5){
		torsion[6*n+ll] = 1.0e5;
		fprintf(stdout, "ARRRR %d = %d %d %d\n ",ll,i,j,k);
		fprintf(stdout, "ARRR XZ %e %e \n", vel_x[n+dzp], vel_x[n+dzm]);
		fprintf(stdout, "ARRR XY %e %e \n", vel_x[n+dyp], vel_x[n+dym]);
		fprintf(stdout, "ARRR YX %e %e \n", vel_y[n+dxp], vel_y[n+dxm]);
		fprintf(stdout, "ARRR YZ %e %e \n", vel_y[n+dzp], vel_y[n+dzm]);
	      }
	      if(torsion[6*n+ll]<-1.0e5){
		torsion[6*n+ll] = -1.0e5;
	      }

	    }
	    */
	  }
      }
}

    density->HasTorsion=1;
    fprintf(stdout, "Finished to compute torsion on %d grid size\n", density->Nx);    
    return 1;
}


//transforms a density_field into a density_grid
//tables is *NOT* reallocated, so don't free density_f ... 
int DensityField2DensityGrid(density_field density_f,density_grid *density_g)
{
    int i;

    density_g->Nx=density_f.Nx;
    density_g->Ny=density_f.Ny;
    density_g->Nz=density_f.Nz;
    density_g->NNodes=density_f.NNodes;

    density_g->x0=-density_f.header.BoxSize/2.;
    density_g->y0=-density_f.header.BoxSize/2.;
    density_g->z0=-density_f.header.BoxSize/2.;
    density_g->dx=density_f.delta;
    density_g->dy=density_f.delta;
    density_g->dz=density_f.delta;

    density_g->grid = (FLOAT *) malloc (sizeof(FLOAT)*density_g->NNodes);
    
    for (i=0;i<density_g->NNodes;i++) 
	density_g->grid[i]=(FLOAT)density_f.grid[i];

    free(density_f.grid);

    return 0;
    
}

int Compute_CIC_Density(float *Pos,float *Val,int ValDim,int ValPad,long long npart,density_grid *grid,
			int Nx,int Ny,int Nz,
			float x0,float y0,float z0,
			float xs,float ys,float zs)
{
    long long i,j;

    int xp,yp,zp;
    int dx,dy,dz;

    int n[8];
    double p[8];
    int xdir,ydir,zdir;
    float lx,ly,lz,tmp;
    float vol;
    float point_val = 1.0;
    float min_val, max_val;
    float min_val_x, max_val_x, delta_x;
    float min_val_y, max_val_y, delta_y;
    float min_val_z, max_val_z, delta_z;
    min_val = 1e10;
    max_val = -1e10;
    min_val = 1e10;    max_val = -1e10;  
    min_val = 1e10;    max_val = -1e10; 
    min_val = 1e10;    max_val = -1e10;



    printf ("CIC xs= %f\n",xs);

    fprintf (stdout, "Computing CIC density\n");fflush(stdout);
    if(Val!=NULL){
      fprintf(stdout, "Using Weights in the Interpolation\n");
    }

    grid->Nx=Nx;
    grid->Ny=Ny;
    grid->Nz=Nz;
    
    grid->dx=xs/Nx;
    grid->dy=ys/Ny;
    grid->dz=zs/Nz;

    grid->x0=x0;//+grid->dx/2.;
    grid->y0=y0;//+grid->dy/2.;
    grid->z0=z0;//+grid->dz/2.;

    grid->NNodes = Nx*Ny*Nz;

    grid->HasGradiant=grid->HasHessian=0;
    grid->maxindex=0;

    /*compute the miniumn and maximum of the particle positions*/
    for (i=0;i<((long long)(3*npart));i+=3)
    {

	if(Pos[i]>max_val_x)	    max_val_x = Pos[i];
	if(Pos[i]<min_val_x)	    min_val_x = Pos[i];
	if(Pos[i+1]>max_val_y)	    max_val_y = Pos[i+1];
	if(Pos[i+1]<min_val_y)	    min_val_y = Pos[i+1];
	if(Pos[i+2]>max_val_z)	    max_val_z = Pos[i+2];
	if(Pos[i+2]<min_val_z)	    min_val_z = Pos[i+2];
    }
    fprintf(stdout, "min max (x_pos) %g %g \n", min_val_x, max_val_x);
    fprintf(stdout, "min max (y_pos) %g %g \n", min_val_y, max_val_y);
    fprintf(stdout, "min max (z_pos) %g %g \n", min_val_z, max_val_z);
    delta_x = max_val_x - min_val_x;
    delta_y = max_val_y - min_val_y;
    delta_z = max_val_z - min_val_z; 


    //Now compute the CIC stuff ...
    if (grid->grid!=NULL) free(grid->grid);
    grid->grid = calloc (grid->NNodes,sizeof(FLOAT));
    printf ("(%d particles) ... \n",npart);

    dx = 1;
    dy = dx*grid->Nx;
    dz = dy*grid->Ny;

    vol = (grid->dx*grid->dy*grid->dz);
    printf ("x0=%f y0=%f z0=%f ds=%f vol=%f npar*3 %lld \n",grid->x0, grid->y0, grid->z0, grid->dx, vol, 3*(long long)npart);

    for (i=0;i<grid->NNodes;i++){
	grid->grid[i]=0;
    }
    for (i=0;i<((long long)(3*npart));i+=3)
    {
      if (Val!=NULL) {
	    point_val =Val[ValDim*(i/3)+ValPad];
	    //	    fprintf(stdout, "VAL %d %d %f\n",i, ValDim*(i/3)+ValPad,point_val);
      }
	if (Pos[i]>=x0+xs) Pos[i]-=xs;
	if (Pos[i+1]>=y0+ys) Pos[i+1]-=ys;
	if (Pos[i+2]>=z0+zs) Pos[i+2]-=zs;

	if (Pos[i]<x0) Pos[i]+=xs;
	if (Pos[i+1]<y0) Pos[i+1]+=ys;
	if (Pos[i+2]<z0) Pos[i+2]+=zs;
	
	xp = (Pos[i]-grid->x0)/grid->dx;    
	yp = (Pos[i+1]-grid->y0)/grid->dy;
	zp = (Pos[i+2]-grid->z0)/grid->dz;

	if (xp>=grid->Nx) xp=0;
	if (yp>=grid->Ny) yp=0;
	if (zp>=grid->Nz) zp=0;

	if (xp<0) xp=grid->Nx-1;
	if (yp<0) yp=grid->Ny-1;
	if (zp<0) zp=grid->Nz-1;

	
	lx=Pos[i]-xp*grid->dx-grid->x0;
	tmp=grid->x0+(xp+1)*grid->dx-Pos[i];
	if (lx<tmp) 
	{
	    xdir = -1;
	    lx+=grid->dx/2;
	}
	else 
	{
	    xdir = 1;
	    lx = tmp + grid->dx/2;
	}

	ly=Pos[i+1]-yp*grid->dy-grid->y0;
	tmp=grid->y0+(yp+1)*grid->dy-Pos[i+1];
	if (ly<tmp) 
	{
	    ydir = -1;
	    ly+=grid->dy/2;
	}
	else 
	{
	    ydir = 1;
	    ly = tmp + grid->dy/2;
	}

	lz=Pos[i+2]-zp*grid->dz-grid->z0;
	tmp=grid->z0+(zp+1)*grid->dz-Pos[i+2];
	if (lz<tmp) 
	{
	    zdir = -1;
	    lz+=grid->dz/2;
	}
	else 
	{
	    zdir = 1;
	    lz = tmp + grid->dz/2;
	}


	n[0]=dx*xp + dy*yp +dz*zp;
	n[1]=dx*xp + dy*yp +dz*MOD((zp+zdir),(grid->Nz-1));
	n[2]=dx*xp + dy*MOD((yp+ydir),(grid->Ny-1)) +dz*zp;
	n[3]=dx*xp + dy*MOD((yp+ydir),(grid->Ny-1)) +dz*MOD((zp+zdir),(grid->Nz-1));
	n[4]=dx*MOD((xp+xdir),(grid->Nx-1)) + dy*yp +dz*zp;
	n[5]=dx*MOD((xp+xdir),(grid->Nx-1)) + dy*yp +dz*MOD((zp+zdir),(grid->Nz-1));
	n[6]=dx*MOD((xp+xdir),(grid->Nx-1)) + dy*MOD((yp+ydir),(grid->Ny-1)) +dz*zp;
	n[7]=dx*MOD((xp+xdir),(grid->Nx-1)) + dy*MOD((yp+ydir),(grid->Ny-1)) +dz*MOD((zp+zdir),(grid->Nz-1));


	p[0]=(lx*ly*lz)/vol;
	p[1]=(lx*ly*(grid->dz-lz))/vol;
	p[2]=(lx*(grid->dy-ly)*lz)/vol;
	p[3]=(lx*(grid->dy-ly)*(grid->dz-lz))/vol;
	p[4]=((grid->dx-lx)*ly*lz)/vol;
	p[5]=((grid->dx-lx)*ly*(grid->dz-lz))/vol;
	p[6]=((grid->dx-lx)*(grid->dy-ly)*lz)/vol;
	p[7]=((grid->dx-lx)*(grid->dy-ly)*(grid->dz-lz))/vol;

	
	if (Val==NULL){
	  for (j=0;j<8;j++) grid->grid[n[j]] += p[j];
	}else{
	  for (j=0;j<8;j++) grid->grid[n[j]] += point_val * p[j];
	}
/*
	if(p[0]>1.0){
	    fprintf(stdout, "%d -> %f %f %f [%f]\n", i, Pos[i], Pos[i+1], Pos[i+2], p[0]);
	    exit(1);
	}
*/
	//for (j=0;j<8;j++) printf ("%f ",p[j]);
	//for (j=0;j<8;j++) printf ("%d ",n[j]);printf("\n");
	//printf ("\nP=%f\n",p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7]);
	//printf ("%d %d \n\n",MOD((-1),3),2%2);
    }
    //for (i=0;i<npart;i++) grid->grid[i]=log(grid->grid[i]);


    grid->x0+=grid->dx/2.0;
    grid->y0+=grid->dy/2.0;
    grid->z0+=grid->dz/2.0;

    /*get the minimum and maximum*/
    /*
    for(i=0;i<grid->NNodes;i++){
	if(grid->grid[i]<min_val)
	    min_val = grid->grid[i];
	if(grid->grid[i]>max_val)
	    max_val = grid->grid[i];
    }
    fprintf(stdout, "min max in CIC interpolation : %g %g\n", min_val, max_val);

    float mean_value=0.0;
    for (i=0;i<grid->NNodes;i++)
    {
	mean_value += grid->grid[i];
    }
    mean_value = mean_value/grid->NNodes;

    for (i=0;i<grid->NNodes;i++)
    {
	grid->grid[i] = (grid->grid[i]/mean_value) - 1.0;
    }
    */
    printf ("done.\n");
    return 0;
}


//#if (USE_FFTW==1)

int GaussSmooth(density_grid *grid,float sigma_p, int item)
{

    
    //float sigma = sigma_p*grid->dx;
    
    int i,j,k,l,m;
    long long n_nodes;
    size_t Nalloc = 2*(grid->Nx/2+1)*grid->Ny*grid->Nz;


    FLOAT *dummy;
    FLOAT *mygrid;
    

    FFTW_COMPLEX *cpl_grid;
    FFTW_PLAN p1,p3;
    mygrid=NULL;
    dummy=NULL;
    printf ("Smoothing with sigma=%f ... ",sigma_p);fflush(0);

    //Resort data for inplace fft ...
    if(item!=DENS){
      if(item==VEL_X){grid->grid_vx = realloc(grid->grid_vx,Nalloc*sizeof(FLOAT));}
      if(item==VEL_Y){grid->grid_vy = realloc(grid->grid_vy,Nalloc*sizeof(FLOAT));}
      if(item==VEL_Z){grid->grid_vz = realloc(grid->grid_vz,Nalloc*sizeof(FLOAT));}
    }else{
      grid->grid = realloc(grid->grid,Nalloc*sizeof(FLOAT));
    }
    dummy = calloc (Nalloc,sizeof(FLOAT));
    fprintf(stderr, "N to alloc %d\n", Nalloc);
    /*
    if (sizeof(FLOAT) != sizeof(double))
      {
	mygrid = calloc(Nalloc,sizeof(FLOAT));
	for (i=0;i<Nalloc;i++)
	  mygrid[i]=grid->grid[i];
	free (grid->grid);
      }
    else
*/
    if(item!=DENS){
      if(item==VEL_X){mygrid = grid->grid_vx;fprintf(stdout, "doing VELX\n");}
      if(item==VEL_Y){mygrid = grid->grid_vy;fprintf(stdout, "doing VELY\n");}
      if(item==VEL_Z){mygrid = grid->grid_vz;fprintf(stdout, "doing VELZ\n");}
    }else{
      mygrid = grid->grid;
    }
    //memcpy(dummy,grid->grid,Nalloc*sizeof(double));

    for (i=0;i<Nalloc;i++)
      dummy[i]=mygrid[i];

    n_nodes = grid->Nz * grid->Ny * grid ->Nx;
    //Padds ...
    for (k=0,l=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,l++,m++)
	    {
		mygrid[m]=dummy[l];
	    }
    
    free(dummy);

#ifdef FFTW_NTHREADS
    //printf ("Using fftw3 with %d threads.\n",FFTW_NTHREADS);
    FFTWNAME(init_threads)();
    FFTWNAME(plan_with_nthreads)(FFTW_NTHREADS);
#endif

    cpl_grid = (FFTW_COMPLEX *) mygrid;

    printf ("\rSmoothing with sigma=%f ... (grid FFT)",sigma_p);fflush(0);



    p1=FFTWNAME(plan_dft_r2c_3d)(grid->Nz,grid->Ny,grid->Nx,mygrid,(FFTW_COMPLEX *)mygrid,FFTW_ESTIMATE);
    FFTWNAME(execute)(p1);
    
    printf ("\rSmoothing with sigma=%f ... (convolving)",sigma_p);fflush(0);

    for (k=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++)
	    for (i=0;i<grid->Nx/2+1;i++,m++)
	    {
	      double k2;

	      k2 = pow(IndGen_3D(i,j,k,grid->Nx),2);
	      //	      cpl_grid[m] *= sqrt(PI*2*sigma_p*sigma_p)*exp(-PI*PI*k2*2*sigma_p*sigma_p)/grid->NNodes;
	      cpl_grid[m] *= sqrt(PI*2*sigma_p*sigma_p)*exp(-PI*PI*k2*2*sigma_p*sigma_p)/(1.0*n_nodes);
	    }


    printf ("\rSmoothing with sigma=%f ... (inverse FFT)",sigma_p);fflush(0);


    p3=FFTWNAME(plan_dft_c2r_3d)(grid->Nz,grid->Ny,grid->Nx,(FFTW_COMPLEX *)mygrid,mygrid,FFTW_ESTIMATE);
    FFTWNAME(execute)(p3);
    
#ifdef FFTW_NTHREADS
    FFTWNAME(cleanup_threads)();
#endif


    for (k=0,l=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,l++,m++)
	    {
	      //		mygrid[l]=mygrid[m]/(grid->NNodes);
		mygrid[l]=mygrid[m]/(1.0*n_nodes);
	    }
    
    mygrid = realloc(mygrid,grid->NNodes*sizeof(FLOAT));


    printf ("\rSmoothing with sigma=%f(%f) ... done.                     \n",sigma_p,sigma_p*grid->dx);
  
    return 0;
}

int FilterSmooth(density_grid *grid,double (*filterfunc)(double,double,double,void *),const char *filtername, void *parms)
{
    FLOAT *filter;
    
    //float sigma = sigma_p*grid->dx;
    
    int i,j,k,l,m;
    double u,v,w;
    size_t Nalloc = 2*(grid->Nx/2+1)*grid->Ny*grid->Nz;
    double mean;
    FLOAT *dummy;
    FLOAT *mygrid;
    
    FFTW_COMPLEX *cpl_filter;
    FFTW_COMPLEX *cpl_grid;
    FFTW_PLAN p1,p2,p3;

    printf ("Filtering with '%s' ... ",filtername);fflush(0);

    //Resort data for inplace fft ...
    grid->grid = realloc(grid->grid,Nalloc*sizeof(FLOAT));
    dummy = calloc (Nalloc,sizeof(FLOAT));
/*
    if (sizeof(FLOAT) != sizeof(double))
      {
	mygrid = calloc(Nalloc,sizeof(FLOAT));
	for (i=0;i<Nalloc;i++)
	  mygrid[i]=grid->grid[i];
	free (grid->grid);
      }
    else
*/
      mygrid = grid->grid;
	  
    //memcpy(dummy,grid->grid,Nalloc*sizeof(double));

    for (i=0;i<Nalloc;i++)
      dummy[i]=mygrid[i];

    //Padds ...
    for (k=0,l=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,l++,m++)
	    {
		mygrid[m]=dummy[l];
	    }
    
    free(dummy);

#ifdef FFTW_NTHREADS
    printf ("Using fftw3 with %d threads.\n",FFTW_NTHREADS);
    FFTWNAME(init_threads)();
    FFTWNAME(plan_with_nthreads)(FFTW_NTHREADS);
#endif

    cpl_grid = (FFTW_COMPLEX *) mygrid;

    printf ("\rFiltering with '%s' ... (grid FFT)",filtername);fflush(0);

    p1=FFTWNAME(plan_dft_r2c_3d)(grid->Nz,grid->Ny,grid->Nx,mygrid,(FFTW_COMPLEX *)mygrid,FFTW_ESTIMATE);
    FFTWNAME(execute)(p1);
    

    printf ("\rFiltering with '%s' ... (filter FFT)",filtername);fflush(0);
    filter = malloc (Nalloc*sizeof(FLOAT));
    mean=0;
    for (k=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,m++)
	    {
		u=(double)grid->Nx/2 - fabs(i - (double)grid->Nx/2);//*grid->dx;
		v=(double)grid->Ny/2 - fabs(j - (double)grid->Ny/2);//*grid->dy;
		w=(double)grid->Nz/2 - fabs(k - (double)grid->Nz/2);//*grid->dz;
		filter[m]= filterfunc(u,v,w,parms);
		//Norm * exp (-1/2.*(u*u+v*v+w*w)/(sigma_p*sigma_p));
		
		//printf ("%d %d %d ||| %f %f %f:      %le\n",i,j,k,u,v,w,filter[m]);
		mean +=  filter[m];
	    }

    p2=FFTWNAME(plan_dft_r2c_3d)(grid->Nz,grid->Ny,grid->Nx,filter,(FFTW_COMPLEX *)filter,FFTW_ESTIMATE);
    FFTWNAME(execute)(p2);

    //printf ("mean=%le\n",mean);
    cpl_filter = (FFTW_COMPLEX *) filter;
      
    printf ("\rFiltering with '%s' ... (convolving)",filtername);fflush(0);
    
    for (i=0;i<Nalloc/2;i++)
    {
	//printf ("%le+i*%le\n",creal(cpl_filter[i]),cimag(cpl_filter[i]));
	cpl_grid[i] *= (cpl_filter[i]/grid->NNodes);
    }
    free(filter);

    p3=FFTWNAME(plan_dft_c2r_3d)(grid->Nz,grid->Ny,grid->Nx,(FFTW_COMPLEX *)mygrid,mygrid,FFTW_ESTIMATE);
    FFTWNAME(execute)(p3);
    
#ifdef FFTW_NTHREADS
    FFTWNAME(cleanup_threads)();
#endif


    for (k=0,l=0,m=0;k<grid->Nz;k++)
	for (j=0;j<grid->Ny;j++,m+=2)
	    for (i=0;i<grid->Nx;i++,l++,m++)
	    {
		mygrid[l]=mygrid[m]/(grid->NNodes);
	    }
    
    mygrid = realloc(mygrid,grid->NNodes*sizeof(FLOAT));

    printf ("\rFiltering with '%s' ... done.                     \n",filtername);
  
    return 0;
}

//#endif
