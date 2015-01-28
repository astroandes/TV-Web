#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "fftw3.h"
#include "fft_diff.h"

FFTW_COMPLEX Ipow(int n)
{
  int i;

  FFTW_COMPLEX a;

  i=abs(n)&3;
  if (i==0) a = +1.+0.*I;
  else if (i==1) a = ((n>0)?(+0.+1.*I):(+0.-1.*I));
  else if (i==2) a = -1.+0.*I;
  else if (i==3) a = ((n>0)?(+0.-1.*I):(+0.+1.*I));

  return a;
}

//generate properly ordered index for wave vector
FFTW_COMPLEX *fft_indgen(int n,double factor_p,int order)
{
    int j;
    double factor;
    FFTW_COMPLEX *val;
    FFTW_COMPLEX c;

    if (factor_p==0) 
      factor = 1;
    else
      factor=factor_p;

    val = (FFTW_COMPLEX *) FFTWNAME(malloc) (n*sizeof(FFTW_COMPLEX));
    
    if ((order==0)||(n==1))
      {
	for (j=0;j<n;j++) val[j]=1.;
      }
    else
      {
	c=Ipow(order)*pow((double)2.*PI/((double)n*factor),order);
	for (j=0;j<n;j++)
	  {
	   
	    if (j>n/2)
		val[j]=c*pow((double)(j-n),order);
	    else
		val[j]=c*pow((double)j,order);
	    //printf ("%f\n",val[i]);
	  }

	if (order<0) val[0]=0;//1.E12+0*val[1]/pow((double)n,-order);
      }

    //printf ("pow=%lg\n",pow(-1.5,-2));

    return val;
}

//tab_in must be periodic !!!

//f is for floats ...
//nx,ny and nz are number of elements along x/y/z
//fx,fy,fz are the variables value as a function of x,y,z (x'=fx*x ...)
//this is used to normalize properly 
//tab_out must be allocated before calling (or NULL)
int FFT_diff(FLOAT *tab_in,int nx,int ny,int nz,double fx,double fy,double fz,int flags,FLOAT **tab_out)
{
  FFTW_COMPLEX *kx;
  FFTW_COMPLEX *ky;
  FFTW_COMPLEX *kz;
  FFTW_COMPLEX *k;

  FFTW_PLAN p;
  FFTW_COMPLEX *out;

  int nout=0;
  int nxout,nyout,nzout;

  int x,y,z,n;
  int order;
  FLOAT *mytab;

  //FFTW_COMPLEX *tmp;

  //printf ("%d %d %d\n",nx,ny,nz);

  //tmp = (double *) malloc(sizeof(double)*nx*ny*nz);
  
  if ((nx<1)||(ny<1)||(nz<1)) 
    {
      fprintf (stderr,"Error in FFT_diff_d(), nx,ny and nz must be >=1\n");
      return -1;
    }
  /*
  if (nz!=1) {nout = nx*ny*(nz/2+1);nzout=nz/2+1;}
  else if (ny!=1) {nout = nx*(ny/2+1);nyout=ny/2+1;}
  else {nxout=nout=nx/2+1;}
  */
   
  //Could improve all that ...

  nout=nx*ny*(nz);
  nxout=nx;
  nyout=ny;
  nzout=nz;


#ifdef FFTW_NTHREADS
  FFTWNAME(init_threads)();
  FFTWNAME(plan_with_nthreads)(FFTW_NTHREADS);
#endif
  
  //printf ("nout=%d\n",nout);exit(0);

  out = FFTWNAME(malloc)(sizeof(FFTW_COMPLEX)*nout);

  for (x=0;x<nout;x++) out[x]=tab_in[x];

  if (*tab_out==NULL) 
    {
      *tab_out = (FLOAT *) malloc(sizeof(FLOAT)*nx*ny*nz);
      mytab=*tab_out;
    }
  else mytab=*tab_out;

  k=(FFTW_COMPLEX *)FFTWNAME(malloc) (sizeof(FFTW_COMPLEX)*nout);
  for (n=0;n<nout;n++) k[n]=1./(FLOAT)(nx*ny*nz);

  if ((flags&FLAG_DD)||(flags&FLAG_II))
    {
      FFTW_COMPLEX max;

      kx=fft_indgen(nx,fx,2);
      ky=fft_indgen(ny,fy,2);
      kz=fft_indgen(nz,fz,2);

      kx[0]=ky[0]=kz[0]=0.;

      if (flags&FLAG_DD)
	for (z=0,n=0;z<nzout;z++)
	  for (y=0;y<nyout;y++)  
	    for (x=0;x<nxout;x++,n++)
	      k[n]*=(kx[x]+ky[y]+kz[z]);
      else
	{
	  if (ny==1) max=0;//k[0]/(kx[1]/(nx*nx));
	  else if (nz==1) max=0;//k[0]/(kx[1]/(nx*nx)+ky[1]/(ny*ny));
	  else max=0;//1.E12+0*k[0]/(kx[1]/(nx*nx)+ky[1]/(ny*ny)+kz[1]/(nz*nz));

	  //printf ("max = %f\n");

	  for (z=0,n=0;z<nzout;z++)
	    for (y=0;y<nyout;y++)  
	      for (x=0;x<nxout;x++,n++)
		if (kx[x]+ky[y]+kz[z]!=0)
		  k[n]/=(kx[x]+ky[y]+kz[z]);
		else
		  k[n]=max;
	}

      //for (x=0;x<nout;x++) printf ("k[%d]=%lg+%lg i\n",x,creal(k[x]),cimag(k[x]));
      
      FFTWNAME(free)(kx);FFTWNAME(free)(ky);FFTWNAME(free)(kz);
    }

  order=0;
  if (flags&FLAG_DX) order+=1;if (flags&FLAG_D2X) order+=2;if (flags&FLAG_D3X) order+=3;
  if (flags&FLAG_IX) order-=1;if (flags&FLAG_I2X) order-=2;if (flags&FLAG_I3X) order-=3;
  kx = fft_indgen(nx,fx,order);

  order=0;
  if (flags&FLAG_DY) order+=1;if (flags&FLAG_D2Y) order+=2;if (flags&FLAG_D3Y) order+=3;
  if (flags&FLAG_IY) order-=1;if (flags&FLAG_I2Y) order-=2;if (flags&FLAG_I3Y) order-=3;
  ky = fft_indgen(ny,fy,order);

  order=0;
  if (flags&FLAG_DZ) order+=1;if (flags&FLAG_D2Z) order+=2;if (flags&FLAG_D3Z) order+=3;
  if (flags&FLAG_IZ) order-=1;if (flags&FLAG_I2Z) order-=2;if (flags&FLAG_I3Z) order-=3;
  kz = fft_indgen(nz,fz,order);   

  //for (x=0,n=0;x<nx;x++) printf ("%lg %lg+%lg i\n",fx,creal(kx[x]),cimag(kx[x]));
  //printf (" %lg+%lg i,  %lg+%lg i\n",creal(ky[0]),cimag(ky[0]),creal(kz[0]),cimag(kz[0]));
  
  for (z=0,n=0;z<nzout;z++)
    for (y=0;y<nyout;y++)  
      for (x=0;x<nxout;x++,n++)
      {
	  k[n]*=kx[x]*ky[y]*kz[z];
	  //printf ("%d: %lg+i*%lg  --  %lg %lg %lg\n",x,creal(k[n]),cimag(k[n]),cimag(kx[x]),cimag(ky[y]),cimag(kz[z]));
      }

  //printf ("kz[%d] -> %lg+i*%lg\n",nzout,creal(kx[nzout]),cimag(kx[nzout]));exit(0);
  //for (x=0;x<nout;x++) printf ("k[%d]=%lg+%lg i\n",x,creal(k[x]),cimag(k[x]));
  FFTWNAME(free)(kx);FFTWNAME(free)(ky);FFTWNAME(free)(kz);
  
  //Now the serious stuff
  //p = fftw_plan_dft_r2c_3d(nx,ny,nz,tab_in, out,FFTW_ESTIMATE);
  p = FFTWNAME(plan_dft_3d)(nx,ny,nz,out,out,-1,FFTW_ESTIMATE);
  //for (x=0;x<nout;x++) out[x]=0*I+tab_in[x];
  FFTWNAME(execute)(p); 
  //fftw_destroy_plan(p);
  
  for (n=0;n<nout;n++)
  {
      out[n]*=k[n];
      //printf("%lg+%lgi\n",creal(out[n]),cimag(out[n])); 
  }
  
  
  //printf ("nel=%d (101*101*51)=%d\n",sizeof(out)/(2*sizeof(double)),(101*101*51));exit(0);
  
  //for (n=0;n<nout;n++) printf ("%lg+%lg i\n",creal(out[n]),cimag(out[n]))
  
  FFTWNAME(free)(k);
  
  //p = fftw_plan_dft_c2r_3d(nx,ny,nz,out,mytab,FFTW_ESTIMATE);
  p = FFTWNAME(plan_dft_3d)(nx,ny,nz,out,out,1,FFTW_ESTIMATE);
  FFTWNAME(execute)(p); 
  //fftw_destroy_plan(p);
    
  for (x=0;x<nout;x++) mytab[x]=creal(out[x]);

  FFTWNAME(free)(out);

#ifdef FFTW_NTHREADS
  FFTWNAME(cleanup_threads)();
#endif
  //printf ("%d %d %d\n",nx,ny,nz);

  return 0;
}
