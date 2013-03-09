#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sph.h"

//Monaghan and lattanzio 1985
double SPH_Kernel(double d,double h)
{
  double val=(double)1./(PI*h*h*h);
  double v=d/h;

  //printf("v=%f (%f %f)\n",v,d,h);

  if (v<1)
    val *= 1-1.5*(v*v)+0.75*(v*v*v);
  else if (v<2)
    val *= 0.25*pow(2-v,3.);
  else return 0;
  //printf ("valk=%e\n",val);
  return val;
}

//Used to compute d/dx_i()
double SPH_DiffKernel(double dv[3],double d,int i,double h)
{
  double val=(double)1./(PI*h*h*h);
  double v;

  v = d/h;

  if (v<1)
    {
      val *= (dv[i]*(9./4.*v-3.))/(h*h);
    }
  else if (v<2)
    val *= (-pow(2-v,2)*dv[i]/(v*h*h) * 3./4.);
  else 
      return 0;

  return val;
}

double SPH_ComputeDensity(Nei_List *list,int N, double (*kernel)(double,double),double h,int include)
{
  double d=0;
  int i;

  if (include) d = kernel(0,h);
  for (i=0;i<N;i++)
    d+=kernel(list[i].d,h);
    
  //printf ("sd=%f %f -> %lg\n",list[0].d,h,SPH_Kernel(list[0].d,h));
  return d;
}

//Estimate h
//d returns the density
double SPH_Estimateh(Nei_List *list,int N,double(*kernel)(double,double),double *sph_d,int include)
{
  double h;
  int i=0;
  double hi;
  double oldh;
  double d;
  //double deltah=1.E-3;
  double deltah=0.1;

  
  
  h=list[N-1].d/2;
  hi=h;
  deltah*=h;
  /*
  for (i=0;i<N;i++)
    printf ("d(%d)=%f\n",i,list[i].d);
  return 0;
  */
  
  do
    {
      i++;
      d=SPH_ComputeDensity(list,N,kernel,h,include);
      //printf ("h=%f\n",h);
      oldh=h;
      h = 1./2.*pow(3./(4*PI) * N/d,1./3.);
    } while (fabs(oldh-h)>deltah);

  

  if (sph_d!=NULL) 
    *sph_d = (float)d;
    
  //printf ("%d (%f -> %f)\n",i,hi,h);

  return h;
}

