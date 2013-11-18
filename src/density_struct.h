#ifndef __DENSITY_STRUCT__
#define __DENSITY_STRUCT__

#include "struct.h"

#ifdef SIMPLEPRECISION
#define FLOAT float
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN fftwf_plan
#define FFTWNAME(x) fftwf_##x
#else
#define FLOAT double
#define FFTW_COMPLEX fftw_complex
#define FFTW_PLAN fftw_plan
#define FFTWNAME(x) fftw_##x
#endif

#define VEL_X 1
#define VEL_Y 2
#define VEL_Z 3
#define DENS 0


typedef struct density_grid_str
{
  int Nx,Ny,Nz;//Number of nodes along x,y and z.
  long long NNodes;//number of nodes
  
  float x0,y0,z0;//start coordinates
  float dx,dy,dz;//grid spacing
  
  FLOAT *grid;//value at every node ...
   
  int HasGradiant;//is gradiant computed ?
  FLOAT *grad;//gradient at every node (gx1,gy1,gz1,gx2,gy2,gz2,...)
  int HasHessian;
  FLOAT *hessian;//Hessian (xx1,xy1,xz1,yy1,yz1,zz1,xx2 ...)
  int maxindex;//index of the maximum value of grid (computed with the gradiant)
  float *eigenvalue_1;//eigenvalues of the hessian 
  float *eigenvalue_2;//eigenvalues of the hessian 
  float *eigenvalue_3;//eigenvalues of the hessian 
  FLOAT *trace;//sum of the eigenvalues
  FLOAT *grid_vx; // velocity x-component
  FLOAT *grid_vy; // velocity y-component
  FLOAT *grid_vz; // velocity z-component
  FLOAT *torsion;
  float *eigenvector_1;
  float *eigenvector_2;
  float *eigenvector_3;
  int HasTorsion;
  int *environment; //evironment for each point

  float redshift;
  float smoothing;
} density_grid;

typedef struct density_field_str
{
    snapshot_header header;//gadget header

    long long NNodes;//number of nodes on the grid
    float delta;//distance between two nodes (should be <= sigma)
    int Nx,Ny,Nz;//nodes along X,Y and Z
    
    float sigma;//grid spacing and gaussian sigma
    float nsigma;//gaussian is computed up to r=nsigma*sigma
    double Norm,B;//Gauss=Norm*exp(B*r²)
    float Mass;//Mass of a particle (0 or 1 for galaxies)

    float *grid;//density value at a node
    float *pos;//positions of the nodes(only when delta==0)
    int *NSpec;//Number of spectra at a given node if (strlen(header.fill)==2)
    int **Spec;//Spectra indexes if (strlen(header.fill)==2)
    float **SpecWeight;//The probability this spec has to be attributed (btwn 0..1)
} density_field;

#endif
