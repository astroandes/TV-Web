#ifndef UTILS_H
#define UTILS_H

#include "struct.h"
#include "mcubes.h"
#include "density_struct.h"

#define SIGN(x) ((x>0)?1:((x==0)?0:-1))
#define  MAXMIN(x1, x2, min) (x1 >= x2 ? (min = x2, x1) : (min = x1, x2))

#ifndef PI
#define PI 3.1415926535897932384
#endif

int Diagonalise3x3(double *M,double *L,double *vec);
int Line2plane_Inter(float *p1,float *p2, float *p3, float *n,float *inter);
int Find_2Faces_Inter(float *p_p1, float *p_q1, float *p_r1,float *p_p2, float *p_q2, float *p_r2,float *seg);
int SegIntersect(float *seg1,float *seg2,float *ipt);
double GetPointDensity(density_grid *density,float x,float y,float z);
int InterpolateVector(FLOAT *val,double *result,int size,int dx,int dy,int dz,double ddu,double ddv,double ddw);
int ExtractDensityBlock(density_grid *density,density_grid **p_sub_density,int imin,int imax,int jmin,int jmax,int kmin,int kmax);
int Find_3Surf_Inter_Points(MC_Object *s1,MC_Object *s2,MC_Object *s3,float **itrlist);

#endif
