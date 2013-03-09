#ifndef __Struct_IO__
#define __Struct_IO__

typedef struct snapshot_header_str
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];  /* fills to 256 Bytes */
} snapshot_header;


typedef struct snapshot_data_str 
{
    snapshot_header header;
    float  *Pos;
    float  *Vel;
    float  *Mass;
    int    *Id;
    char   *Type;
    
    float  *Rho, *U, *Temp, *Ne;
    long long N;
    int flags;
} snapshot_data;

typedef struct cluster_data_str
{
    char FileName[256];
    int N;
    int NGroups;
    float Eps;
    //int *InGroup;
    int **InGroup;
    int *NInGroup;

} cluster_data;

typedef struct SurveyMask_struct
{
    float almin,almax;
    float abmin,abmax;
    int Nal,Nab;
    float dal,dab;

    unsigned char *Mask;//what pixel is observed or not (true/false))
} SurveyMask;


typedef struct SurveyCone_struct
{
    int N;//Number of galaxies in the survey
    SurveyMask Mask;//The mask of the survey

    int al,ab;//index corresponding to where is al,ab, ... in Data
    int v,z,d;
    int absmag[10],appmag[10];

    float h;//hubble constant (0.6<h<0.85)

    char filterName[10][30];//Name of the filter for mag
    int NFilters;
    char fromfile[256];//the file it was loaded from
    
    int NRows;//Number of raws in data
    
    float **Data;//The values
    float *Pos;
    float *Vel;
  //char *DataName[64];//A string to describe data

} SurveyCone;

#include "tree_struct.h"
#include "density_struct.h"

#endif
