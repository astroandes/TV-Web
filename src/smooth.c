#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <hdf5.h>
#include "io.h"
#include "field.h"
#include "READ_ART.h"
#include "tools.h"

#define NORM(a) (sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]))

#ifndef PI
#define PI 3.1415926535897932384
#endif

typedef struct sort_str_typ
{
  float pos[3];
  float vel[3];
  int id;
} sort_str;

int Comp_sort_str(const void *pa,const void *pb)
{
  sort_str *a=(sort_str *)pa;
  sort_str *b=(sort_str *)pb;

  if (a->id > b->id) return 1;
  else if (a->id < b->id) return -1;
  else return 0;

}

int Comp_sort_float(const void *a, const void *b)
{
  const float *val_a = (const float *)a;
  const float *val_b = (const float *)b;

  if (*val_a > *val_b) return 1;
  else if (*val_a < *val_b) return -1;
  else return 0;

}

char *CutName(char *Name)
{
    int i,j;

    for (i=0,j=-1;i<strlen(Name);i++) 
	if (Name[i]=='/') j=i;

    if (j!=-1)
	return &Name[j+1];
    else
	return Name;
}

typedef struct EllGaussianParms_str{
  double sigma;
  double ex;
  double ey;
  double ez;
} EllGaussianParms;

double filter_EllGaussian(double u_p,double v_p,double w_p,void *parms)
{
  EllGaussianParms *p=(EllGaussianParms *)parms;
  double u,v,w,sigma;

  u=u_p*p->ex;
  v=v_p*p->ey;
  w=w_p*p->ez;
  sigma=p->sigma;

  return 
    1./pow(sqrt(2*PI)*sigma,3.) * 
    exp (-1/2.*(u*u+v*v+w*w)/(sigma*sigma));
}

void Usage(char *ExecName)
{
    printf("Usage:  %s FileName | -d DensityFileName  \n",ExecName);
    printf("        [-n <grid size>] [-s sigma] [-center] [-overd]\n");
    printf("        [-outdir <dir name>] [-degrade < 0.0<p<1.0 >]\n");
    printf("        [-eigenname <file name>]\n");
    printf("        [-contrast <c>] [-genmask] [-vel] [-test <fname>]\n");
    printf("        [-filter <filter name> <params list>]\n");
    printf("        [-artfile <art_header_file>] [-artnfiles <nfiles>][-myart <nfile>][-t threshold]\n");
    printf("  -mass: reads or fills the masses from the gadget snapshot");
    printf("  -n: the number of nodes on the grid for CIC.\n");
    printf("  -t: the threshold value for the environment.\n");
    printf("  -s: the value of sigma for smoothing in grid units.\n");
    printf("  -degrade: 1= keep everything, 0=reject everything\n");
    printf("  -contrast: density -> density^c\n");
    printf("  -vel: smooth velocity field \n");
    printf("  -genmask: Only generates a mask of the location where galaxy are present  \n");
    printf("   the resulting file is passed to skelex when analysing a non cubic survey\n");
    printf("  -filter: used to convolve with a filter instead of a gaussian.\n");
    printf("  -artfile: art_header_file, needed if the provided datafile is ART\n");
    printf("   Possible filters: EllGauss    ->   Parameters: (Ex Ey Ez Sigma)\n");
    exit(0);
}

int main(int argc,char **argv)
{
    int i,j;
    char FileName[256];
    char ART_FileName[256];
    char Opt_EigenName[256];
    char Density_FileName[256];
    int Opt_HasDensity=0;
    int Opt_HasSnap=0;
    int Opt_GridSize=-1;
    int Opt_Degrade=0;
    int Opt_Contrast=0;
    int Opt_UseFilter=0;
    
    int Opt_nx=-1;
    int Opt_ny=-1;
    int Opt_nz=-1;

    int Opt_genmask=0;

    char OutFileName[256];

    int Opt_NSmooth=0;
    float Opt_SmoothSig[100];

    int Opt_Center=0;
    int Opt_Overd=0;
    


    char Opt_OutDir[256];

    float Opt_Degval=1.0;
    float Opt_Contval=1.0;

    int Opt_Vel = 0;
    
    int Opt_ART = 0;

    int Opt_HDF5 = 0;

    int Opt_ART_NFILES = 0;

    int Opt_ART_MYFILE = -1;

    int Opt_Gadget_Mass = 0;




    int Opt_EigenVectors=0;

    char FilterName[255];
    double Total=0;
    density_grid Density;
    long long n;
    long long n_to_alloc;
    void *FilterParams;
    double (*filterfunc)(double,double,double,void *);
    EllGaussianParms EGparams;

    strcpy(Opt_OutDir,"./");
    memset(&Density,0,sizeof(Density));

    snapshot_data *Simu = calloc (1,sizeof(snapshot_data));

    Density.grid=NULL;

    for (i=1;i<argc;)
    {
	if (argv[i][0] != '-')
	{
	    strcpy(FileName,argv[i]);
	    i++;
	    Opt_HasSnap=1;
	}
	else if (!strcmp(argv[i],"-n"))
	{
	    i++;
	    Opt_GridSize = atoi(argv[i]);
	    Opt_nx=Opt_ny=Opt_nz = Opt_GridSize;
	    i++;    
	}
	else if (!strcmp(argv[i],"-d"))
	{
	    i++;
	    strcpy(Density_FileName,argv[i]);
	    Opt_HasDensity=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-s"))
	{
	    i++;
	    Opt_SmoothSig[Opt_NSmooth]=atof(argv[i]);
	    Opt_NSmooth=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-center"))
	{
	    Opt_Center=1;
	    i++;
	}     
	else if (!strcmp(argv[i],"-eigenvector"))
	  {
	    Opt_EigenVectors=1;
	    i++;
	  }       
	else if (!strcmp(argv[i],"-overd"))
	  {
	    Opt_Overd=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-genmask"))
	{
	    Opt_genmask=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-outdir"))
	{
	    i++;
	    strcpy(Opt_OutDir,argv[i]);;
	    i++;
	}
	else if (!strcmp(argv[i],"-eigenname"))
	{
	    i++;
	    strcpy(Opt_EigenName,argv[i]);;
	    i++;
	}
	else if (!strcmp(argv[i],"-degrade"))
	{
	    Opt_Degrade=1;
	    i++;
	    Opt_Degval=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-filter"))
	{
	    Opt_UseFilter=1;
	    i++;
	    strcpy(FilterName,argv[i]);
	    i++;
	    if (!strcmp(FilterName,"EllGauss"))
	      {
		EGparams.ex=atof(argv[i++]);
		EGparams.ey=atof(argv[i++]);
		EGparams.ez=atof(argv[i++]);
		EGparams.sigma=atof(argv[i++]);
		FilterParams = (void *)&EGparams;
		filterfunc=&filter_EllGaussian;
	      }
	    else
	      {
		fprintf (stderr,"\nInvalid filter name: %s\n\n",FilterName);
		Usage(argv[0]);
	      }
	    Opt_NSmooth=1;
	    
	}
	else if (!strcmp(argv[i],"-contrast"))
	{
	    Opt_Contrast=1;
	    i++;
	    Opt_Contval=atof(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-vel"))
	{
	    Opt_Vel=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-mass"))
	{
	    Opt_Gadget_Mass=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-artfile"))
	{
	    i++;
	    strcpy(ART_FileName,argv[i]);
	    Opt_ART=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-hdf5"))
	{
	    Opt_HDF5=1;
	    i++;
	}
	else if (!strcmp(argv[i],"-artnfiles"))
	{
	    i++;
	    Opt_ART_NFILES=atoi(argv[i]);
	    i++;
	}
	else if (!strcmp(argv[i],"-myart"))
	{
	    i++;
	    Opt_ART_MYFILE=atoi(argv[i]);
	    i++;
	}
	else
	{
	    fprintf (stderr,"\nUnrecognized option: %s\n\n",argv[i]);
	    Usage(argv[0]);
	}
    }
    
    if ((!Opt_HasSnap)&&(!Opt_HasDensity))
    {
	fprintf(stderr,"Please provide snapshot name or a density file name ('-d').\n");
	Usage(argv[0]);
    }


    if (Opt_OutDir[strlen(Opt_OutDir)-1]!='/')
    {
	i=strlen(Opt_OutDir);
	Opt_OutDir[i]='/';
	Opt_OutDir[i+1]='\0';
    }

    //CIC
    if ((Opt_HasSnap)&&(!Opt_Vel))
    {
	if (Opt_GridSize == -1)
	{
	    fprintf(stderr,"Please provide a grid size (option '-n') or a density file name.\n");
	    exit(0);
	}

	if (IsGadgetFile(FileName)&&!(Opt_ART)&&!(Opt_HDF5))
	{
	  if (Opt_Degrade){
		ReadGadget(FileName,Simu,FLAG_POS);
	  }else{ 
	    if(Opt_Gadget_Mass){
	      ReadGadget(FileName,Simu,FLAG_POS|FLAG_MASS);
	    }else{
	      ReadGadget(FileName,Simu,FLAG_POS);
	    }
	  }
	}
	else if(Opt_ART)
	{
	  ReadART(ART_FileName, FileName, Simu, Opt_ART_NFILES, Opt_ART_MYFILE);
	}
	else if(Opt_HDF5&& (!Opt_ART)){
	  ReadHDF5File(FileName, Simu, FLAG_POS);
	}
	else	
	{
	    if (Opt_Degrade)
		ReadSIMPLE2Gadget(FileName,Simu,FLAG_ALL);
	    else ReadSIMPLE2Gadget(FileName,Simu,FLAG_POS);
	}
	if (Simu->header.BoxSize ==0)
	    Simu->header.BoxSize =50000.;

	if (Opt_Degrade)
	{
	    FILE *file;
	    int dum;
	    char TmpFileName[256];
	    sort_str *sorted_snap;

	    n=Simu->N;
	    sorted_snap = malloc (n*sizeof(sort_str));
	    printf ("size:%d\n",(int)sizeof(sort_str));
	    printf ("Degrading data (%lld -> %d) ... ",n,(int)(Opt_Degval*n));fflush(0); 
	    srand(1);
			       
	    for (i=0;i<n;i++)
	      {
		sorted_snap[i].id = Simu->Id[i];
		memcpy(sorted_snap[i].pos,&Simu->Pos[3*i],3*sizeof(float));
		memcpy(sorted_snap[i].vel,&Simu->Vel[3*i],3*sizeof(float));
	      }
	    printf ("Sorting ...");fflush(0);
	    qsort (sorted_snap,n,sizeof(sort_str),Comp_sort_str);
	
	    printf (" Cutting ...");fflush(0);

	    for(i=0;i<(1.-Opt_Degval)*n;i++)
	    {
	      j=(Simu->N-1) *( ((double)rand())/RAND_MAX );
	      if (sorted_snap[j].id!=sorted_snap[Simu->N-1].id)
		  {
		    memcpy(&sorted_snap[j],&sorted_snap[Simu->N-1],sizeof(sort_str));
		  }
		Simu->N--;
	    }
	    for (i=0;i<Simu->N;i++)
	      {
		memcpy(&Simu->Pos[3*i],sorted_snap[i].pos,3*sizeof(float));
		memcpy(&Simu->Vel[3*i],sorted_snap[i].vel,3*sizeof(float));
		Simu->Id[i] = sorted_snap[i].id;
	      }
	    free (sorted_snap);
	    Simu->Pos = realloc (Simu->Pos,3*sizeof(float)*Simu->N);

	    printf (" Saving ...");fflush(0);
	    Simu->header.npartTotal[1] = Simu->N;
	    Simu->header.npart[1] = Simu->N;
	    Simu->header.mass[1] /= Opt_Degval;
	    
	    sprintf (TmpFileName,"./%s.degrade",CutName(FileName));
	    file=fopen(TmpFileName,"w");
	    dum=256;
	    fwrite (&dum,sizeof(int),1,file);
	    fwrite (&Simu->header,256,1,file);
	    fwrite (&dum,sizeof(int),1,file);
	    dum=Simu->N*3*sizeof(float);
	    fwrite (&dum,sizeof(int),1,file);
	    fwrite (Simu->Pos,sizeof(float),3*Simu->N,file);
	    fwrite (&dum,sizeof(int),1,file);
	    fwrite (&dum,sizeof(int),1,file);
	    fwrite (Simu->Vel,sizeof(float),3*Simu->N,file);
	    fwrite (&dum,sizeof(int),1,file);
	    fwrite (&dum,sizeof(int),1,file);
	    fwrite (Simu->Id,sizeof(int),Simu->N,file);
	    fwrite (&dum,sizeof(int),1,file);
	    fclose(file);
	    printf ("done.\n");
	    

	    return 0;
	}

	/*
	if (Opt_genmask)
	  {
	    
	    Mask->Nx=Opt_nx;
	    Mask->Ny=Opt_ny;
	    Mask->Nz=Opt_nz;
	    
	    Mask->dx=Simu->header.BoxSize/Opt_nx;
	    Mask->dy=Simu->header.BoxSize/Opt_ny;
	    Mask->dz=Simu->header.BoxSize/Opt_nz;
	    
	    Mask->x0=0;//+Mask->dx/2.;
	    Mask->y0=0;//+Mask->dy/2.;
	    Mask->z0=0;//+Mask->dz/2.;

	    if (Opt_Center)
	      {
		Mask.x0-=Simu->header.BoxSize/2.;
		Mask.y0-=Simu->header.BoxSize/2.;
		Mask.z0-=Simu->header.BoxSize/2.;
	      }
	    
	    Mask->NNodes = Nx*Ny*Nz;

	    Mask.Mask = calloc (Mask->nnodes,sizeof(FLOAT));

	    
	  }
	*/
	if(Opt_Gadget_Mass){
	  Compute_CIC_Density(Simu->Pos,Simu->Mass,1,0,Simu->N,&Density,Opt_nx,Opt_ny,Opt_nz,0,0,0,
			      Simu->header.BoxSize,Simu->header.BoxSize,Simu->header.BoxSize);
	}else{
	  Compute_CIC_Density(Simu->Pos,NULL,0,0,Simu->N,&Density,Opt_nx,Opt_ny,Opt_nz,0,0,0,
			      Simu->header.BoxSize,Simu->header.BoxSize,Simu->header.BoxSize);
	}
	
	Density.redshift = Simu->header.redshift;

	if (Opt_Center)
	{
	    Density.x0-=Simu->header.BoxSize/2.;
	    Density.y0-=Simu->header.BoxSize/2.;
	    Density.z0-=Simu->header.BoxSize/2.;
	}
	
	if (Opt_genmask)
	  {
	    sprintf(OutFileName,"%s%s.n%d.mask",Opt_OutDir,CutName(FileName),Opt_nx);
	    SaveDensityGrid(OutFileName,&Density);
	    return 0;
	  }

	if (Opt_Overd)
	{
	  double lo_mean=0;

	    for (i=0;i<Density.NNodes;i++)
	      {
		lo_mean += Density.grid[i];
	      }
	    lo_mean /= Density.NNodes;

	    for (i=0;i<Density.NNodes;i++)
	      {
		Density.grid[i]/=lo_mean;
	      }
	}

	if (Opt_Contrast)
	  {
	    printf ("Changing density contrast (d -> d^%f) ... ",Opt_Contval);
	    for (i=0;i<Density.NNodes;i++)
	      {
		Density.grid[i] = pow(Density.grid[i],Opt_Contval);
	      }
	    printf ("done.\n");
	  }

	sprintf(OutFileName,"%s%s.CIC",Opt_OutDir,CutName(FileName));
	SaveDensityGrid(OutFileName,&Density);
    }
    else if ((Opt_HasSnap)&&(Opt_Vel))
    {
	int k;
	float *velnorm;

	if (IsGadgetFile(FileName)&&!(Opt_ART))
	    ReadGadget(FileName,Simu,FLAG_POS|FLAG_VEL);
	else if(Opt_ART)
	  ReadART(ART_FileName, FileName, Simu, Opt_ART_NFILES, Opt_ART_MYFILE);
	else
	  ReadSIMPLE2Gadget(FileName,Simu,FLAG_POS|FLAG_VEL);


	/*CIC interpolation of the velocities*/
	for(i=0;i<3;i++){
	  Compute_CIC_Density(Simu->Pos,Simu->Vel,3,i,Simu->N, &Density,Opt_nx,Opt_ny,Opt_nz,0,0,0,
			      Simu->header.BoxSize,Simu->header.BoxSize,Simu->header.BoxSize);
	  sprintf(OutFileName,"%s%s.v%d.CIC",Opt_OutDir,CutName(FileName),i);
	  SaveDensityGrid(OutFileName,&Density);
	}

	Compute_CIC_Density(Simu->Pos,NULL,0,0,Simu->N,&Density,Opt_nx,Opt_ny,Opt_nz,0,0,0,
			    Simu->header.BoxSize,Simu->header.BoxSize,Simu->header.BoxSize);
	sprintf(OutFileName,"%s%s.CIC",Opt_OutDir,CutName(FileName));

	SaveDensityGrid(OutFileName,&Density);

	for (j=0;j<Opt_NSmooth;j++)
	{
	    for (i=0;i<4;i++)
	    {
		if (i==0)
		{
		    velnorm = malloc (Simu->N*sizeof(float));
		    for (k=0;k<Simu->N;k++)
			velnorm[k] = NORM((&Simu->Vel[3*k]));
		    Compute_CIC_Density(Simu->Pos,velnorm,1,0,Simu->N,&Density,Opt_nx,Opt_ny,Opt_nz,0,0,0,
					Simu->header.BoxSize,Simu->header.BoxSize,Simu->header.BoxSize);
		    free(velnorm);
		}
		else
		    Compute_CIC_Density(Simu->Pos,Simu->Vel,3,i-1,Simu->N,&Density,Opt_nx,Opt_ny,Opt_nz,0,0,0,
					Simu->header.BoxSize,Simu->header.BoxSize,Simu->header.BoxSize);

		if (Opt_Center)
		{
		    Density.x0-=Simu->header.BoxSize/2.;
		    Density.y0-=Simu->header.BoxSize/2.;
		    Density.z0-=Simu->header.BoxSize/2.;
		}
		if (!Opt_UseFilter)
		  GaussSmooth(&Density,Opt_SmoothSig[j], 0);
		else
		  FilterSmooth(&Density,filterfunc,FilterName,FilterParams);

		Density.smoothing = Opt_SmoothSig[j];

		if (!Opt_UseFilter)
		  sprintf(OutFileName,"%s%s.v%d.CIC.s%.2f",Opt_OutDir,CutName(FileName),i,Opt_SmoothSig[j]);
		else
		  sprintf(OutFileName,"%s%s.v%d.CIC.%s",Opt_OutDir,CutName(FileName),i,FilterName);

		SaveDensityGrid(OutFileName,&Density);
/*
		for (k=0;k<100;k++) 
		    printf ("%e ",Density.grid[Density.NNodes*rand()/RAND_MAX]);
		printf ("\n");
*/
	    }
	}
    }else if ((!Opt_HasSnap)&&(Opt_Vel)){




      sprintf(OutFileName,"%s%s.CIC",Opt_OutDir,CutName(Density_FileName));
      LoadDensityGrid(OutFileName,&Density);


      
      sprintf(OutFileName,"%s%s.v%d.CIC",Opt_OutDir,CutName(Density_FileName),0);
      LoadVelocityGrid(OutFileName, &Density, 0);


      sprintf(OutFileName,"%s%s.v%d.CIC",Opt_OutDir,CutName(Density_FileName),1);
      LoadVelocityGrid(OutFileName, &Density, 1);

      sprintf(OutFileName,"%s%s.v%d.CIC",Opt_OutDir,CutName(Density_FileName),2);
      LoadVelocityGrid(OutFileName, &Density, 2);





      if(Opt_NSmooth){
	puts("about to do the  gaussian smoothing");
	GaussSmooth(&Density,Opt_SmoothSig[0], DENS);
	GaussSmooth(&Density,Opt_SmoothSig[0], VEL_X);
	GaussSmooth(&Density,Opt_SmoothSig[0], VEL_Y);
	GaussSmooth(&Density,Opt_SmoothSig[0], VEL_Z);
	puts("Finished all the smoothing");
      }

      float my_tiny;
      if(Opt_NSmooth){
	my_tiny=1.0E-6;
      }else{
	my_tiny=1.0E-3;
      }

      for (n=0;n<Density.NNodes;n++)
	{
	  if(Density.grid[n]>my_tiny){
	    Density.grid_vx[n] = Density.grid_vx[n]/Density.grid[n];
	    Density.grid_vy[n] = Density.grid_vy[n]/Density.grid[n];
	    Density.grid_vz[n] = Density.grid_vz[n]/Density.grid[n];
	  }else{
	    Density.grid_vx[n] = Density.grid_vx[n]/my_tiny;
	    Density.grid_vy[n] = Density.grid_vy[n]/my_tiny;
	    Density.grid_vz[n] = Density.grid_vz[n]/my_tiny;
	  }
	}
    
      
      ComputeVelocityTorsion(&Density);      

      free(Density.grid_vx);
      free(Density.grid_vy);
      free(Density.grid_vz);

      fprintf(stdout, "The memory from the velocity is free\n");

	    /*allocate the memory*/
      if(Opt_EigenVectors){
	n_to_alloc = sizeof(float)*Density.NNodes*3;
	if(!(Density.eigenvector_1 = malloc(n_to_alloc))){
	  puts("Problem with eigenvector allocation");
	  exit(1);
	}
	if(!(Density.eigenvector_2 = malloc(n_to_alloc))){
	  puts("Problem with eigenvector allocation");
	  exit(1);
	}
	if(!(Density.eigenvector_3 = malloc(n_to_alloc))){
	  puts("Problem with eigenvector allocation");
	  exit(1);
	}
      }
      if(!Opt_EigenVectors){
	n_to_alloc = sizeof(float)*Density.NNodes;
	if(!(Density.eigenvalue_1 =  malloc(n_to_alloc)))
	  {
	    puts("problem with the eigenvalue allocation");
	    exit(1);
	  }
	
	if(!(Density.eigenvalue_2 =  malloc(n_to_alloc)))
	  {
	    puts("problem with the eigenvalue allocation");
	    exit(1);
	  }
	
	if(!(Density.eigenvalue_3 =  malloc(n_to_alloc)))
	  {
	    puts("problem with the eigenvalue allocation");
	    exit(1);
	  }
      }
	    if(!(Density.trace =  malloc(n_to_alloc)))
	    {
		puts("problem with the trace allocation");
		exit(1);
	    }

	    double Vp[9];
	    double vp[3];
	    //	    double pos[3];
	    double norm_v_web;
	    double tmp;
	    //	    double tmp_vec[3];
	    norm_v_web = -(1.0/100.0)*(1.0/(Density.dx/1000.0));
	    fprintf(stdout, "Norm v-web %f\n", norm_v_web);
	    for (n=0;n<Density.NNodes;n++)
	    {

		Diagonalise3x3(&Density.torsion[6*n],vp,Vp);
		/*the order is inverted because the sign is inverted*/
		vp[0] = vp[0] * norm_v_web;
		vp[1] = vp[1] * norm_v_web;
		vp[2] = vp[2] * norm_v_web;

		tmp = vp[0];
		vp[0] = vp[2];
		vp[2] = tmp;

		if(!Opt_EigenVectors){
		  Density.eigenvalue_1[n] = vp[0];
		  Density.eigenvalue_2[n] = vp[1];
		  Density.eigenvalue_3[n] = vp[2];
		  Density.trace[n] = vp[2] + vp[1] + vp[0];
		}



		if(Opt_EigenVectors){
		/*the order is inverted because the sign in the eigenvalues is inverted*/
		  Density.eigenvector_1[3*n+0] = Vp[6];
		  Density.eigenvector_1[3*n+1] = Vp[7];
		  Density.eigenvector_1[3*n+2] = Vp[8];
		  
		  Density.eigenvector_2[3*n+0] = Vp[3];
		  Density.eigenvector_2[3*n+1] = Vp[4];
		  Density.eigenvector_2[3*n+2] = Vp[5];
		  
		  Density.eigenvector_3[3*n+0] = Vp[0];
		  Density.eigenvector_3[3*n+1] = Vp[1];
		  Density.eigenvector_3[3*n+2] = Vp[2];
		}
		
		/*
				fprintf(stdout, "SAL %d %e %e %e %e %e %e %e\n", 
			n, Density.grid_vx[n], Density.grid_vy[n], Density.grid_vz[n], 
			vp[0], vp[1], vp[2], Density.trace[n]);	
				fprintf(stdout, "%e %e %e\n", Density.torsion[6*n], Density.torsion[6*n+1], Density.torsion[6*n+2]);
				fprintf(stdout, "%e %e %e\n", Density.torsion[6*n+1], Density.torsion[6*n+3], Density.torsion[6*n+4]);
				fprintf(stdout, "%e %e %e\n", Density.torsion[6*n+2], Density.torsion[6*n+4], Density.torsion[6*n+5]);
		*/
	    }


	    if(!Opt_EigenVectors){
	      if (Opt_HasDensity){
		if (Opt_NSmooth)
		  sprintf(OutFileName,"%s%s.s%.2f.eigen",Opt_OutDir, 
			  CutName(Density_FileName), Opt_SmoothSig[0]);
		else
		  sprintf(OutFileName,"%s%s.eigen",Opt_OutDir,
			  CutName(Density_FileName));
	      }else{
		if (Opt_NSmooth)
		  sprintf(OutFileName,"%s%s.s%.2f.eigen",Opt_OutDir,
			  CutName(FileName),Opt_SmoothSig[0]);
		else
		  sprintf(OutFileName,"%s%s.eigen",Opt_OutDir,
			  CutName(FileName));
	      }
	      
	      SaveEigenvalueGrid(OutFileName,&Density,1);	    	    	    
	      SaveEigenvalueGrid(OutFileName,&Density,2);	    	    	    
	      SaveEigenvalueGrid(OutFileName,&Density,3);

	      if (Opt_HasDensity)
		if (Opt_NSmooth)
		  sprintf(OutFileName,"%s%s.s%.2f.trace",Opt_OutDir,
			  CutName(Density_FileName), Opt_SmoothSig[0]);
		else
		  sprintf(OutFileName,"%s%s.trace",Opt_OutDir,
			  CutName(Density_FileName));
	      else
		if (Opt_NSmooth)
		  sprintf(OutFileName,"%s%s.s%.2f.trace",Opt_OutDir,
			  CutName(FileName),Opt_SmoothSig[0]);
		else
		  sprintf(OutFileName,"%s%s.trace",Opt_OutDir,
			  CutName(FileName));

	      
	      SaveTraceGrid(OutFileName,&Density);	    	      
	    }

	    if(Opt_EigenVectors){
	      if (Opt_HasDensity){
		if (Opt_NSmooth)
		  sprintf(OutFileName,"%s%s.s%.2f.eigenvec",Opt_OutDir, 
			  CutName(Density_FileName), Opt_SmoothSig[0]);
		else
		  sprintf(OutFileName,"%s%s.eigenvec",Opt_OutDir,
			  CutName(Density_FileName));
	      }else{
		if (Opt_NSmooth)
		  sprintf(OutFileName,"%s%s.s%.2f.eigenvec",Opt_OutDir,
			  CutName(FileName),Opt_SmoothSig[0]);
		else
		  sprintf(OutFileName,"%s%s.eigenvec",Opt_OutDir,
			  CutName(FileName));
	      }
	      
	      SaveEigenvectorGrid(OutFileName,&Density,1);	    	    	    
	      SaveEigenvectorGrid(OutFileName,&Density,2);	    	    	    
	      SaveEigenvectorGrid(OutFileName,&Density,3);	      
	    }

    }
    else{
      sprintf(OutFileName,"%s%s.CIC",Opt_OutDir,CutName(Density_FileName));
      LoadDensityGrid(OutFileName,&Density);
    }

    for (i=0,Total=0;i<Density.NNodes;i++){
	Total += Density.grid[i];
    }
    fprintf(stdout, "00 Number of nodes %lld\n", Density.NNodes);


    if (!Opt_Vel){
	//smooth

      if(Opt_NSmooth){
	puts("about to do the smoothing0");
	GaussSmooth(&Density,Opt_SmoothSig[0], DENS);
      }

      puts("getting the overdensity1");
      double lo_mean=0;
      
      for (j=0;j<Density.NNodes;j++){
	lo_mean += Density.grid[j];	
      }	
      lo_mean /= Density.NNodes;

      for (j=0;j<Density.NNodes;j++){	
	Density.grid[j]/=lo_mean;
	Density.grid[j]-=1.0;
      }
      fprintf(stdout, "got the overdensity1\n");
      fflush(stdout);

      if (Opt_HasDensity)
	if (Opt_NSmooth)
	  sprintf(OutFileName,"%s%s.s%.2f.DELTA",Opt_OutDir,CutName(Density_FileName),Opt_SmoothSig[0]);
	else
	  sprintf(OutFileName,"%s%s.DELTA",Opt_OutDir,CutName(Density_FileName));
      else
	if (Opt_NSmooth)
	  sprintf(OutFileName,"%s%s.s%.2f.DELTA",Opt_OutDir,CutName(FileName),Opt_SmoothSig[0]);
	else
	  sprintf(OutFileName,"%s%s.DELTA",Opt_OutDir,CutName(FileName));

      /*
      for (n=0;n<Density.NNodes;n++){
		Density.grid[n] = Density.grid[n]/(1.0*Density.NNodes);
	}
      */
      if (Opt_NSmooth){
	Density.smoothing = Opt_SmoothSig[0];
	SaveDensityGrid(OutFileName,&Density);
      }
      

      fprintf(stdout, "computing potential 2\n");
      fflush(stdout);
      ComputePotential(&Density);	   
      fprintf(stdout, "computed the potential 2\n");            
      fflush(stdout);
      
      for (n=0;n<Density.NNodes;n++){
	//works for 256
	Density.grid[n] = Density.grid[n] * pow((1.0*Density.NNodes), 0.75);
      }
      
      fprintf(stdout, "computing hessian density\n");

      ComputeDensityHessian(&Density, USE_FD);
      
      
	    
      if(Opt_EigenVectors){
	n_to_alloc = sizeof(float)*Density.NNodes*3;
	if(!(Density.eigenvector_1 = malloc(n_to_alloc))){
	  puts("Problem with eigenvector allocation");
	  exit(1);
	}
	if(!(Density.eigenvector_2 = malloc(n_to_alloc))){
	  puts("Problem with eigenvector allocation");
	  exit(1);
	}
	if(!(Density.eigenvector_3 = malloc(n_to_alloc))){
	  puts("Problem with eigenvector allocation");
	  exit(1);
	}
      }
      
      if(!Opt_EigenVectors){
	n_to_alloc = sizeof(float)*Density.NNodes;
	/*allocate the memory*/
	if(!(Density.eigenvalue_1 =  malloc(n_to_alloc))){	    
	  puts("problem with the eigenvalue allocation");
	  exit(1);
	}
	
	if(!(Density.eigenvalue_2 =  malloc(n_to_alloc))){
	  puts("problem with the eigenvalue allocation");
	  exit(1);
	}
	
	if(!(Density.eigenvalue_3 =  malloc(n_to_alloc))){
	  puts("problem with the eigenvalue allocation");
	  exit(1);
	}
	
	if(!(Density.trace =  malloc(n_to_alloc))){
	  puts("problem with the trace allocation");
	  exit(1);
	}
      }
      
      /*now make the calculation of the igenvalues*/
      double Vp[9];
      double vp[3];
      double norm_t_web;
      float *ratio;

      n_to_alloc = sizeof(float)*Density.NNodes;
      if(!(ratio=malloc(n_to_alloc))){
	fprintf(stdout, "Problem with ratio allocation\n");
	exit(1);
      }
      
      /*gets the normalization as the mean of 'density to trace' ratio*/
      /*
      for (n=0;n<Density.NNodes;n++)
	{
	  Diagonalise3x3(&Density.hessian[6*n],vp,Vp);
	  ratio[n] = fabs(Density.grid[n]/(vp[2]+vp[1]+vp[0]));
	}      
      qsort (ratio,Density.NNodes,sizeof(float), Comp_sort_float);

      norm_t_web = ratio[(int)(Density.NNodes/2)];
      
      fprintf(stdout, "%f %f %f %f \n", ratio[0], ratio[1], ratio[2], ratio[Density.NNodes-1]);
      fprintf(stdout, "Norm T Web = %f\n", norm_t_web);
      free(ratio);
      */
      norm_t_web = 1.0;
      for (n=0;n<Density.NNodes;n++)
	{
	  Diagonalise3x3(&Density.hessian[6*n],vp,Vp);
	  
	  vp[2] = vp[2] * norm_t_web;
	  vp[1] = vp[1] * norm_t_web;
	  vp[0] = vp[0] * norm_t_web;
	  
	  if(Opt_EigenVectors){
	    Density.eigenvector_1[3*n+0] = Vp[0];
	    Density.eigenvector_1[3*n+1] = Vp[1];
	    Density.eigenvector_1[3*n+2] = Vp[2];
	    
	    Density.eigenvector_2[3*n+0] = Vp[3];
	    Density.eigenvector_2[3*n+1] = Vp[4];
	    Density.eigenvector_2[3*n+2] = Vp[5];
	    
	    Density.eigenvector_3[3*n+0] = Vp[6];
	    Density.eigenvector_3[3*n+1] = Vp[7];
	    Density.eigenvector_3[3*n+2] = Vp[8];
	    
	  }



		/*here I try to calculate the positions*/
	  if(!Opt_EigenVectors){
	    Density.eigenvalue_1[n] = vp[0];
	    Density.eigenvalue_2[n] = vp[1];
	    Density.eigenvalue_3[n] = vp[2];
	    Density.trace[n] = vp[2] + vp[1] + vp[0];
	    
	  }
	}
      
      if(!Opt_EigenVectors){
	if (Opt_HasDensity)
	  if (Opt_NSmooth)
	    sprintf(OutFileName,"%s%s.s%.2f.eigen",Opt_OutDir, 
		    CutName(Density_FileName),Opt_SmoothSig[0]);
	  else
	    sprintf(OutFileName,"%s%s.eigen",Opt_OutDir,
		    CutName(Density_FileName));
	else
	  if (Opt_NSmooth)
	    sprintf(OutFileName,"%s_%s.s%.2f.eigen",Opt_OutDir,
		    CutName(FileName),Opt_SmoothSig[0]);
	  else
	    sprintf(OutFileName,"%s%s.eigen",Opt_OutDir,
		    CutName(FileName));
	
	
	
	SaveEigenvalueGrid(OutFileName,&Density,1);	    	    	    
	SaveEigenvalueGrid(OutFileName,&Density,2);	    	    	    
	SaveEigenvalueGrid(OutFileName,&Density,3);	    	    	    
	
	if (Opt_HasDensity)
	  if (Opt_NSmooth)
	    sprintf(OutFileName,"%s%s.s%.2f.trace",Opt_OutDir, 
		    CutName(Density_FileName),Opt_SmoothSig[0]);
	  else
	    sprintf(OutFileName,"%s%s.trace",Opt_OutDir,
		    CutName(Density_FileName));
	else
	  if (Opt_NSmooth)
	    sprintf(OutFileName,"%s_%s.s%.2f.trace",Opt_OutDir,
		    CutName(FileName),Opt_SmoothSig[i]);
	  else
	    sprintf(OutFileName,"%s%s.trace",Opt_OutDir,
		    CutName(FileName));
	
	
	SaveTraceGrid(OutFileName,&Density);	    	    	    
      }
      
      
      if(Opt_EigenVectors){
	if (Opt_HasDensity){
	  if (Opt_NSmooth)
	    sprintf(OutFileName,"%s%s.s%.2f.eigenvec",Opt_OutDir, 
		    CutName(Density_FileName), Opt_SmoothSig[0]);
	  else
	    sprintf(OutFileName,"%s%s.eigenvec",Opt_OutDir,
		    CutName(Density_FileName));
	}else{
	  if (Opt_NSmooth)
	    sprintf(OutFileName,"%s%s.s%.2f.eigenvec",Opt_OutDir,
		    CutName(FileName),Opt_SmoothSig[0]);
	  else
	    sprintf(OutFileName,"%s%s.eigenvec",Opt_OutDir,
		    CutName(FileName));
	}
	
	SaveEigenvectorGrid(OutFileName,&Density,1);	    	    	    
	SaveEigenvectorGrid(OutFileName,&Density,2);	    	    	    
	SaveEigenvectorGrid(OutFileName,&Density,3);	      
      }
      
      
    }
    
    
    return 0;
}
