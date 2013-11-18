/***********************************************************

  Date        : Wed Feb 20 11:26:51 CET 2008
  Author      : Jaime Forero
  Description : Reads an ART snapshot stored in one file and 
                puts it into a gadget structure.

  Notes       : It is based on the  code READ_ART.f
  Updates     : Sun Oct 23 01:01:52 CEST 2011
                Updates to read ART in splitted files


***********************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "endian.h"
#include "struct.h"
#include "READ_ART.h"

void replace_string(char *original, char *search, char *replace){
  static char buffer[4096];
  char *p;

  if(!(p = strstr(original, search)))  // Is 'orig' even in 'str'?
    return;

  strncpy(buffer, original, p-original); // Copy characters from 'str' start to 'orig' st$
  buffer[p-original] = '\0';

  sprintf(buffer+(p-original), "%s%s", replace, p+strlen(search));
  strcpy(original, buffer);

  return;
}


int read_art_snap(char * namein, char * datain, snapshot_data *P)
{
#define SKIP  {fread(&dummy,sizeof(int),1,filein);}
  FILE * filein;
  FILE * datafile;
  long long Nparticles;
  long long i, jstart, jend, j, k;
  int dummy;
  float aexpn, aexp0, amplt, astep;
  int istep;
  float partw, tintg, ekin, ekin1, ekin2, au0, aeu0;
  int nrowc, ngridc, nspecies, nseed;
  float om0, oml0, hubble, wp5, ocurv;
  float extras[100];
  float ww;
  char header[45];
  long long * lspecies;
  float * wspecies;
  float box;
  float *buffer;
  
  
  /*replaces EQUIVALENCE statement in f77 code*/
  lspecies = (long long *)(&extras[10]);
  wspecies = &extras[0];


  fprintf(stdout, "man sizeof (int) %d\n", sizeof(int));
  fprintf(stdout, "man sizeof (float) %d\n", sizeof(float));

  if(!(filein = fopen(namein, "r")))
    {
      printf("problem opening file  %s \n", namein);
    }

  if(!(datafile = fopen(datain, "r")))
    {
      printf("problem opening file  %s \n", datain);
    }

  SKIP;
  fread(header, sizeof(char), 45, filein);
  fread(&aexpn, sizeof(float), 1, filein); 
  fread(&aexp0, sizeof(float), 1, filein); 
  fread(&amplt, sizeof(float), 1, filein);
  fread(&astep, sizeof(float), 1, filein);
  fread(&istep, sizeof(int), 1, filein);
  fread(&partw, sizeof(float), 1, filein);
  fread(&tintg, sizeof(float), 1, filein);
  fread(&ekin, sizeof(float), 1, filein);
  fread(&ekin1, sizeof(float), 1, filein); 
  fread(&ekin2, sizeof(float), 1, filein);
  fread(&au0, sizeof(float), 1, filein); 
  fread(&aeu0, sizeof(float), 1, filein); 
  fread(&nrowc, sizeof(int), 1, filein); 
  fread(&ngridc, sizeof(int), 1, filein);
  fread(&nspecies, sizeof(int), 1, filein);
  fread(&nseed, sizeof(int), 1, filein);
  fread(&om0, sizeof(float), 1, filein);
  fread(&oml0, sizeof(float), 1, filein); 
  fread(&hubble, sizeof(float), 1, filein); 
  fread(&wp5, sizeof(float), 1, filein);
  fread(&ocurv, sizeof(float), 1, filein);
  fread(extras, sizeof(float), 100, filein);
  SKIP;
#ifdef SWAP	     
  aexpn = swapF(aexpn);
  aexp0 = swapF(aexp0);
  amplt = swapF(amplt);
  astep = swapF(astep);
  istep = swapI(istep);
  partw = swapF(partw);
  tintg = swapF(tintg);
  ekin = swapF(ekin);
  ekin1 = swapF(ekin1);
  ekin2 = swapF(ekin2);
  au0 = swapF(au0);
  aeu0 = swapF(aeu0);
  nrowc = swapI(nrowc);
  ngridc = swapI(ngridc);
  nspecies = swapI(nspecies);
  nseed = swapI(nseed);
  om0 = swapF(om0);
  oml0 = swapF(oml0);
  hubble = swapF(hubble);
  wp5 = swapF(wp5);
  ocurv = swapF(ocurv);
  for(i=0;i<100;i++)
    {
      extras[i] = swapF(extras[i]);
    }
#endif



 if(NROW!=nrowc) 
   {
     printf("NROW %d in paparameter and nrowc %d in file are different\n", NROW, nrowc);
         exit(1);
   }
 
 if(NGRID!=ngridc) 
   {
     printf("NGRID %d in paparameter and ngridc %d in file are different\n", NGRID, ngridc);
         exit(1);
   }
  
  printf("header %s\n", header);
  printf("aexpn %f aexp0 %f amplt %f astep %f istep %d\n", aexpn, aexp0, amplt, astep, istep);
  printf("partw %f tintg %f ekin %f ekin1 %f ekin2 %f\n", partw, tintg, ekin, ekin1, ekin2);
  printf("au0 %f aeu0 %f nrowc %d\n", au0,aeu0,nrowc);
  printf("ngridc %d nspecies %d nseed %d\n", ngridc, nspecies, nseed);
  printf("om0 %f oml0 %f hubble %f wp5 %f ocurv %f\n", om0, oml0, hubble, wp5, ocurv);

  for(i=0;i<100;i++){
    //    printf("extras %f wspecies %f\n", extras[i], wspecies[0]);
  }
  if(extras[99]>0.)
    box =extras[99];
  else
    {
      puts("needs box size in comoving Mpc/h =");
      exit(1);
    }
  printf("box [Mpc/h]%f\n", box);

  /* SetWeights for all particles (multi masses possible)*/
  if(nspecies==0)/*old constant weigths*/
    {
      Nparticles = (NROW*NROW*NROW);    
      if(!(P->Pos = malloc(3*Nparticles*sizeof(float))))
	{
	  printf("problem allocating memory for particles positions (%d particles)\n", Nparticles);
	  exit(1);
	}
      if(!(P->Vel = malloc(3*Nparticles*sizeof(float))))
	{
	  printf("problem allocating memory for particles velocities\n");
	  exit(1);
	}
      if(!(P->Mass = malloc(3*Nparticles*sizeof(float))))
	{
	  printf("problem allocating memory for particles masses\n");
	  exit(1);
	}
      ww = 1.0;
    }
  else
    {
      Nparticles = lspecies[nspecies-1];
      printf("Nparticles %lld\n", Nparticles);
      jstart = 1;
      if(Nparticles<0)
      {
	  printf("wrong number of particles");
	  printf("Npart %d Nspecies %d N %d ", Nparticles, nspecies,(int)(extras[9 + nspecies]));
	  exit(1);
      }

      if(!(P->Pos = malloc(3*Nparticles*sizeof(float))))
      {
	  printf("problem allocating memory for particles positions (%d particles)\n", Nparticles);
	  exit(1);
      }
      if(!(P->Vel = malloc(3*Nparticles*sizeof(float))))
      {
	  printf("problem allocating memory for particles velocities\n");
	  exit(1);
      }
      if(!(P->Mass = malloc(3*Nparticles*sizeof(float))))
      {
	  printf("problem allocating memory for particles masses\n");
	  exit(1);
      }
      
      for(j=1;j<=nspecies;j++)
      {
	  
	jend = lspecies[j-1];
	printf("%d %d\n", jstart, jend);
	for(k=jstart;k<=jend;k++)
	  {
	    P->Mass[k-1] = wspecies[j-1];
	  }
	jstart = jend;
      }            
    }
  printf("Nparticles %lld\n", Nparticles);

  


  float xscale;// Scale for comoving coordinates
  float vscale;// Scale for velocities
  long long npages, n_in_last;
  
   xscale= box/NGRID ;
   vscale= 100.0*xscale/aexpn;
   
   if(nspecies==0)
     {
       npages =NROW;
       n_in_last = NPAGE;       
     }
   else
     {
       Nparticles = lspecies[nspecies-1];
       npages = (Nparticles - 1)/NPAGE  + 1 ;
       n_in_last = Nparticles - NPAGE*(npages - 1);
     }


   printf("Pages %lld Species %d n_in_last %d \n", npages, nspecies, n_in_last);

   /*allocate the buffer to read the data from the file*/
   if(!(buffer = malloc(NRECL*sizeof(float))))
     {
       printf("problem allocation memory for read buffer\n");
       exit(1);
     }
   int ipage;
   long long ntot;
   long long i_in_page;


   ntot = 0;
   for(ipage=1;ipage<=npages;ipage++)
     {
       if(ipage<npages)
	 {
	   i_in_page = NPAGE;
	 }
       else
	 {
	   i_in_page = n_in_last;
	 }	
       
       SKIP;
       fread(buffer, sizeof(float), i_in_page*6, datafile);       
       SKIP;
       printf("read page %d\n",ipage);fflush(0);

       for(i=0;i< i_in_page;i++){
//       if(!(i%100000))printf("%d %d %lld\n", i, i_in_page, (i + 3*NPAGE));fflush(0);
#ifdef SWAP 
	   P->Pos[3*(ntot + i) + 0] = (swapF(buffer[i    ])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 1] = (swapF(buffer[i + NPAGE])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 2] = (swapF(buffer[i + 2*NPAGE])-1.0)*xscale*1000.0;
	   P->Vel[3*(ntot + i) + 0] = swapF(buffer[i + 3*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 1] = swapF(buffer[i + 4*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 2] = swapF(buffer[i + 5*NPAGE])*vscale;
#else
	   P->Pos[3*(ntot + i) + 0] = ((buffer[i    ])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 1] = ((buffer[i + NPAGE])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 2] = ((buffer[i + 2*NPAGE])-1.0)*xscale*1000.0;
	   P->Vel[3*(ntot + i) + 0] = (buffer[i + 3*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 1] = (buffer[i + 4*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 2] = (buffer[i + 5*NPAGE])*vscale;
#endif
       }

       ntot = ntot + i_in_page;                 
       printf("read page %d, ntot %lld\n",ipage, ntot);fflush(0);    
     }



   printf("finished reading file %d particles in total \n", ntot);


         



  memset(P->header.npart,0,6*sizeof(int));
  memset(P->header.npartTotal,0,6*sizeof(int));
  memset(P->header.mass,0,6*sizeof(double));

  P->header.time = aexpn/1.0;
  P->header.redshift = aexpn;
  P->N = P->header.npart[1] = P->header.npartTotal[1] = ntot;
  P->header.mass[1]=27.8782*pow(box,3)/Nparticles
      *pow(hubble/100.,3.)*om0;             //unit: 10^(10) Msol/h
  P->header.BoxSize = box*1000.0;   //box size in kpc/h
  P->header.Omega0 = om0;
  P->header.OmegaLambda = oml0;
  P->header.HubbleParam = hubble;



   printf("starts the test\n");
   /* test */
   float minx=1.E10, miny=1.E10, minz=1.E10, minvx=1.E10, minvy=1.E10, minvz=1.E10;
   float maxx=-1.E10, maxy=-1.E10, maxz=-1.E10, maxvx=-1.E10, maxvy=-1.E10, maxvz=-1.E10;
   for(i=0;i<ntot;i++)
     {
       minx = MIN(minx,P->Pos[3*i + 0]);
       miny = MIN(miny,P->Pos[3*i + 1]);
       minz = MIN(minz,P->Pos[3*i + 2]);
       minvx = MIN(minvx,P->Vel[3*i + 0]);
       minvy = MIN(minvy,P->Vel[3*i + 1]);
       minvz = MIN(minvz,P->Vel[3*i + 2]);

       maxx = MAX(maxx,P->Pos[3*i + 0]);
       maxy = MAX(maxy,P->Pos[3*i + 1]);
       maxz = MAX(maxz,P->Pos[3*i + 2]);
       maxvx = MAX(maxvx,P->Vel[3*i + 0]);
       maxvy = MAX(maxvy,P->Vel[3*i + 1]);
       maxvz = MAX(maxvz,P->Vel[3*i + 2]);
     }   
   
   printf("xmin, xmax = %g %g h-1 Mpc\n", minx, maxx);
   printf("ymin, ymax = %g %g h-1 Mpc\n", miny, maxy);
   printf("zmin, zmax = %g %g h-1 Mpc\n", minz, maxz);

   printf("vxmin, vxmax = %g %g km/s \n", minvx, maxvx);
   printf("vymin, vymax = %g %g km/s\n", minvy, maxvy);
   printf("vzmin, vzmax = %g %g km/s\n", minvz, maxvz);


   fclose(datafile);
   fclose(filein);
  return ntot; 
}

int read_art_snap_multifile(char * namein, char * datain, snapshot_data *P, int n_files, int my_snap)
{
#define SKIP  {fread(&dummy,sizeof(int),1,filein);}
#define FLOAT float
  FILE * filein;
  FILE * datafile;
  long long Nparticles;
  long long i, jstart, jend, j, k;
  int dummy;
  float aexpn, aexp0, amplt, astep;
  int istep;
  float partw, tintg, ekin, ekin1, ekin2, au0, aeu0;
  int nrowc, ngridc, nspecies, nseed;
  float om0, oml0, hubble, wp5, ocurv;
  float extras[100];
  float ww;
  char header[45];
  long long * lspecies;
  float * wspecies;
  float box;
  FLOAT *buffer;
  int i_file;
  int n_p_f;
  int ipage;
  long long ntot;
  long long i_in_page;
  char newnamein[512];
  char number[12];
  float xscale;// Scale for comoving coordinates
  float vscale;// Scale for velocities
  long long npages, n_in_last;

  /*replaces EQUIVALENCE statement in f77 code*/
  lspecies = (long long *)(&extras[10]);
  wspecies = &extras[0];
  
  fprintf(stdout, "Nfiles %d\n", n_files);

  fprintf(stdout, "man sizeof (int) %d\n", sizeof(int));
  fprintf(stdout, "man sizeof (float) %d\n", sizeof(float));

  /*loope over the files*/
  for(i_file=0;i_file<n_files;i_file++){
    sprintf(number, "%d", i_file);
    strcpy(newnamein, datain);
    replace_string(newnamein, "NUM", number);
    fprintf(stdout, "%s \n", newnamein);
  }

   /*allocate the buffer to read the data from the file*/
   if(!(buffer = malloc(NRECL*sizeof(FLOAT))))
     {
       printf("problem allocation memory for read buffer\n");
       exit(1);
     }

   
   

   for(i_file=0;i_file<n_files;i_file++){
     fprintf(stdout, "i_file %d\n", i_file);
     if(my_snap==i_file && my_snap>-1){
       fprintf(stdout, "Inside i_file %d\n", i_file);
       sprintf(number, "%d", i_file);
       strcpy(newnamein, datain);
       replace_string(newnamein, "NUM", number);
       fprintf(stdout, "%s %s %s\n", namein, newnamein, number);
       
       if(!(filein = fopen(namein, "r")))
	 {
	   printf("problem opening file  %s \n", namein);
	   exit(1);
	 }
    
       if(!(datafile = fopen(newnamein, "r")))
	 {
	   printf("problem opening file  %s \n", datain);
	   exit(1);
	 }
       
       SKIP;
       fread(header, sizeof(char), 45, filein);
       fread(&aexpn, sizeof(float), 1, filein); 
       fread(&aexp0, sizeof(float), 1, filein); 
       fread(&amplt, sizeof(float), 1, filein);
       fread(&astep, sizeof(float), 1, filein);
       fread(&istep, sizeof(int), 1, filein);
       fread(&partw, sizeof(float), 1, filein);
       fread(&tintg, sizeof(float), 1, filein);
       fread(&ekin, sizeof(float), 1, filein);
       fread(&ekin1, sizeof(float), 1, filein); 
       fread(&ekin2, sizeof(float), 1, filein);
       fread(&au0, sizeof(float), 1, filein); 
       fread(&aeu0, sizeof(float), 1, filein); 
       fread(&nrowc, sizeof(int), 1, filein); 
       fread(&ngridc, sizeof(int), 1, filein);
       fread(&nspecies, sizeof(int), 1, filein);
       fread(&nseed, sizeof(int), 1, filein);
       fread(&om0, sizeof(float), 1, filein);
       fread(&oml0, sizeof(float), 1, filein); 
       fread(&hubble, sizeof(float), 1, filein); 
       fread(&wp5, sizeof(float), 1, filein);
       fread(&ocurv, sizeof(float), 1, filein);
       fread(extras, sizeof(float), 100, filein);
       SKIP;
#ifdef SWAP	     
    aexpn = swapF(aexpn);
    aexp0 = swapF(aexp0);
    amplt = swapF(amplt);
    astep = swapF(astep);
    istep = swapI(istep);
    partw = swapF(partw);
    tintg = swapF(tintg);
    ekin = swapF(ekin);
    ekin1 = swapF(ekin1);
    ekin2 = swapF(ekin2);
    au0 = swapF(au0);
    aeu0 = swapF(aeu0);
    nrowc = swapI(nrowc);
    ngridc = swapI(ngridc);
    nspecies = swapI(nspecies);
    nseed = swapI(nseed);
    om0 = swapF(om0);
    oml0 = swapF(oml0);
    hubble = swapF(hubble);
    wp5 = swapF(wp5);
    ocurv = swapF(ocurv);
    for(i=0;i<100;i++)
      {
	extras[i] = swapF(extras[i]);
      }
#endif
    
    
    
    if(NROW!=nrowc) 
      {
	printf("NROW %d in paparameter and nrowc %d in file are different\n", NROW, nrowc);
     //         exit(1);
      }
    
    if(NGRID!=ngridc) 
      {
	printf("NGRID %d in paparameter and ngridc %d in file are different\n", NGRID, ngridc);
	exit(1);
      }
    
    printf("header %s\n", header);
    printf("aexpn %f aexp0 %f amplt %f astep %f istep %d\n", aexpn, aexp0, amplt, astep, istep);
    printf("partw %f tintg %f ekin %f ekin1 %f ekin2 %f\n", partw, tintg, ekin, ekin1, ekin2);
    printf("au0 %f aeu0 %f nrowc %d\n", au0,aeu0,nrowc);
    printf("ngridc %d nspecies %d nseed %d\n", ngridc, nspecies, nseed);
    printf("om0 %f oml0 %f hubble %f wp5 %f ocurv %f\n", om0, oml0, hubble, wp5, ocurv);
    
    for(i=0;i<100;i++){
      //    printf("extras %f wspecies %f\n", extras[i], wspecies[0]);
    }
    if(extras[99]>0.)
      box =extras[99];
    else
      {
	puts("needs box size in comoving Mpc/h =");
	exit(1);
      }
    printf("box [Mpc/h]%f\n", box);
    
    if(i_file==my_snap){/*allocate only if it is the first file*/
      /* SetWeights for all particles (multi masses possible)*/
      if(nspecies==0)/*old constant weigths*/
	{
	  Nparticles = (NROW*NROW*NROW);    
	  if(!(P->Pos = malloc(3*Nparticles*sizeof(float))))
	    {
	      printf("problem allocating memory for particles positions (%lld particles)\n", Nparticles);
	      exit(1);
	    }
	  if(!(P->Vel = malloc(3*Nparticles*sizeof(float))))
	    {
	      printf("problem allocating memory for particles velocities\n");
	      exit(1);
	    }
	  if(!(P->Mass = malloc(sizeof(float))))
	    {
	      printf("problem allocating memory for particles masses\n");
	      exit(1);
	    }
	  ww = 1.0;
	}
      else
	{
	  Nparticles = lspecies[nspecies-1];



     
	  xscale= box/NGRID ;
	  vscale= 100.0*xscale/aexpn;
	  
	  if(nspecies==0)
	    {
	      npages =NROW;
	      n_in_last = NPAGE;       
	    }
	  else
	    {
	      Nparticles = lspecies[nspecies-1];
	      npages = (Nparticles - 1)/NPAGE  + 1 ;
	      n_in_last = Nparticles - NPAGE*(npages - 1);
	    }

	  
	  n_p_f = npages/n_files; /*number of pages per file*/
	  if(i_file==0){
	    ntot = 0;
	  }	  
	  if(i_file*n_p_f<npages && i_file<n_files)
	    {
	      i_in_page = NPAGE;
	    }
	  else
	    {
	      i_in_page = n_in_last;
	    }		  
	  printf("Nparticles %lld pages %d n_p_f %d i_in_page %d\n", Nparticles, npages, n_p_f, i_in_page);
	  
	  Nparticles = i_in_page * npages;

	  printf("Nparticles %lld pages %d n_p_f %d i_in_page %d\n", Nparticles, npages, n_p_f, i_in_page);

	  jstart = 1;
	  if(Nparticles<0)
	    {
	      printf("wrong number of particles");
	      printf("Npart %d Nspecies %d N %d ", Nparticles, nspecies,(int)(extras[9 + nspecies]));
	      exit(1);
	    }
	  
	  if(!(P->Pos = malloc(3*Nparticles*sizeof(float))))
	    {
	      printf("problem allocating memory for particles positions (%lld particles)\n", Nparticles);
	      exit(1);
	    }
	  if(!(P->Vel = malloc(3*Nparticles*sizeof(float))))
	    {
	      printf("problem allocating memory for particles velocities\n");
	      exit(1);
	    }
	  if(!(P->Mass = malloc(sizeof(float))))
	    {
	      printf("problem allocating memory for particles masses\n");
	      exit(1);
	    }
	  
	  for(j=1;j<=nspecies;j++)
	    {
	      /*
		jend = lspecies[j-1];	    
		printf("%d %d\n", jstart, jend);
		for(k=jstart;k<=jend;k++)
		{
		P->Mass[k-1] = wspecies[j-1];
		}
		jstart = jend;
	      */
	      P->Mass[0] = wspecies[j-1];
	    }            
	}
    }
     
     printf("Nparticles %lld\n", Nparticles);
  



   printf("Pages %lld Species %d n_in_last %d \n", npages, nspecies, n_in_last);


   n_p_f = npages/n_files; /*number of pages per file*/

     ntot = 0;

   for(ipage=1;ipage<=n_p_f;ipage++)
     {
       if(i_file*n_p_f<npages && i_file<n_files)
	 {
	   i_in_page = NPAGE;
	 }
       else
	 {
	   i_in_page = n_in_last;
	 }	
       
       SKIP;
       fread(buffer, sizeof(FLOAT), i_in_page*6, datafile);       
       SKIP;
       printf("read page %d\n",ipage);fflush(0);

       for(i=0;i< i_in_page;i++){
	 //       if(!(i%100000))printf("%d %d %lld\n", i, i_in_page, (i + 3*NPAGE));fflush(0);

#ifdef SWAP 
	 if(sizeof(FLOAT)==4){
	   P->Pos[3*(ntot + i) + 0] = (swapF(buffer[i    ])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 1] = (swapF(buffer[i + NPAGE])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 2] = (swapF(buffer[i + 2*NPAGE])-1.0)*xscale*1000.0;
	   P->Vel[3*(ntot + i) + 0] = swapF(buffer[i + 3*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 1] = swapF(buffer[i + 4*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 2] = swapF(buffer[i + 5*NPAGE])*vscale;
	 }else{
	   P->Pos[3*(ntot + i) + 0] = (swapD(buffer[i    ])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 1] = (swapD(buffer[i + NPAGE])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 2] = (swapD(buffer[i + 2*NPAGE])-1.0)*xscale*1000.0;
	   P->Vel[3*(ntot + i) + 0] = swapD(buffer[i + 3*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 1] = swapD(buffer[i + 4*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 2] = swapD(buffer[i + 5*NPAGE])*vscale;
	 }
#else
	   P->Pos[3*(ntot + i) + 0] = ((buffer[i    ])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 1] = ((buffer[i + NPAGE])-1.0)*xscale*1000.0;
	   P->Pos[3*(ntot + i) + 2] = ((buffer[i + 2*NPAGE])-1.0)*xscale*1000.0;
	   P->Vel[3*(ntot + i) + 0] = (buffer[i + 3*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 1] = (buffer[i + 4*NPAGE])*vscale;
	   P->Vel[3*(ntot + i) + 2] = (buffer[i + 5*NPAGE])*vscale;
#endif
	   //	   if(!(i%100000))fprintf(stdout, "%f %f %f\n", P->Pos[3*(ntot+i) +0],
	   //				  P->Pos[3*(ntot+i) +1],
	   //				  P->Pos[3*(ntot+i) +2]);
       }

       ntot = ntot + i_in_page;                 
       printf("read page %d, ntot %lld outof %lld \n",ipage, ntot,Nparticles);fflush(0);    
     }
   printf("finished reading file %d and %lld particles out of %lld\n", i_file, ntot, Nparticles);



   fclose(datafile);
   fclose(filein);
     }
  }
   //   ntot = Nparticles;
         



  memset(P->header.npart,0,6*sizeof(int));
  memset(P->header.npartTotal,0,6*sizeof(int));
  memset(P->header.mass,0,6*sizeof(double));

  P->header.time = aexpn/1.0;
  P->header.redshift = aexpn;
  P->N = P->header.npart[1] = P->header.npartTotal[1] = ntot;
  P->header.mass[1]=27.8782*pow(box,3)/Nparticles
      *pow(hubble/100.,3.)*om0;             //unit: 10^(10) Msol/h
  P->header.BoxSize = box*1000.0;   //box size in kpc/h
  P->header.Omega0 = om0;
  P->header.OmegaLambda = oml0;
  P->header.HubbleParam = hubble;



   printf("starts the test\n");
   /* test */
   float minx=1.E10, miny=1.E10, minz=1.E10, minvx=1.E10, minvy=1.E10, minvz=1.E10;
   float maxx=-1.E10, maxy=-1.E10, maxz=-1.E10, maxvx=-1.E10, maxvy=-1.E10, maxvz=-1.E10;
   for(i=0;i<ntot;i++)
     {
       minx = MIN(minx,P->Pos[3*i + 0]);
       miny = MIN(miny,P->Pos[3*i + 1]);
       minz = MIN(minz,P->Pos[3*i + 2]);
       minvx = MIN(minvx,P->Vel[3*i + 0]);
       minvy = MIN(minvy,P->Vel[3*i + 1]);
       minvz = MIN(minvz,P->Vel[3*i + 2]);

       maxx = MAX(maxx,P->Pos[3*i + 0]);
       maxy = MAX(maxy,P->Pos[3*i + 1]);
       maxz = MAX(maxz,P->Pos[3*i + 2]);
       maxvx = MAX(maxvx,P->Vel[3*i + 0]);
       maxvy = MAX(maxvy,P->Vel[3*i + 1]);
       maxvz = MAX(maxvz,P->Vel[3*i + 2]);
     }   
   
   printf("xmin, xmax = %g %g h-1 Mpc\n", minx, maxx);
   printf("ymin, ymax = %g %g h-1 Mpc\n", miny, maxy);
   printf("zmin, zmax = %g %g h-1 Mpc\n", minz, maxz);

   printf("vxmin, vxmax = %g %g km/s \n", minvx, maxvx);
   printf("vymin, vymax = %g %g km/s\n", minvy, maxvy);
   printf("vzmin, vzmax = %g %g km/s\n", minvz, maxvz);



  return ntot; 
}

