#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include "io.h"
#include "READ_ART.h"

/* Default value for line length.  */
static const int line_size = 1024;

ssize_t 
Mygetline (char **lineptr, int *n, FILE *stream)
{
  return getdelim (lineptr, n, '\n', stream);
}

//checks if file s an arepo file
int IsHDF5File(char *fname){
    
    hid_t file;
    //check if file opens
    if((file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT))<0){
      return 1;
    }else{
      return 0;
    }
}


//checks if file is a gadget file
int IsGadgetFile(char *fname)
{
    FILE *fd;
    int i;
   
    //check if file exists
    if(!(fd=fopen(fname,"r")))
	return 1;

    fread(&i,sizeof(int),1,fd);


#ifdef COMP
    if ((i!=8)&&(swapI(i)!=8))
	return 0;
#else
    if ((i!=256)&&(swapI(i)!=256))
	return 0;
#endif

    fclose(fd);

    return 1;
}

/*
This reads the ART snapshot and puts it into the adequate sstructure used in the rest of the code;

*/
int ReadART(char *fname_struct, char *fname_data, snapshot_data *P, int nfiles, int my_file)
{
    int N_part;
    if(nfiles>0){
      N_part = read_art_snap_multifile(fname_struct, fname_data, P, nfiles, my_file);
    }else{
      N_part = read_art_snap(fname_struct, fname_data, P);
    }
    return 0;
}

/*
  reads full simulation in hdf5 format, asuming the structure in arepo
  simulations
*/
int ReadHDF5File(char *fname, snapshot_data *P, int flags){
  int nfiles;
  int n_particles_this_file;
  long long n_particles_total[6];
  long long n_particles_file[6];
  double masses[6];
  char filein[1024];
  hid_t file, group, attr, status_n, space, dset, filespace, memspace;
  int ndim;
  float *tmp_data;
  hsize_t dims[2];
  int rank;
  int i;
  int ifile;
  long long n_items;
  long long i_part;
  long long total_part;

  fprintf(stdout, "Hello, deja el show\n");  
  sprintf(filein, "%s.0.hdf5", fname);

  file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
  group = H5Gopen(file, "/Header", H5P_DEFAULT);
  attr = H5Aopen(group, "NumFilesPerSnapshot", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_INT, &nfiles);
  attr = H5Aopen(group, "NumPart_Total", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_LONG, n_particles_total);
  attr = H5Aopen(group, "NumPart_ThisFile", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_LONG, n_particles_file);
  attr = H5Aopen(group, "MassTable", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, masses);
  attr = H5Aopen(group, "Time", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, &(P->header.time));
  attr = H5Aopen(group, "Redshift", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, &(P->header.redshift));
  attr = H5Aopen(group, "BoxSize", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, &(P->header.BoxSize));
  attr = H5Aopen(group, "Flag_Sfr", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_INT, &(P->header.flag_sfr));
  attr = H5Aopen(group, "Flag_Feedback", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_INT, &(P->header.flag_feedback));
  attr = H5Aopen(group, "Flag_Cooling", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_INT, &(P->header.flag_cooling));
  attr = H5Aopen(group, "Omega0", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, &(P->header.Omega0));
  attr = H5Aopen(group, "OmegaLambda", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, &(P->header.OmegaLambda));
  attr = H5Aopen(group, "HubbleParam", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_DOUBLE, &(P->header.HubbleParam));

  H5Aclose(attr);
  H5Gclose(group);
  H5Fclose(file);

  for(i=0;i<6;i++){
    fprintf(stdout, "\n");
    fprintf(stdout, "Npart Total %lld\n", n_particles_total[i]);
    fprintf(stdout, "Npart File %lld\n", n_particles_file[i]);
    fprintf(stdout, "Mass Table %f\n", masses[i]);
    P->header.npart[i] = n_particles_file[i];
    P->header.npartTotal[i] = n_particles_total[i];
    P->header.mass[i] = masses[i];
  }
  P->header.num_files = nfiles;

  fprintf(stdout, "The total number of files to read:%d\n", nfiles);
  fprintf(stdout, "%.2f millions of particles to allocate\n", n_particles_total[DM_TYPE]/1.0E6);
  fprintf(stdout, "%.2f MB to allocate\n", 6.0*sizeof(float)*n_particles_total[DM_TYPE]/(1024.0*1024.0));

  /* now go for data allocation. Only DM.*/
  if(flags&FLAG_POS){
    if(!(P->Pos = (float *)realloc(P->Pos,3*n_particles_total[DM_TYPE]*sizeof(float))))
      {
	fprintf(stderr,"failed to allocate memory pos.\n");
	exit(0);
      }
    fprintf(stdout, "realloc'd pos\n");
  }
  
  if (flags&FLAG_VEL){
    if(!(P->Vel = (float *)realloc(P->Vel,3*n_particles_total[DM_TYPE]*sizeof(float))))
      {
	fprintf(stderr,"failed to allocate memory vel.\n");
	exit(0);
      }
  }
  
  total_part = 0;  
  /*loop over all files and fill the data*/
  for(ifile=0;ifile<nfiles;ifile++){
    sprintf(filein, "%s.%d.hdf5", fname, ifile);
    file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    /*Get the number of particles for this file*/
    group = H5Gopen(file, "/Header", H5P_DEFAULT);
    attr = H5Aopen(group, "NumPart_ThisFile", H5P_DEFAULT);
    status_n = H5Aread(attr, H5T_NATIVE_LONG, n_particles_file);

    /*Read all the positions*/
    if(flags&FLAG_POS){
      dset = H5Dopen2(file, DM_COORDINATES,H5P_DEFAULT);
      fprintf(stdout, "getting dataset %s\n", DM_COORDINATES);
      filespace = H5Dget_space (dset);
      
      rank = H5Sget_simple_extent_ndims (filespace);
      
      status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
      fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));
      
    
      n_items = (long long)dims[0] * (long long)dims[1];
      if(!(tmp_data=malloc(sizeof(float) * n_items))){
	fprintf(stderr, "problem with data allocation\n");
      }
      memspace = H5Screate_simple (rank,dims,NULL);  
      status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
			  H5P_DEFAULT, tmp_data);
      
      if(n_particles_file[DM_TYPE]!=dims[0]){
	fprintf(stderr, "Inconsistent number of particles\n");
	exit(1);
      }
   
      fprintf(stdout, "total part reading %lld\n", total_part);
      for(i_part=0;i_part<n_particles_file[DM_TYPE];i_part++){
	P->Pos[(total_part*3) + i_part*3 + 0] = tmp_data[i_part*3 + 0];
	P->Pos[(total_part*3) + i_part*3 + 1] = tmp_data[i_part*3 + 1];
	P->Pos[(total_part*3) + i_part*3 + 2] = tmp_data[i_part*3 + 2];
      }
      free(tmp_data);
    }
    
    /*Read all the velocities*/
    if(flags&FLAG_VEL){
      dset = H5Dopen2(file, DM_COORDINATES,H5P_DEFAULT);
      fprintf(stdout, "getting dataset %s\n", DM_VELOCITIES);
      filespace = H5Dget_space (dset);
      
      rank = H5Sget_simple_extent_ndims (filespace);
      
      status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
      fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));
      
    
      n_items = (long long)dims[0] * (long long)dims[1];
      if(!(tmp_data=malloc(sizeof(float) * n_items))){
	fprintf(stderr, "problem with data allocation\n");
      }
      
      memspace = H5Screate_simple (rank,dims,NULL);  
      status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
			  H5P_DEFAULT, tmp_data);
      
      if(n_particles_file[DM_TYPE]!=dims[0]){
	fprintf(stderr, "Inconsistent number of particles\n");
	exit(1);
      }
           
      for(i_part=0;i_part<n_particles_file[DM_TYPE];i_part++){
	P->Vel[(total_part*3) + i_part*3 + 0] = tmp_data[i_part*3 + 0];
	P->Vel[(total_part*3) + i_part*3 + 1] = tmp_data[i_part*3 + 1];
	P->Vel[(total_part*3) + i_part*3 + 2] = tmp_data[i_part*3 + 2];
      }
      free(tmp_data);
    }
    
    total_part  += n_particles_file[DM_TYPE];    
    H5Fclose(file);        
  }
  P->N = total_part;
  
  
  //  exit(1);
  return 0;
}


/* POINTERS in P MUST BE SET TO NULL BEFORE CALLING */
/* OR TO PREVIOUSLY ALLOCATED DATA  */
/* SHOULD CGANGE REALLOC-> MALLOC OR CALLOC */
/*
 * reads gadget format (multiple files) 
 * if filename is "" stdin is read
 */
int ReadGadget(char *fname, snapshot_data *P, int flags)
{
    FILE *fd;
    char buf[200];
    char name[4];
    int i,k,ntot_withmasses, ifile;
    long long  n,pc,pc_new,pc_sph;
#ifdef JUANK
    long long my_pc;
#endif
    int Ngas;
    long long NumPart=0;
    int files;
    long long i_part;
    long long n_part;

    //guess how many files there are
    sprintf(buf,"%s",fname);
    if(!(fd=fopen(buf,"r")))
    {
	sprintf(buf,"%s.%d",fname,0);
	files=0;
	while ((fd=fopen(buf,"r")))
	{
	    fclose(fd);
	    sprintf(buf,"%s.%d",fname,++files);
	}
	if(files<2)
	{
	    printf("can't open file '%s' nor `%s`\n",fname,buf);
	    exit(0);
	}
    }
    else 
    {
	files=1;
	fclose(fd);
    }
    fprintf(stdout, "files %d\n", files);
    for(ifile=0, pc=0; ifile<files; ifile++, pc=pc_new)
    {
	if(files>1)
	    sprintf(buf,"%s.%d",fname,ifile);
	else
	    sprintf(buf,"%s",fname);
	
	if (strlen(buf)==0) fd=stdin;
	else
	    if(!(fd=fopen(buf,"r")))
	    {
		printf("can't open file `%s`\n",buf);
		exit(0);
	    }
	


	printf("reading %s ...",buf); fflush(stdout);

#ifdef NEWFORMAT
 {
	    int skpdummy;
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    fread(name,sizeof(char),4,fd);
	    fprintf(stdout, "Name %s\n", name);
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
 }
#endif	
	//SKIP(fd);
	//Guess if endianness changes or not 
	{
	    int skpdummy;
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    fprintf(stdout, "dummy %d\n", skpdummy);
	    if (skpdummy==256) 
	    {
		flags &= (~FLAG_SWAPENDIAN);
	    }
		else
	    {
		flags |= FLAG_SWAPENDIAN;
		printf ("(Endianness differs)");fflush(0);
	    }
	}

         fread(&P->header, sizeof(P->header), 1, fd);
 

	 for (i=0;i<6;i++){
	   printf ("Size=%d Header Mass=%f\n",P->header.npart[i], P->header.mass[i]);
	 }
        //eventually swap endianness
	if (flags&FLAG_SWAPENDIAN)
	{
	    Dswap4BArr(P->header.npart,6);
	    Dswap4BArr(P->header.npartTotal,6);
	    Dswap8BArr(P->header.mass,6);

	    Dswap4B(&P->header.flag_sfr);Dswap4B(&P->header.flag_feedback);Dswap4B(&P->header.flag_cooling);
	    Dswap4B(&P->header.num_files);
	    Dswap8B(&P->header.time);Dswap8B(&P->header.redshift);Dswap8B(&P->header.BoxSize);
	    Dswap8B(&P->header.Omega0);Dswap8B(&P->header.OmegaLambda);Dswap8B(&P->header.HubbleParam);
	}
	SKIP(fd);

	fprintf(stdout, "redshift %e\n", P->header.redshift);
	fprintf(stdout, "num files %d\n", P->header.num_files);
	fprintf(stdout, "Omega_0 %g\n", P->header.Omega0);

	if(files==1)
	{
	    for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
		NumPart+= P->header.npart[k];
	    Ngas= P->header.npart[0];
	}
	else
	{
	    for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
		NumPart+= P->header.npartTotal[k];
	    Ngas= P->header.npartTotal[0];
	}
	
      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(P->header.mass[k]==0)
	      ntot_withmasses+= P->header.npart[k];
	}

      P->N = NumPart;

#ifdef NEWFORMAT
 {
	    int skpdummy;
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    fread(name,sizeof(char),4,fd);
	    fprintf(stdout, "Name %s\n", name);
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
 }
#endif	      
 
      fprintf(stdout, "Particle number %lld %lld %lld \n", 3*NumPart, NumPart, (long long) P->header.npart[0]);
      fprintf(stdout, "N tot with masses %d\n", ntot_withmasses);
      if(ifile==0)
      {
	  if(flags&FLAG_POS){
	      if(!(P->Pos=(float *)realloc(P->Pos,3*NumPart*sizeof(float))))
	      {
		  fprintf(stderr,"failed to allocate memory pos.\n");
		  exit(0);
	      }
	      fprintf(stdout, "realloc'd pos\n");
	  }
	  
	  if (flags&FLAG_VEL)
	      if(!(P->Vel=(float *)realloc(P->Vel,3*NumPart*sizeof(float))))
	      {
		  fprintf(stderr,"failed to allocate memory vel.\n");
		  exit(0);
	      }

	  if (flags&FLAG_ID)
	      if(!(P->Id=(int *)realloc(P->Id,NumPart*sizeof(int))))
	      {
		  fprintf(stderr,"failed to allocate memory id.\n");
		  exit(0);
	      }

	  if (flags&FLAG_TYPE)
	      if(!(P->Type=(char *)realloc(P->Type,NumPart*sizeof(char))))
	      {
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
	      }
	  

	  if (P->header.npartTotal[0])
	  {
	      if (flags&FLAG_GAS)
	      {
		  if(!(P->Rho=(float *)realloc(P->Rho,P->header.npartTotal[0]*sizeof(float))))
		  {
		      fprintf(stderr,"failed to allocate memory.\n");
		      exit(0);
		  }
		  
		  if(!(P->U=(float *)realloc(P->U,P->header.npartTotal[0]*sizeof(float))))
		  {
		      fprintf(stderr,"failed to allocate memory.\n");
		      exit(0);
		  }
		  
		  if(!(P->Temp=(float *)realloc(P->Temp,P->header.npartTotal[0]*sizeof(float))))
		  {
		      fprintf(stderr,"failed to allocate memory.\n");
		      exit(0);
		  }
		  
		  if(!(P->Ne=(float *)realloc(P->Ne,P->header.npartTotal[0]*sizeof(float))))
		  {
		      fprintf(stderr,"failed to allocate memory.\n");
		      exit(0);
		  }
	      }
	  }

	  if ((ntot_withmasses)&&(flags&FLAG_MASS)){
	    fprintf(stderr, "Allocating masses (ntot_withmasses %d\n)", ntot_withmasses);
	    if(!(P->Mass=(float *)realloc(P->Mass,NumPart*sizeof(float))))
	      {
		fprintf(stderr,"failed to allocate memory.\n");
		exit(0);
	      }	   
	    
	  }else{
	    if(flags&FLAG_MASS){
	      fprintf(stderr, "Allocating masses (ntot_withmasses %d) for constant values\n", ntot_withmasses);
	      if(!(P->Mass=(float *)realloc(P->Mass,NumPart*sizeof(float))))
		{
		  fprintf(stderr,"failed to allocate memory.\n");
		  exit(0);
		}	   
	      
	      n_part = 0;
	      for(k=0;k<6;k++){	      
		for(i_part=n_part;i_part<(n_part + P->header.npart[k]);i_part++){
		  P->Mass[i_part] = P->header.mass[k];
		}
		n_part = i_part;
		fprintf(stdout, "Filled with mass %f down to %lld particles now\n", P->header.mass[k], n_part);
	      }
	      fprintf(stdout, "Filled masses for %lld particles\n", n_part);
	    }
	  }
      }
	  
      SKIP(fd);
      
#ifdef JUANK
      my_pc = 0;
#endif
      for(k=0;k<6;k++)
      {

	  if (flags&FLAG_POS && P->header.npart[k])
	  {
#ifdef JUANK
	    
	      fread(&P->Pos[3*my_pc], sizeof(float), 3*P->header.npart[k], fd);
	      fprintf(stdout, "read now %lld, accumulated %lld %e\n", P->header.npart[k], 3*my_pc, P->Pos[3*my_pc]);
	      if (flags&FLAG_SWAPENDIAN)
		  Dswap4BArr(&P->Pos[3*pc],3*P->header.npart[k]);
	      my_pc = my_pc + P->header.npart[k];
#else	      	      
	      fread(&P->Pos[3*pc], sizeof(float), 3*P->header.npart[k], fd);
	      fprintf(stdout, "read now %lld, accumulated %lld %e\n", P->header.npart[k], 3*pc, P->Pos[3*pc]);
	      if (flags&FLAG_SWAPENDIAN)
		  Dswap4BArr(&P->Pos[3*pc],3*P->header.npart[k]);	      
#endif
	  }
	  else
	      fseek(fd,3*P->header.npart[k]*sizeof(float),SEEK_CUR);
      }
      SKIP(fd);

#ifdef NEWFORMAT
 {
	    int skpdummy;
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
	    fread(name,sizeof(char),4,fd);
	    fprintf(stdout, "Name %s\n", name);
	    fread(&skpdummy,sizeof(skpdummy),1,fd);
//	    fread(&skpdummy,sizeof(skpdummy),1,fd);
 }
#endif	      

      SKIP(fd);
      for(k=0;k<6;k++)
      {
	  if (flags&FLAG_VEL)
	  {
	      fread(&P->Vel[3*pc], sizeof(float), 3*P->header.npart[k], fd);
	      if (flags&FLAG_SWAPENDIAN)
		  Dswap4BArr(&P->Vel[3*pc],3*P->header.npart[k]);
	  }
	  else
	      fseek(fd,3*P->header.npart[k]*sizeof(float),SEEK_CUR);
      }
      SKIP(fd);

      SKIP(fd);
      for(k=0;k<6;k++)
      {
	  if (flags&FLAG_ID)
	  {
	      //printf ("pc=%d\n",pc);
	      fread(&P->Id[pc], sizeof(int), P->header.npart[k], fd);
	      if (flags&FLAG_SWAPENDIAN)
		  Dswap4BArr(&P->Id[pc],P->header.npart[k]);
	  }
	  else
	      fseek(fd,P->header.npart[k]*sizeof(int),SEEK_CUR);
      }
      SKIP(fd);

      if(ntot_withmasses>0)
	  SKIP(fd);

      for(k=0, pc_new=pc; k<6; k++)
      {
	  for(n=0;n<P->header.npart[k];n++)
	  {
	      if (flags&FLAG_TYPE)
		  P->Type[pc_new]=(char)k;
	      
	      if (flags&FLAG_MASS)
	      {
		  if (P->header.mass[k]==0)
		  {
		    fprintf(stdout, "READING MASS\n");
		      fread(&P->Mass[pc_new], sizeof(float), 1, fd);
		      if (flags&FLAG_SWAPENDIAN)
			  Dswap4B(&P->Mass[pc_new]);
		  }
		  else if(ntot_withmasses>0)
		      P->Mass[pc_new]= P->header.mass[k];
	      }
	      else if (P->header.mass[k]==0)
		  fseek(fd,sizeof(float),SEEK_CUR);

	      pc_new++;
	  }
      }
      
      if(ntot_withmasses>0)
	  SKIP(fd);

      if ((P->header.npart[0]>0)&&(flags&FLAG_GAS))
	{
	    SKIP(fd);
	    fread(&P->U[pc], sizeof(float), P->header.npart[0], fd);
	    if (flags&FLAG_SWAPENDIAN)
		Dswap4BArr(&P->U[pc],P->header.npart[0]);
	    SKIP(fd);
	    
	    SKIP(fd);
	    fread(&P->Rho[pc], sizeof(float), P->header.npart[0], fd);
	    if (flags&FLAG_SWAPENDIAN)
		Dswap4BArr(&P->Rho[pc],P->header.npart[0]);
	    SKIP(fd);

	    if(P->header.flag_cooling)
	    {
		SKIP(fd);
		fread(&P->Ne[pc], sizeof(float), P->header.npart[0], fd);
		if (flags&FLAG_SWAPENDIAN)
		    Dswap4BArr(&P->Ne[pc],P->header.npart[0]);
		SKIP(fd);
	    }
	    else
		for(n=0, pc_sph=pc; n<P->header.npart[0];n++)
		{
		    P->Ne[pc_sph]= 1.0;
		    pc_sph++;
		}
	}
      
      fclose(fd);
      printf(" done.\n",buf); fflush(stdout);
    }
}

/* Reads SIMPLE format to a gadget structure */
/* Can only read position and velocities so far */
int ReadSIMPLE2Gadget(char *fname, snapshot_data *P, int flags)
{

    struct {
	int n;
	float tnow;
	float alpha;
	float massp;
	float aexp;
	float Omega0;
	float Omegat;
	float OmegaLt;
	float OmegaL0;
	float hubble;
    } SMP_header;

    struct {
	float ai;// a initial
	float Omegaf;//Omega matter final
	float OmegaLf;// Omega Lambda final
	float af;//a final
	float Lf;//Final length of the box in Mpc
	float Hf;//final hubble constant
    } input_dec;
    
    typedef long long INT8B;
    INT8B ikey;

    int i,j,k,l;
    float temp[(1<<18)*3];
    FILE *fd;
    char buf[256];
    char buf2[256];
    char *line = NULL;

    //First read cosmo parmaeters from input_dec.dat, same directory
    
    for (k=0,i=0;k<strlen(fname);k++) if (fname[k]=='/') i=k;
    strcpy(buf,fname);
    if (i==0) buf[0]='\0'; else buf[i+1]='\0';
    sprintf(buf,"%s%s",buf,"input_dec.dat");
    
    if(!(fd=fopen(buf,"r")))
    {
	printf("can't open file `%s`\n",buf);fflush(0);
	exit(0);
    }
    memset(&input_dec,0,sizeof(input_dec));
    input_dec.ai=1.;
    input_dec.Omegaf=0.3333;
    input_dec.OmegaLf=0.6667;
    while (Mygetline(&line,&i,fd)!=-1)
    {
	sscanf(line,"%s = %s",buf,buf2);
	if (!strcmp(buf,"af")) input_dec.af=atof(buf2);
	if (!strcmp(buf,"ai")) input_dec.ai=atof(buf2);
	if (!strcmp(buf,"H_f")) input_dec.Hf=atof(buf2);
	if (!strcmp(buf,"lbox")) input_dec.Lf=atof(buf2);
	if (!strcmp(buf,"omega_f")) input_dec.Omegaf=atof(buf2);
	if (!strcmp(buf,"lambda_f")) input_dec.OmegaLf=atof(buf2);
    }
    if ((input_dec.af==0)||(input_dec.ai==0)||(input_dec.Hf==0)||(input_dec.Lf==0))
    {
	printf("Missing parameters in input_dec.dat\n");
    }
    
    fclose(fd);
    
    if(!(fd=fopen(fname,"r")))
    {
	printf("can't open file `%s`\n",fname);fflush(0);
	exit(0);
    }

    printf("reading %s ...",fname); fflush(0);
    
    //Guess if endianness changes or not 
    {
	int skpdummy;
	fread(&skpdummy,sizeof(skpdummy),1,fd);
	if (skpdummy==40) 
	{
	    flags &= (~FLAG_SWAPENDIAN);
	}
	else
	{
	    flags |= FLAG_SWAPENDIAN;
	    printf ("(Endianness differs)");fflush(0);
	}
    }

    fread(&ikey,sizeof(INT8B),1,fd);
    if (flags&FLAG_SWAPENDIAN) Dswap8BArr(&ikey,1);

    SKIP(fd);

    if (ikey>0)
    {
	ikey =-110;
	fseek(fd,4,SEEK_SET);
	printf ("(Old file format)");fflush(0);
    }

    fread (&SMP_header,sizeof(SMP_header),1,fd);
    if (flags&FLAG_SWAPENDIAN) Dswap4BArr((int *)&(SMP_header),sizeof(SMP_header)/sizeof(int));

    SKIP(fd);
    memset(P->header.npart,0,6*sizeof(int));
    memset(P->header.npartTotal,0,6*sizeof(int));
    memset(P->header.mass,0,6*sizeof(double));

    //Got everything to fill Gadget header ...
    P->header.time=SMP_header.aexp/input_dec.af;
    P->header.redshift=1./P->header.time-1.;
    P->N=P->header.npart[1]=P->header.npartTotal[1]=SMP_header.n;
    P->header.mass[1]=27.8782*pow(input_dec.Lf,3)/SMP_header.n
	*pow(input_dec.Hf/100.,3.)*input_dec.Omegaf;             //unit: 10^(10) Msol/h
    P->header.BoxSize = input_dec.Hf/100. * input_dec.Lf * 1000.;   //box size in kpc/h
    P->header.Omega0 = input_dec.Omegaf;
    P->header.OmegaLambda = input_dec.OmegaLf;
    P->header.HubbleParam = input_dec.Hf/100.0;

    
    //Read the masses here
    //if (((INT8B)(-ikey)!=(INT8B)(-ikey/10)*10)&&(flags&FLAG_MASS))

    //now read data
    //Positons, ikey is not an int but we don t care here, only need velocity and position ... 
    if (((INT8B)(-ikey/10)!=(INT8B)(-ikey/100)*10)&&(flags&FLAG_POS))
    {

	if (!(P->Pos = (float *) malloc(3*sizeof(float)*SMP_header.n)))
	{
	    fprintf(stderr,"failed to allocate memory for Positions in ReadSIMPLE2Gadget.\n");
	    exit(0);
	}

	for (i=0,j=0;i<SMP_header.n/(1<<18);i++)
	{
	    SKIP(fd);
	    fread(temp,3*sizeof(float)*(1<<18),1,fd);
	    SKIP(fd);
	    
	    if (flags&FLAG_SWAPENDIAN)
	    {
		for (k=0;k<(1<<18);k++)
		    for (l=0;l<3;l++)
			P->Pos[3*(k+j)+l]=(swapF(temp[k+l*(1<<18)])+0.5)*P->header.BoxSize;
	    }
	    else
	    {    
		for (k=0;k<(1<<18);k++)
		    for (l=0;l<3;l++)
			P->Pos[3*(k+j)+l]=(temp[k+l*(1<<18)]+0.5)*P->header.BoxSize;
	    }

	    j+=(1<<18); 
	}
	//if number of particles is not a multiple of (1<<18)
	i=SMP_header.n%(1<<18);
	if (i)
	{
	    SKIP(fd);
	    fread(temp,3*sizeof(float)*i,1,fd);
	    SKIP(fd);
	    
	    if (flags&FLAG_SWAPENDIAN)
	    {
		for (k=0;k<i;k++)
		    for (l=0;l<3;l++)
			P->Pos[3*(k+j)+l]=(swapF(temp[k+l*(1<<18)])+0.5)*P->header.BoxSize;
	    }
	    else
	    {
		for (k=0;k<i;k++)
		    for (l=0;l<3;l++)
			P->Pos[3*(k+j)+l]=(temp[k+l*(1<<18)]+0.5)*P->header.BoxSize;
	    }
	    j+=i;
	}
	if (j!=SMP_header.n) 
	    printf("There is a bug in ReadSIMPLE2Gadget, only %d particules read over %d total",j,P->header.npart[1]);
    }
    else if ((INT8B)(-ikey/10)!=(INT8B)(-ikey/100)*10)
    {
	fseek(fd,SMP_header.n*3*sizeof(float),SEEK_CUR);
    }



    if (((INT8B)(-ikey/100)!=(INT8B)(-ikey/1000)*100)&&(flags&FLAG_VEL))
    {
	float v_fact;

	//factor to convert simple 2 gadget units (km/s/sqrt(a))
	v_fact = P->header.BoxSize*P->header.HubbleParam/sqrt(P->header.time)* 
	    sqrt(P->header.Omega0*pow(P->header.time,3)+P->header.OmegaLambda);
	    
	
	//Velocities, have the right conversion !!!!
	if (!(P->Vel = (float *) malloc(3*sizeof(float)*SMP_header.n)))
	{
	    fprintf(stderr,"failed to allocate memory for Velocities in ReadSIMPLE2Gadget.\n");
	    exit(0);
	}

	

	for (i=0,j=0;i<SMP_header.n/(1<<18);i++)
	{
	    SKIP(fd);
	    fread(temp,3*sizeof(float)*(1<<18),1,fd);
	    SKIP(fd);
	    
	    if (flags&FLAG_SWAPENDIAN)
	    {
		for (k=0;k<(1<<18);k++)
		    for (l=0;l<3;l++)
			P->Vel[3*(k+j)+l]=v_fact*swapF(temp[k+l*(1<<18)]);
	    }
	    else
	    {
		for (k=0;k<(1<<18);k++)
		    for (l=0;l<3;l++)
			P->Vel[3*(k+j)+l]=v_fact*temp[k+l*(1<<18)];
	    }

	    j+=(1<<18); 
	}
	i=SMP_header.n%(1<<18);
	if (i)
	{
	    SKIP(fd);
	    fread(temp,3*sizeof(float)*i,1,fd);
	    SKIP(fd);
	    
	    if (flags&FLAG_SWAPENDIAN)
	    {
		for (k=0;k<i;k++)
		    for (l=0;l<3;l++)
			P->Vel[3*(k+j)+l]=swapF(temp[k+l*(1<<18)]);
	    }
	    else
	    {
		for (k=0;k<i;k++)
		    for (l=0;l<3;l++)
			P->Vel[3*(k+j)+l]=temp[k+l*(1<<18)];
	    }

	    j+=i;
	}
	if (j!=SMP_header.n) 
	    printf("There is a bug in ReadSIMPLE2Gadget, only %d particules read over %d total",j,P->header.npart[1]);
    }
    else if ((INT8B)(-ikey/100)!=(INT8B)(-ikey/1000)*100)
    {
	fseek(fd,SMP_header.n*3*sizeof(float),SEEK_CUR);
    }

    fclose(fd);
    
    printf(" done.\n",buf); fflush(stdout);

    return 0;
}

int LoadSurveyCone(SurveyCone *Cone,const char *fname)
{
  FILE *f;
  int i;
  int j;

  printf("Loading file %s ...",fname);fflush(0);

  if ((f=fopen(fname,"r"))==NULL)
    {
      fprintf(stderr,"\nCould not open \"%s\" for reading in \"LoadSurveyCone\"\n",fname);
      exit(0);
      return 1;
    }
  
  fread(&i,sizeof(int),1,f);
  fread(&Cone->N,sizeof(int),1,f);
  fread(&Cone->al,sizeof(int),1,f);
  fread(&Cone->ab,sizeof(int),1,f);
  fread(&Cone->v,sizeof(int),1,f);
  fread(&Cone->z,sizeof(int),1,f);
  fread(&Cone->d,sizeof(int),1,f);
  for (i=0;i<10;i++)
    {
	fread(&Cone->absmag[i],sizeof(int),1,f);
	fread(&Cone->appmag[i],sizeof(int),1,f);
    }
  fread(&Cone->d,sizeof(int),1,f);
  fread(&Cone->h,sizeof(float),1,f);
  for (i=0;i<10;i++)
      fread(Cone->filterName[i],sizeof(char)*30,1,f);
  fread(&Cone->NFilters,sizeof(int),1,f);
  fread(Cone->fromfile,sizeof(char)*256,1,f);
  fread(&Cone->NRows,sizeof(int),1,f);
  fread(&i,sizeof(int),1,f);

  printf("(N=%d)",Cone->N);

  Cone->Data=(float **)malloc(sizeof(float*)*Cone->NRows);
  for (j=0;j<Cone->NRows;j++)
    if ((Cone->Data[j]=(float *)malloc(sizeof(float)*Cone->N))==NULL)
      {
	fprintf(stderr,"\nCould not allocate Cone->data to load %s.\n",fname);
	exit(0);
      }

  fread(&i,sizeof(int),1,f);
  for (j=0;j<Cone->NRows;j++)
    fread(Cone->Data[j],sizeof(float),Cone->N,f);
  fread(&i,sizeof(int),1,f);
  
  if ((Cone->Pos=(float *)malloc(sizeof(float)*3*Cone->N))==NULL)
    {
      fprintf(stderr,"\nCould not allocate Cone->Pos to load %s.\n",fname);
      exit(0);
    }
  
  fread(&i,sizeof(int),1,f);
  fread(Cone->Pos,sizeof(float),3*Cone->N,f);
  fread(&i,sizeof(int),1,f);
  
  fclose (f);
  
  printf(" done.\n");
  
  return 0;
}

int LoadEigenvalueGrid(char *fname,density_grid *density)
{
    FILE *f;
    unsigned int i;

    char test[30];
    int swap = 0;
    
    //printf ("Loading density grid from %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }

    fread(&i,sizeof(int),1,f);
    fread(test,sizeof(char)*30,1,f);
    fread(&i,sizeof(int),1,f);

    if (i!=30*sizeof(char)) swap=1;

    if (strcmp(test,"Eigenvalue grid file")) 
	{
	    printf("-->");
	    return -1;
	}
    else printf ("Loading density grid from %s ... ",fname);fflush(0);

    if (swap) printf ("(Swapping)");fflush(0);
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->Nx),sizeof(int),1,f,swap);
    fread_sw(&(density->Ny),sizeof(int),1,f,swap);
    fread_sw(&(density->Nz),sizeof(int),1,f,swap);
    fread_sw(&(density->NNodes),sizeof(long long),1,f,swap);

    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);

    if (!(density->eigenvalue_1=malloc((size_t)density->NNodes*sizeof(float))))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->eigenvalue_1[0]),sizeof(float),density->NNodes,f,swap);
    fread(&i,sizeof(int),1,f);
    
    i=sizeof(float);
    fread(&i,sizeof(int),1,f);
    fread_sw(&density->redshift,sizeof(float),1,f,swap);
    fread_sw(&density->smoothing,sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);
    
/*
    for (i=0;i<density->NNodes;i++)
	density->grid[i]=density->grid[i]*10000+10;
    printf ("NNodes = %d\n",density->NNodes);
    for (i=0;i<100;i++)
    {
	j=i+0*(int)(density->NNodes*((double)rand()/RAND_MAX));
	printf ("density[%d] = %e\n",j,density->grid[j]);
    }
    printf ("dx=%f %f %f\n",density->dx,density->dy,density->dz);
*/
    //density->dx=density->dy=density->dz=1000;
    density->HasGradiant=0;
     printf ("done.\n");
    return 0;
}

int LoadEnvGrid(char *fname,density_grid *density)
{
    FILE *f;
    unsigned int i;

    char test[30];
    int swap = 0;
    
    //printf ("Loading density grid from %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }

    fread(&i,sizeof(int),1,f);
    fread(test,sizeof(char)*30,1,f);
    fread(&i,sizeof(int),1,f);

    if (i!=30*sizeof(char)) swap=1;

    if (strcmp(test,"Environment grid file")) 
	{
	    printf("-->");
	    return -1;
	}
    else printf ("Loading density grid from %s ... ",fname);fflush(0);

    if (swap) printf ("(Swapping)");fflush(0);
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->Nx),sizeof(int),1,f,swap);
    fread_sw(&(density->Ny),sizeof(int),1,f,swap);
    fread_sw(&(density->Nz),sizeof(int),1,f,swap);
    fread_sw(&(density->NNodes),sizeof(long long),1,f,swap);

    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);

    if (!(density->environment=(int *)malloc((size_t)density->NNodes*sizeof(FLOAT))))
    {
	fprintf(stderr,"Not enough memory for density->environment while loading\n");
	exit(0);
    }
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->environment[0]),sizeof(int),density->NNodes,f,swap);
    fread(&i,sizeof(int),1,f);
    
    i=sizeof(float);
    fread(&i,sizeof(int),1,f);
    fread_sw(&density->redshift,sizeof(float),1,f,swap);
    fread_sw(&density->smoothing,sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);
    
/*
    for (i=0;i<density->NNodes;i++)
	density->grid[i]=density->grid[i]*10000+10;
    printf ("NNodes = %d\n",density->NNodes);
    for (i=0;i<100;i++)
    {
	j=i+0*(int)(density->NNodes*((double)rand()/RAND_MAX));
	printf ("density[%d] = %e\n",j,density->grid[j]);
    }
    printf ("dx=%f %f %f\n",density->dx,density->dy,density->dz);
*/
    //density->dx=density->dy=density->dz=1000;
    density->HasGradiant=0;
     printf ("done.\n");
    return 0;
}



int LoadDensityGrid(char *fname, density_grid *density)
{
    FILE *f;
    unsigned int i;

    char test[30];
    int swap = 0;
    
    //printf ("Loading density grid from %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }

    fread(&i,sizeof(int),1,f);
    fread(test,sizeof(char)*30,1,f);
    fread(&i,sizeof(int),1,f);

    if (i!=30*sizeof(char)) swap=1;
 
    if (strcmp(test,"Density grid file")) 
	{
	    printf("-->");
	    return -1;
	}
    else printf ("Loading density grid from %s ... ",fname);fflush(0);

    if (swap) printf ("(Swapping)");fflush(0);
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->Nx),sizeof(int),1,f,swap);    
    fread_sw(&(density->Ny),sizeof(int),1,f,swap);
    fread_sw(&(density->Nz),sizeof(int),1,f,swap);
    fread_sw(&(density->NNodes),sizeof(long long),1,f,swap);

    fprintf(stderr, "Nx Ny Nz : %d %d %d\n", density->Nx, density->Ny, density->Nz);

    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);

    if (!(density->grid=(FLOAT *)malloc((size_t)density->NNodes*sizeof(FLOAT))))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->grid[0]),sizeof(FLOAT),density->NNodes,f,swap);
    fread(&i,sizeof(int),1,f);
    
    i=sizeof(float);
    fread(&i,sizeof(int),1,f);
    fread_sw(&density->redshift,sizeof(float),1,f,swap);
    fread_sw(&density->smoothing,sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);
    
/*
    for (i=0;i<density->NNodes;i++)
	density->grid[i]=density->grid[i]*10000+10;
    printf ("NNodes = %d\n",density->NNodes);
    for (i=0;i<100;i++)
    {
	j=i+0*(int)(density->NNodes*((double)rand()/RAND_MAX));
	printf ("density[%d] = %e\n",j,density->grid[j]);
    }
    printf ("dx=%f %f %f\n",density->dx,density->dy,density->dz);
*/
    //density->dx=density->dy=density->dz=1000;
    density->HasGradiant=0;
     printf ("done.\n");
    return 0;
}

int LoadVelocityGrid(char *fname, density_grid *density, int component)
{
    FILE *f;
    unsigned int i;
    long long n_to_alloc;
    char test[30];
    int swap = 0;
    
    //printf ("Loading density grid from %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }

    fread(&i,sizeof(int),1,f);
    fread(test,sizeof(char)*30,1,f);
    fread(&i,sizeof(int),1,f);

    if (i!=30*sizeof(char)) swap=1;
 
    if (strcmp(test,"Density grid file")) 
	{
	    printf("-->");
	    return -1;
	}
    else printf ("Loading velocity grid grid from %s ... ",fname);fflush(0);

    if (swap) printf ("(Swapping)");fflush(0);
    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->Nx),sizeof(int),1,f,swap);    
    fread_sw(&(density->Ny),sizeof(int),1,f,swap);
    fread_sw(&(density->Nz),sizeof(int),1,f,swap);
    fread_sw(&(density->NNodes),sizeof(long long),1,f,swap);

    fprintf(stderr, "Nx Ny Nz : %d %d %d\n", density->Nx, density->Ny, density->Nz);

    fread_sw(&(density->x0),sizeof(float),1,f,swap);
    fread_sw(&(density->y0),sizeof(float),1,f,swap);
    fread_sw(&(density->z0),sizeof(float),1,f,swap);
    
    fread_sw(&(density->dx),sizeof(float),1,f,swap);
    fread_sw(&(density->dy),sizeof(float),1,f,swap);
    fread_sw(&(density->dz),sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);

    n_to_alloc = density->NNodes*sizeof(FLOAT);
    if(component==0){
    if (!(density->grid_vx=(FLOAT *)malloc((size_t)n_to_alloc)))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->grid_vx[0]),sizeof(FLOAT),density->NNodes,f,swap);
    fread(&i,sizeof(int),1,f);
    }

    if(component==1){
    if (!(density->grid_vy=(FLOAT *)malloc((size_t)n_to_alloc)))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->grid_vy[0]),sizeof(FLOAT),density->NNodes,f,swap);
    fread(&i,sizeof(int),1,f);
    }

    if(component==2){
    if (!(density->grid_vz=(FLOAT *)malloc((size_t)n_to_alloc)))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }    
    fread(&i,sizeof(int),1,f);
    fread_sw(&(density->grid_vz[0]),sizeof(FLOAT),density->NNodes,f,swap);
    fread(&i,sizeof(int),1,f);
    }

    i=sizeof(float);
    fread(&i,sizeof(int),1,f);
    fread_sw(&density->redshift,sizeof(float),1,f,swap);
    fread_sw(&density->smoothing,sizeof(float),1,f,swap);
    fread(&i,sizeof(int),1,f);

    density->HasGradiant=0;
     printf ("done.\n");
    return 0;
}

int LoadDensity(char *fname,density_field *density)
{
    FILE *f;
    unsigned int i;
    int j,k;
    int swap = 0;

    printf ("Loading density from %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"r")))
    {
	fprintf(stderr,"File %s does not exist.\n",fname);
	return 1;
    }

    fread(&i,sizeof(int),1,f);

    if (i!=sizeof(snapshot_header)+sizeof(float)*4+sizeof(int)*4+sizeof(double)*2) swap = 1;

    if (swap) printf ("(Swapping)");fflush(0);

    fread(&(density->header),sizeof(snapshot_header),1,f);
    fread(&(density->NNodes),sizeof(long long),1,f);
    fread(&(density->delta),sizeof(float),1,f);
    fread(&(density->Nx),sizeof(int),1,f);
    fread(&(density->Ny),sizeof(int),1,f);
    fread(&(density->Nz),sizeof(int),1,f);
    fread(&(density->sigma),sizeof(float),1,f);
    fread(&(density->nsigma),sizeof(float),1,f);
    fread(&(density->Norm),sizeof(double),1,f);
    fread(&(density->B),sizeof(double),1,f);
    fread(&(density->Mass),sizeof(float),1,f);
    fread(&i,sizeof(int),1,f);

    if (swap)
    {
	Dswap4BArr(density->header.npart,6);
	Dswap4BArr(density->header.npartTotal,6);
	Dswap8BArr(density->header.mass,6);
	Dswap4B(&density->header.flag_sfr);
	Dswap4B(&density->header.flag_feedback);
	Dswap4B(&density->header.flag_cooling);
	Dswap4B(&density->header.num_files);
	Dswap8B(&density->header.time);Dswap8B(&density->header.redshift);Dswap8B(&density->header.BoxSize);
	Dswap8B(&density->header.Omega0);Dswap8B(&density->header.OmegaLambda);Dswap8B(&density->header.HubbleParam);
	Dswap8B(&(density->NNodes));
	Dswap4B(&(density->delta));
	Dswap4B(&(density->Nx));
	Dswap4B(&(density->Ny));
	Dswap4B(&(density->Nz));
	Dswap4B(&(density->sigma));
	Dswap4B(&(density->nsigma));
	Dswap8B(&(density->Norm));
	Dswap8B(&(density->B));
	Dswap4B(&(density->Mass));
    }



    //printf("(M=%1.3f)",density->Mass);
    if (!(density->grid=(float *)malloc((size_t)density->NNodes*sizeof(float))))
    {
	fprintf(stderr,"Not enough memory for density->grid while loading\n");
	exit(0);
    }

    fread(&i,sizeof(int),1,f);
    /*
    for (j=0;j<8;j++)
      fread(&(density->grid[j*(density->NNodes/8)]),sizeof(float),density->NNodes/8,f);
    if (density->NNodes%8)
      fread(&(density->grid[j*(density->NNodes/8)]),sizeof(float),density->NNodes%8,f);
    */
    fread(&(density->grid[0]),sizeof(float),density->NNodes,f);
    fread(&i,sizeof(int),1,f);
    if (swap) Dswap4BArr(&(density->grid[0]),density->NNodes);
    
    if (density->delta==0)
    {
	if (!(density->pos=(float *)malloc(3*density->NNodes*sizeof(float))))
	{
	    fprintf(stderr,"Not enough memory for density->pos while loading\n");
	    exit(0);
	}
	fread(&i,sizeof(int),1,f);
	fread(&(density->pos[0]),sizeof(float),3*density->NNodes,f);
	fread(&i,sizeof(int),1,f);
	if (swap) Dswap4BArr(&(density->pos[0]),3*density->NNodes);
    }

    if (strlen(density->header.fill)==2)
    {
	if (!(density->NSpec=(int *)malloc(density->NNodes*sizeof(int))))
	{
	    fprintf(stderr,"Not enough memory for density->pos while loading\n");
	    exit(0);
	}
	fread(&i,sizeof(int),1,f);
	fread(density->NSpec,sizeof(int),density->NNodes,f);
	fread(&i,sizeof(int),1,f);
	if (swap) Dswap4BArr(density->NSpec,density->NNodes);

	if (
	    (!(density->Spec=(int **)malloc(density->NNodes*sizeof(int*))))||
	    (!(density->SpecWeight=(float **)malloc(density->NNodes*sizeof(float*))))
	    )
	{
	    fprintf(stderr,"Not enough memory for density->pos while loading\n");
	    exit(0);
	}
	for (i=0;i<density->NNodes;i++)
	    if (
		(!(density->Spec[i]=(int *)malloc(density->NSpec[i]*sizeof(int))))||
		(!(density->SpecWeight[i]=(float *)malloc(density->NSpec[i]*sizeof(float))))
		)
	    {
		fprintf(stderr,"Not enough memory for density->pos while loading\n");
		exit(0);
	    }
	
	fread(&i,sizeof(int),1,f);
	for (i=0;i<density->NNodes;i++)
	{
	    fread(density->Spec[i],sizeof(int),density->NSpec[i],f);
	    if (swap) Dswap4BArr(density->Spec[i],density->NSpec[i]);
	}
	fread(&i,sizeof(int),1,f);

	fread(&i,sizeof(int),1,f);
	for (i=0;i<density->NNodes;i++)
	{
	    fread(density->SpecWeight[i],sizeof(float),density->NSpec[i],f);
	    if (swap) Dswap4BArr(density->SpecWeight[i],density->NSpec[i]);
	}
	fread(&i,sizeof(int),1,f);
    }
    
    
    fclose (f);
    //printf("N=%d\n",density->NNodes);
    printf ("done.\n");fflush(0);

    return 0;
}

int SaveDensityGrid(char *fname,density_grid *density)
{
    FILE *f;
    unsigned int i,j;
    long long l;
    char test[30];

    printf ("Saving density to %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"w")))
    {
	fprintf(stderr,"Could not open %s for writing.\n",fname);
	return 1;
    }

    strcpy(test,"Density grid file");

    i=30*sizeof(char);
    fwrite(&i,sizeof(int),1,f);
    fwrite(test,30*sizeof(char),1,f);
    fwrite(&i,sizeof(int),1,f);

    i=sizeof(int)*4+sizeof(float)*6;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->Nx),sizeof(int),1,f);
    fwrite(&(density->Ny),sizeof(int),1,f);
    fwrite(&(density->Nz),sizeof(int),1,f);
    fwrite(&(density->NNodes),sizeof(long long),1,f);
    printf("N nodes to write %lld\n", density->NNodes);
    fwrite(&(density->x0),sizeof(float),1,f);
    fwrite(&(density->y0),sizeof(float),1,f);
    fwrite(&(density->z0),sizeof(float),1,f);
    
    fwrite(&(density->dx),sizeof(float),1,f);
    fwrite(&(density->dy),sizeof(float),1,f);
    fwrite(&(density->dz),sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    i=sizeof(double)*density->NNodes;
    if(i==0)
	i=4;
    fprintf(stdout, "pad Save %d\n", i);
    
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->grid[0]),sizeof(FLOAT),density->NNodes,f);
    fwrite(&i,sizeof(int),1,f);
    

    i=2*sizeof(float);
    fwrite(&i,sizeof(int),1,f);
    fwrite(&density->redshift,sizeof(float),1,f);
    fwrite(&density->smoothing,sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);

    printf ("done.\n");
}

int SaveEigenvalueGrid(char *fname, density_grid *density, int n_eigen)
{
    FILE *f;
    unsigned int i,j;
    char test[30];
    char new_name[1000];
    
    sprintf(new_name,"%s_%d",fname, n_eigen);

    printf ("Saving eigenvalues to %s ... ",new_name);fflush(0);
    
    if(!(f = fopen(new_name,"w")))
    {
	fprintf(stderr,"Could not open %s for writing.\n",new_name);
	return 1;
    }

    strcpy(test,"Eigenvalue grid file");

    i=30*sizeof(char);
    fwrite(&i,sizeof(int),1,f);
    fwrite(test,30*sizeof(char),1,f);
    fwrite(&i,sizeof(int),1,f);

    i=sizeof(int)*4+sizeof(float)*6;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->Nx),sizeof(int),1,f);
    fwrite(&(density->Ny),sizeof(int),1,f);
    fwrite(&(density->Nz),sizeof(int),1,f);
    fwrite(&(density->NNodes),sizeof(long long),1,f);
    printf("N nodes to write %d\n", density->NNodes);
    fwrite(&(density->x0),sizeof(float),1,f);
    fwrite(&(density->y0),sizeof(float),1,f);
    fwrite(&(density->z0),sizeof(float),1,f);
    
    fwrite(&(density->dx),sizeof(float),1,f);
    fwrite(&(density->dy),sizeof(float),1,f);
    fwrite(&(density->dz),sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    i=sizeof(float)*density->NNodes;
    fwrite(&i,sizeof(int),1,f);
    if(n_eigen==1)
	fwrite(&(density->eigenvalue_1[0]),sizeof(float),density->NNodes,f);
    if(n_eigen==2)
	fwrite(&(density->eigenvalue_2[0]),sizeof(float),density->NNodes,f);
    if(n_eigen==3)
	fwrite(&(density->eigenvalue_3[0]),sizeof(float),density->NNodes,f);
    fwrite(&i,sizeof(int),1,f);

    i=2*sizeof(float);
    fwrite(&i,sizeof(int),1,f);
    fwrite(&density->redshift,sizeof(float),1,f);
    fwrite(&density->smoothing,sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);

    printf ("done.\n");
}

int SaveEigenvectorGrid(char *fname, density_grid *density, int n_eigen)
{
    FILE *f;
    unsigned int i,j;
    char test[30];
    char new_name[1000];
    
    sprintf(new_name,"%s_%d",fname, n_eigen);

    printf ("Saving eigenvectors to %s ... ",new_name);fflush(0);
    
    if(!(f = fopen(new_name,"w")))
    {
	fprintf(stderr,"Could not open %s for writing.\n",new_name);
	return 1;
    }

    strcpy(test,"Eigenvector grid file");

    i=30*sizeof(char);
    fwrite(&i,sizeof(int),1,f);
    fwrite(test,30*sizeof(char),1,f);
    fwrite(&i,sizeof(int),1,f);

    i=sizeof(int)*4+sizeof(float)*6;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->Nx),sizeof(int),1,f);
    fwrite(&(density->Ny),sizeof(int),1,f);
    fwrite(&(density->Nz),sizeof(int),1,f);
    fwrite(&(density->NNodes),sizeof(long long),1,f);
    printf("N nodes to write %ld\n", density->NNodes*3);
    fwrite(&(density->x0),sizeof(float),1,f);
    fwrite(&(density->y0),sizeof(float),1,f);
    fwrite(&(density->z0),sizeof(float),1,f);
    
    fwrite(&(density->dx),sizeof(float),1,f);
    fwrite(&(density->dy),sizeof(float),1,f);
    fwrite(&(density->dz),sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    i=sizeof(float)*density->NNodes*3;
    fwrite(&i,sizeof(int),1,f);
    if(n_eigen==1)
	fwrite(&(density->eigenvector_1[0]),sizeof(float),3*density->NNodes,f);
    if(n_eigen==2)
	fwrite(&(density->eigenvector_2[0]),sizeof(float),3*density->NNodes,f);
    if(n_eigen==3)
	fwrite(&(density->eigenvector_3[0]),sizeof(float),3*density->NNodes,f);
    fwrite(&i,sizeof(int),1,f);

    i=2*sizeof(float);
    fwrite(&i,sizeof(int),1,f);
    fwrite(&density->redshift,sizeof(float),1,f);
    fwrite(&density->smoothing,sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);

    printf ("done.\n");
    return(0);
}


int SaveTraceGrid(char *fname,density_grid *density)
{
    FILE *f;
    unsigned int i,j;
    char test[30];

    printf ("Saving eigenvalues to %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"w")))
    {
	fprintf(stderr,"Could not open %s for writing.\n",fname);
	return 1;
    }

    strcpy(test,"Eigenvalue grid file");

    i=30*sizeof(char);
    fwrite(&i,sizeof(int),1,f);
    fwrite(test,30*sizeof(char),1,f);
    fwrite(&i,sizeof(int),1,f);

    i=sizeof(int)*4+sizeof(float)*6;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->Nx),sizeof(int),1,f);
    fwrite(&(density->Ny),sizeof(int),1,f);
    fwrite(&(density->Nz),sizeof(int),1,f);
    fwrite(&(density->NNodes),sizeof(long long),1,f);
    printf("N nodes to write %d\n", density->NNodes);
    fwrite(&(density->x0),sizeof(float),1,f);
    fwrite(&(density->y0),sizeof(float),1,f);
    fwrite(&(density->z0),sizeof(float),1,f);
    
    fwrite(&(density->dx),sizeof(float),1,f);
    fwrite(&(density->dy),sizeof(float),1,f);
    fwrite(&(density->dz),sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    i=sizeof(double)*density->NNodes;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->trace[0]),sizeof(FLOAT),density->NNodes,f);
    fwrite(&i,sizeof(int),1,f);

    i=2*sizeof(float);
    fwrite(&i,sizeof(int),1,f);
    fwrite(&density->redshift,sizeof(float),1,f);
    fwrite(&density->smoothing,sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    printf ("done.\n");
}

int SaveEnvGrid(char *fname,density_grid *density)
{
    FILE *f;
    unsigned int i,j;
    char test[30];

    printf ("Saving environment to %s ... ",fname);fflush(0);
    
    if(!(f = fopen(fname,"w")))
    {
	fprintf(stderr,"Could not open %s for writing.\n",fname);
	return 1;
    }

    strcpy(test,"Environment grid file");

    i=30*sizeof(char);
    fwrite(&i,sizeof(int),1,f);
    fwrite(test,30*sizeof(char),1,f);
    fwrite(&i,sizeof(int),1,f);

    i=sizeof(int)*4+sizeof(float)*6;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->Nx),sizeof(int),1,f);
    fwrite(&(density->Ny),sizeof(int),1,f);
    fwrite(&(density->Nz),sizeof(int),1,f);
    fwrite(&(density->NNodes),sizeof(int),1,f);
    printf("N nodes to write %d\n", density->NNodes);
    fwrite(&(density->x0),sizeof(float),1,f);
    fwrite(&(density->y0),sizeof(float),1,f);
    fwrite(&(density->z0),sizeof(float),1,f);
    
    fwrite(&(density->dx),sizeof(float),1,f);
    fwrite(&(density->dy),sizeof(float),1,f);
    fwrite(&(density->dz),sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    i=sizeof(double)*density->NNodes;
    fwrite(&i,sizeof(int),1,f);
    fwrite(&(density->environment[0]),sizeof(int),density->NNodes,f);
    fwrite(&i,sizeof(int),1,f);

    i=2*sizeof(float);
    fwrite(&i,sizeof(int),1,f);
    fwrite(&density->redshift,sizeof(float),1,f);
    fwrite(&density->smoothing,sizeof(float),1,f);
    fwrite(&i,sizeof(int),1,f);
    
    printf ("done.\n");
    return 0;
}

int WriteGadget(char *fname,snapshot_data *snap,int flags)
{
 int i,j,k;
  FILE *f;

  if ((f=fopen(fname,"w"))==NULL)
    {
      fprintf (stderr,"Could not open file %s for writing\n",fname);
      return -1;
    }

  printf ("Saving gadget file %s (N=%d)",fname,snap->N);fflush(0);

  i=256;
  fwrite(&i,sizeof(int),1,f);
  fwrite(&snap->header,i,1,f);
  fwrite(&i,sizeof(int),1,f);

  if ((flags&FLAG_POS)&&(snap->Pos!=NULL))
    {
      i=3*sizeof(float)*snap->N;
      printf ("(POS:%.1f)",(float)i/1024/1024);fflush(0);
      fwrite(&i,sizeof(int),1,f);
      fwrite(snap->Pos,sizeof(float),3*snap->N,f);
      fwrite(&i,sizeof(int),1,f);
    }
  if ((flags&FLAG_VEL)&&(snap->Vel!=NULL))
    {
      i=3*sizeof(float)*snap->N;
      printf ("(VEL:%.1f)",(float)i/1024/1024);fflush(0);
      fwrite(&i,sizeof(int),1,f);
      fwrite(snap->Vel,sizeof(float),3*snap->N,f);
      fwrite(&i,sizeof(int),1,f);
    }
  if ((flags&FLAG_ID)&&(snap->Id!=NULL))
    {
      i=sizeof(int)*snap->N;
      printf ("(ID:%.1f)",(float)i/1024/1024);fflush(0);
      fwrite(&i,sizeof(int),1,f);
      fwrite(snap->Id,sizeof(int),snap->N,f);
      fwrite(&i,sizeof(int),1,f);
    }

  fclose(f);
  printf (" done.\n");
  return 0;
}








