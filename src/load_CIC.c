int LoadDensityGrid(char *fname,density_grid *density)
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
    fread_sw(&(density->NNodes),sizeof(int),1,f,swap);

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


    return 0;
}
