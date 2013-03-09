#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include "io.h"

#define TMPFILENAME "__155266488855212__.tmp"

int countcol(const char *str)
{
  int was_space=1;
  int n=0;
  int i;

  for (i=0;i<strlen(str);i++)
    {
      
      if ((str[i]!=' ')&&(was_space)) 
	{n++;was_space=0;}
      
      if (str[i]==' ') was_space=1;

    }

  return n;

}

int main (int argc,char **argv)
{
  int i,j,k,nf;
  snapshot_data* data=NULL;
  FILE *f;
  int nx,ny,nz;

  int tab1[100][100];
  int tab2[100][100];
  int tab3[100][100];

  for (i=0;i<100;i++)
    for (j=0;j<100;j++)
      {
	tab1[i][j]=tab2[i][j]=tab3[i][j]=0;
      }

  data=calloc(1,sizeof(snapshot_data));

  ReadGadget(argv[1],data,FLAG_POS);
  /*
  f=fopen("ftest.dat","w");

  for (i=0;i<data->N;i++)
    {
      if ((data->Pos[3*i+2]<0.2)&&(data->Pos[3*i+2]>0.15))
	{
	  fprintf (f,"%e %e\n",data->Pos[3*i+0],data->Pos[3*i+1]);
	}
    }

    fclose(f);
  */
  

  for (i=0;i<3*data->N;i+=3)
    {
      //nx=(int)(data->Pos[i]/data->header.BoxSize * 100);
      //ny=(int)(data->Pos[i+1]/data->header.BoxSize*100.);
      //nz=(int)(data->Pos[i+2]/data->header.BoxSize*100.);
      nx=(int)(data->Pos[i] * 100);
      ny=(int)(data->Pos[i+1]*100.);
      nz=(int)(data->Pos[i+2]*100.);
      if (nx>=100) nx-=100;
      if (ny>=100) ny-=100;
      if (nz>=100) nz-=100;

      if (nx<0) nx+=100;
      if (ny<0) ny+=100;
      if (nz<0) nz+=100;
      
      tab1[nx][ny]++;
      tab2[ny][nz]++;
      tab3[nx][nz]++;
    }     

  printf("projecting %d particles.\n",data->N);

  f=fopen("fxy.dat","w");
  for (i=0;i<100;i++)
    for (j=0;j<100;j++)
      {
	fprintf (f,"%d\n",tab1[i][j]);
      }
  fclose(f);
  
  f=fopen("fyz.dat","w");
  for (i=0;i<100;i++)
    for (j=0;j<100;j++)
      {
	fprintf (f,"%d\n",tab2[i][j]);
      }
  fclose(f);
  
  f=fopen("fxz.dat","w");
  for (i=0;i<100;i++)
    for (j=0;j<100;j++)
      {
	fprintf (f,"%d\n",tab3[i][j]);
      }
  fclose(f);

  return 0;
}
