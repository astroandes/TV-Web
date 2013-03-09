/***********************************************************

  Date        : Reads an ART snapshot stored in one file.
  Author      : Jaime Forero
  Description : .h file for READ_ART.c

  Notes       : It is based on the  code READ_ART.f

***********************************************************/

#define MIN(x,y) (x>y?y:x)
#define MAX(x,y) (x>y?x:y)
#define NROW 4096 //maximum number of particles in 1D//Bolshoi
//#define NROW 1024 //maximum number of particles in 1D//BigBolshoi
#define NGRID 256 //zero-level mesh in 1D
#define NP    (1024*1024*1024) //max number of particles, Np must be less than NROW**3
#define NPAGE (NROW*NROW)//number of particles in a record
#define NRECL (NPAGE*6)//number of words in a record

typedef struct part{
  float x, y, z, vx, vy, vz, iweigth;
} part;

int read_art_snap(char *, char *, snapshot_data *);
int read_art_snap_multifile(char * namein, char * datain, snapshot_data *P, int n_files, int my_file);
void replace_string(char *original, char *search, char *replace);
