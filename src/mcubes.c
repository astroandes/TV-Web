#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "field.h"
#include "mcubes.h"
#include "mctable.h"

#define MC_NONE  0
#define MC_LEFT  1
#define MC_RIGHT 2
#define MC_UP    3
#define MC_DOWN  4
#define MC_REAR  5
#define MC_FRONT 6

#define MC_REAL FLOAT

#define PI 3.14159265358979323846

//float MC_Level;
//MC_Object *MC_Obj;
//density_field *MC_Density;

char MC_FirstCall=1;
//char *MC_Edges[256];
//char MC_NEdges[256];
//unsigned char *MC_VertexPos[256][12];//gives the position of vertex j in config i, 0 if not used
//unsigned char MC_NVertexPos[256][12];
unsigned char MC_NFaces[256];

int free_MC_Object(MC_Object **Obj)
{
    MC_Object *o=*Obj;
    int i;

    free(o->Vertex);free(o->Face);
    free(o->Normal);free(o->Color);
    free(o->alpha);free(o->BoxIndex);
    free(o->grad);free(o->gradNorm);
    free(o->NNeighbours);
    for (i=0;i<28;i++)
	free(o->VertexNeighbour[i]);

    free(o);
	  	    
    *Obj=NULL;
}

int InitMCubes()
{
    int i,j,k;

    MC_FirstCall=0;

    for (i=0;i<256;i++)
	for (j=0;j<16;j++)
	{
	    if (triTable[i][j]==-1)
	    {
		MC_NFaces[i]=j/3;
		j=16;
	    }
	}
    return 0;
}


inline int GetCubicInter(MC_REAL level,MC_REAL f0,MC_REAL f1,MC_REAL fp0,MC_REAL fp1,MC_REAL *alpha)
{
    MC_REAL a=fp1-2*f1+fp0+2*f0;
    MC_REAL b=-fp1+3*f1-2*fp0-3*f0;
    MC_REAL c=fp0;
    MC_REAL d=f0-level;
    
    MC_REAL p,q,r,phi,D;

    int nsol=0;

    p=(3*a*c-b*b)/(9*a*a);
    q=((2*b*b*b)/(27*a*a*a)-(b*c)/(3*a*a)+d/a)/2.;
    

    r=((q<0)?-1:1)*sqrt(fabs(p));

    //printf ("%lf %lf %lf\n",p,q,D);
    if (p<0)
    {
	D=p*p*p+q*q;
	if (D>0)
	{
	    
	    //printf ("AAA*******");
	    phi = acosh(q/(r*r*r))/3.;
	    p = -2*r*cosh(phi)- b/(3*a);
	    if ((p>0)&&(p<1)) alpha[nsol++]=p;
	}
	else
	{
	    phi = acos(q/(r*r*r))/3.;
	    //printf ("phi = %lf (%lf %lf -> %lf)\n",phi,q,r,q/(r*r*r));
	    p = -2*r*cos(phi) - b/(3*a);
	    if ((p>0)&&(p<1)) alpha[nsol++]=p;

	    p = 2*r*cos(PI/3-phi)-b/(3*a);
	    if ((p>0)&&(p<1)) alpha[nsol++]=p;

	    p = 2*r*cos(PI/3+phi)-b/(3*a);
	    if ((p>0)&&(p<1)) alpha[nsol++]=p;
	}
    }
    else
    {
	
	//printf ("BBB*******");
	phi = asinh(q/(r*r*r))/3.;
	p = -2*r*sinh(phi)- b/(3*a);
	if ((p>0)&&(p<1)) alpha[nsol++]=p;
    }
    //printf ("%d sol.\n",nsol);
    return nsol;
}

//returns the value of interpolated at point 0<alpha<1
//using cubic interpolation
//alpha=0 returns f0 and alpha=1 returns f1
/*
inline MC_REAL CubicInter(MC_REAL f0,MC_REAL f1,MC_REAL fp0,MC_REAL fp1,MC_REAL alpha)
{
    MC_REAL a=fp1-2*f1+fp0+2*f0;
    MC_REAL a=-fp1+3*f1-2*fp0-3*f0;
    MC_REAL a=fp0;
    MC_REAL a=f0;

    return d+alpha*(c+alpha*(b+a*alpha));     
}
*/

//Object should be initialized with calloc (pointer must be NULL or allocated)
int GetIsoSurface(density_grid *density,float level,MC_Object *Object,int p_flags)
{
    unsigned char *MC_Box; 

    MC_REAL x,y,z;
    MC_REAL u,v,w;
    MC_REAL a,b,c;
    int i,j,k,l;
    int n,n3,np,nf;
    int p1,p2,p3;
    int val;
    int dx,dy,dz;
    int dx3,dy3,dz3;
    MC_REAL deltax,deltay,deltaz;
    MC_REAL alpha;
    MC_REAL alpha3f[3];
    float *Point;
    int *Face;
    float *Normal;
    MC_REAL *DGrad;
    float *OGrad;
    int edge;
    int sp=0;
    int *NPoints;//Number of points up to node i
    MC_REAL minpos;
    MC_REAL norm;
    int tmp;
    int flags;
    double dv;

    //printf ("flags: %d\n",flags);
    //printf("Called\n");
    //printf ("%d %d %d\n",density->Nx,density->Ny,density->Nz);
    if (MC_FirstCall) InitMCubes();
    if (p_flags&MCFLAG_CORRECT) 
	flags = p_flags|MCFLAG_GETGRAD;
    else
	flags = p_flags;

    MC_Box = (unsigned char *) calloc(density->NNodes,sizeof(char));
    NPoints = (int *) malloc(density->NNodes*sizeof(int));
    //printf("%d\n",NPoints);

    dx = 1;
    dy = density->Nx;
    dz = density->Ny*density->Nx;

    dx3=3*dx;dy3=3*dy;dz3=3*dz;

/*
    //first give every cell a name (0..255)
    for (k=1,n=dz;k<density->Nz-1;k++,n+=dy)
	for (j=1,n+=dy;j<density->Ny-1;j++,n++)
	    for (i=1,n++;i<density->Nx-1;i++,n++)
	    {
		//printf ("%d_",density->NNodes-n);
		if (density->grid[n]>level)
		{
		    //Marks the corresponding corner of every concerned box
		    MC_Box[n]|=(1<<0);
		    MC_Box[n-dx]|=(1<<1);
		    MC_Box[n-dx-dz]|=(1<<2);
		    MC_Box[n-dz]|=(1<<3);
		    MC_Box[n-dy]|=(1<<4);
		    MC_Box[n-dx-dy]|=(1<<5);
		    MC_Box[n-dx-dz-dy]|=(1<<6);
		    MC_Box[n-dz-dy]|=(1<<7);
		}
	    }
*/
    //first give every cell a name (0..255)

    for (k=MC_DELTA,n=MC_DELTA*dz;k<=density->Nz-MC_DELTA;k++,n+=(MC_DELTA-1)*dy)
	for (j=MC_DELTA,n+=MC_DELTA*dy;j<=density->Ny-MC_DELTA;j++,n+=(MC_DELTA-1))
	    for (i=MC_DELTA,n+=MC_DELTA;i<=density->Nx-MC_DELTA;i++,n++)
	    {
		if ((i==density->Nx-MC_DELTA)&&(density->grid[n-dx]>level)) MC_Box[n-dx]|=(1<<1);
		else if ((j==density->Ny-MC_DELTA)&&(density->grid[n-dy]>level)) MC_Box[n-dy]|=(1<<4);
		else if ((k==density->Nz-MC_DELTA)&&(density->grid[n-dz]>level)) MC_Box[n-dz]|=(1<<3);
		
		//printf ("%d_",density->NNodes-n);
		else if (density->grid[n]>level)
		{
		    //Marks the corresponding corner of every concerned box
		    MC_Box[n]|=(1<<0);
		    MC_Box[n-dx]|=(1<<1);
		    MC_Box[n-dx-dz]|=(1<<2);
		    MC_Box[n-dz]|=(1<<3);
		    MC_Box[n-dy]|=(1<<4);
		    MC_Box[n-dx-dy]|=(1<<5);
		    MC_Box[n-dx-dz-dy]|=(1<<6);
		    MC_Box[n-dz-dy]|=(1<<7);
		}
	    }

    Object->VolumeUp=0;
    Object->VolumeDown=0;

    if (flags&MCFLAG_GETVOLUME)
    {
	dv = density->dx*density->dy*density->dz/8.;

	for (k=MC_DELTA,n=MC_DELTA*dz;k<density->Nz-1-MC_DELTA;k++,n+=(1+MC_DELTA)*dy)
	    for (j=MC_DELTA,n+=MC_DELTA*dy;j<density->Ny-1-MC_DELTA;j++,n+=1+MC_DELTA)
		for (i=MC_DELTA,n+=MC_DELTA;i<density->Nx-1-MC_DELTA;i++,n++)
		{
		    val = MC_Box[n];
		
		    for (l=0;l<8;l++) 
			if (val&(1<<l)) 
			    Object->VolumeUp+=dv;
			else
			    Object->VolumeDown+=dv;
		}
	//printf ("Volume: %g/%g\n",Object->VolumeUp,Object->VolumeDown);
    }

    //Compute the number of faces
    Object->NFaces=0;
    for (k=MC_DELTA,n=MC_DELTA*dz;k<density->Nz-1-MC_DELTA;k++,n+=(1+MC_DELTA)*dy)
	for (j=MC_DELTA,n+=MC_DELTA*dy;j<density->Ny-1-MC_DELTA;j++,n+=1+MC_DELTA)
	    for (i=MC_DELTA,n+=MC_DELTA;i<density->Nx-1-MC_DELTA;i++,n++)
	    {
		val = MC_Box[n];
		if (MC_NFaces[val])
		{
		    Object->NFaces+=MC_NFaces[val];
		}
	    }

   
    //Estimate Object->NVertices --> NPoints~<NFaces 
    //every point is shared with at least 3 faces (isosurfaces are closed) and 1 face has 3 points ...
    //printf ("%d ->\n",Object->NFaces);
    //printf ("%d ->\n",Object->Vertex);



    if (flags&MCFLAG_CORRECT) 
    {
	Object->Invalid=(char *) realloc(Object->Invalid,density->NNodes*sizeof(char));
	Object->Color=(float *) calloc(3*10000,sizeof(float));
	memset (Object->Invalid,0,density->NNodes*sizeof(char));
    }
    
    Object->Vertex=(float *) realloc(Object->Vertex,(1000+3*Object->NFaces)*sizeof(float));
    Object->Face=(int *) realloc(Object->Face,3*Object->NFaces*sizeof(int));
    //printf ("%d ->\n",Object->Vertex);
    //printf("NFaces=%d\n",Object->NFaces);fflush(0);
    Object->NVertices=0;
    np=0;i=0;
    deltax=density->dx;
    deltay=density->dy;
    deltaz=density->dz;
    Point=Object->Vertex;
    
    //Now get the points
    //if we compute the gradiant of the field too... 
    if (flags&MCFLAG_GETGRAD)
    {
	int ttt=-1;
	int correction = flags&MCFLAG_CORRECT;
	double ooo1=0,ooo2=0,ooo3=0;
	int nu=0,nv=0,nw=0;
	Object->grad = (float *) realloc(Object->grad,3000+6*Object->NFaces*sizeof(float));
	//with gradiant of the field
	if (!density->HasGradiant) ComputeDensityGradiant(density,USE_FD);
	
	DGrad = density->grad;
	OGrad = Object->grad;
	
	for (z=density->z0+MC_DELTA*deltaz,k=MC_DELTA,n=MC_DELTA*dz;k<density->Nz-MC_DELTA;z+=deltaz,k++,n+=MC_DELTA*dy)
	    for (y=density->y0+MC_DELTA*deltay,j=MC_DELTA,n+=MC_DELTA*dy;j<density->Ny-MC_DELTA;y+=deltay,j++,n+=MC_DELTA)
		for (x=density->x0+MC_DELTA*deltax,i=MC_DELTA,n+=MC_DELTA;i<density->Nx-MC_DELTA;x+=deltax,i++,n++)
		{
		    
		    //printf ("%d_",n);
		    val = MC_Box[n];
		    edge = edgeTable[val];
		    NPoints[n]=np;//printf ("%d_",NPoints[n]);
		    n3=3*n;

		    //We only compute edges 8,3,0 on every cube (with interpolated gradiant)
		    if ((edge)||(correction))
		    {
		    if (edge&(1<<8))
		    {
			tmp = GetCubicInter(level,density->grid[n],density->grid[n+dy],DGrad[n3+1],DGrad[(n3+dy3)+1],alpha3f);
			
			if (tmp!=1)
			{
			    alpha=(level-density->grid[n])/(density->grid[n+dy]-density->grid[n]);
			    if (correction) {
				Object->Invalid[n]=1;ttt++;
				Object->Color[3*ttt] = x+deltax/2;
				Object->Color[3*ttt+1] = y+deltay/2;
				Object->Color[3*ttt+2] = z+deltaz/2;
			    }
			    //printf ("Bugged box(y,%d) ... alpha = %f\n",tmp,alpha);
			}
			else 
			    alpha=alpha3f[0];
	
			OGrad[np]=DGrad[n3]+alpha*(DGrad[n3+dy3]-DGrad[n3]);
			OGrad[np+1]=DGrad[n3+1]+alpha*(DGrad[(n3+dy3)+1]-DGrad[n3+1]);
			OGrad[np+2]=DGrad[n3+2]+alpha*(DGrad[(n3+dy3)+2]-DGrad[n3+2]);

			Point[np++]=x;
			Point[np++]=y+alpha*deltay;
			Point[np++]=z;

		    }
		    
		    else if ((correction)&&(j!=density->Ny-MC_DELTA-1))
		    {
			tmp = GetCubicInter(level,density->grid[n],density->grid[n+dy],DGrad[n3+1],DGrad[(n3+dy3)+1],alpha3f);
						
			if (tmp!= 0)
			{
			    Object->Invalid[n]=1;ttt++;
			    
			    Object->Color[3*ttt] = x+deltax/2;
			    Object->Color[3*ttt+1] = y+deltay/2;
			    Object->Color[3*ttt+2] = z+deltaz/2;
			}
		    }
	

		    if (edge&(1<<3))
		    {
			tmp = GetCubicInter(level,density->grid[n],density->grid[n+dz],DGrad[n3+2],DGrad[(n3+dz3)+2],alpha3f);
			
			if (tmp!=1)
			{
			    alpha=(level-density->grid[n])/(density->grid[n+dz]-density->grid[n]);
			    if (correction) {
				Object->Invalid[n]=1;ttt++;
				Object->Color[3*ttt] = x+deltax/2;
				Object->Color[3*ttt+1] = y+deltay/2;
				Object->Color[3*ttt+2] = z+deltaz/2;
			    }
			    //printf ("Bugged box(z,%d) ...alpha = %f, %f:(%f<>%f)\n",tmp,alpha,level,density->grid[n],density->grid[n+dz]);
			}
			else alpha=alpha3f[0];
			
			//printf ("alpha=%f\n",alpha);

			OGrad[np]=DGrad[n3]+alpha*(DGrad[(n3+dz3)]-DGrad[n3]);
			OGrad[np+1]=DGrad[n3+1]+alpha*(DGrad[(n3+dz3)+1]-DGrad[n3+1]);
			OGrad[np+2]=DGrad[n3+2]+alpha*(DGrad[(n3+dz3)+2]-DGrad[n3+2]);

			Point[np++]=x;
			Point[np++]=y;
			Point[np++]=z+alpha*deltaz;

		    }
		   
	
		    else if ((correction)&&(k!=density->Nz-MC_DELTA-1))
		    {
			tmp = GetCubicInter(level,density->grid[n],density->grid[n+dz],DGrad[n3+2],DGrad[(n3+dz3)+2],alpha3f);
			if (tmp!= 0)
			{
			    Object->Invalid[n]=1;ttt++;
			    Object->Color[3*ttt] = x+deltax/2;
			    Object->Color[3*ttt+1] = y+deltay/2;
			    Object->Color[3*ttt+2] = z+deltaz/2;
			}
		    }

		    if (edge&(1<<0))
		    {
			
			tmp = GetCubicInter(level,density->grid[n],density->grid[n+dx],DGrad[n3],DGrad[(n3+dx3)],alpha3f);
		
			if (tmp!=1)
			{
			    alpha=(level-density->grid[n])/(density->grid[n+dx]-density->grid[n]);
			    if (correction) {
				Object->Invalid[n]=1;ttt++;
				Object->Color[3*ttt] = x+deltax/2;
				Object->Color[3*ttt+1] = y+deltay/2;
				Object->Color[3*ttt+2] = z+deltaz/2;
			    }
			    //printf ("Bugged box(x,%d) ...alpha = %f (%f)\n",tmp,alpha,alpha3f[0]);
			}
			else 
			    alpha=alpha3f[0];
			

			OGrad[np]=DGrad[n3]+alpha*(DGrad[(n3+dx3)]-DGrad[n3]);
			OGrad[np+1]=DGrad[n3+1]+alpha*(DGrad[(n3+dx3)+1]-DGrad[n3+1]);
			OGrad[np+2]=DGrad[n3+2]+alpha*(DGrad[(n3+dx3)+2]-DGrad[n3+2]);

			Point[np++]=x+alpha*deltax;
			Point[np++]=y;
			Point[np++]=z;
		    }
		
		    
		    else if ((correction)&&(i!=density->Nx-MC_DELTA-1))
		    {
			tmp = GetCubicInter(level,density->grid[n],density->grid[n+dx],DGrad[n3],DGrad[(n3+dx3)],alpha3f);
			if (tmp!= 0)
			{
			   Object->Invalid[n]=1; ttt++;
			   Object->Color[3*ttt] = x+deltax/2;
			   Object->Color[3*ttt+1] = y+deltay/2;
			   Object->Color[3*ttt+2] = z+deltaz/2;
			}			
		    }
		    
		    }
		}
	Object->grad = (float *) realloc(Object->grad,np*sizeof(float));
	if (correction) printf ("%d invalid edges ...\n",ttt);
    }
    else
    {
	//without gradiant of the field
	for (z=density->z0+MC_DELTA*deltaz,k=MC_DELTA,n=MC_DELTA*dz;k<density->Nz-MC_DELTA;z+=deltaz,k++,n+=MC_DELTA*dy)
	    for (y=density->y0+MC_DELTA*deltay,j=MC_DELTA,n+=MC_DELTA*dy;j<density->Ny-MC_DELTA;y+=deltay,j++,n+=MC_DELTA)
		for (x=density->x0+MC_DELTA*deltax,i=MC_DELTA,n+=MC_DELTA;i<density->Nx-MC_DELTA;x+=deltax,i++,n++)
		{
		    val = MC_Box[n];
		    edge = edgeTable[val];
		    NPoints[n]=np;
		    
		    //We only compute edges 8,3,0 on every cube
		    if (edge&(1<<8))
		    {
			alpha=(level-density->grid[n])/(density->grid[n+dy]-density->grid[n]);
			
			Point[np++]=x;
			Point[np++]=(double)y+alpha*deltay;
			Point[np++]=z;
		    }
		    
		    if (edge&(1<<3))
		    {
			alpha=(level-density->grid[n])/(density->grid[n+dz]-density->grid[n]);
			
			Point[np++]=x;
			Point[np++]=y;
			Point[np++]=(double)z+alpha*deltaz;
		    }
		    
		    if (edge&(1<<0))
		    {
			alpha=(level-density->grid[n])/(density->grid[n+dx]-density->grid[n]);
			
			Point[np++]=(double)x+alpha*deltax;
			Point[np++]=y;
			Point[np++]=z;
		    }
		}
    }
    //printf ("%d\n",np);
    //printf ("BBBBBBBBBBBBBB\n");
    
    Object->Vertex=realloc(Object->Vertex,np*sizeof(float));
    
    //printf ("BBBBBBBBBBBBBB\n");
    if (flags&MCFLAG_GETVERTEXTREE)
    {
	Object->NNeighbours = (int *) realloc(Object->NNeighbours,np/3*sizeof(int));
	memset(Object->NNeighbours,0,np/3*sizeof(int));
	
	for (i=0;i<28;i++)
	    Object->VertexNeighbour[i] = (int *) realloc(Object->VertexNeighbour[i],np/3*sizeof(int));
    }

    Object->NVertices=np/3;
    //printf("NPoints=%d\n",Object->NVertices);
    Face=Object->Face;
    //Still need to get the faces
    nf=0;

    if (flags&MCFLAG_MARKFACES)
	Object->BoxIndex = (int *) realloc(Object->BoxIndex,Object->NFaces*sizeof(int));

    //Now link the faces to their vertices    
    for (k=MC_DELTA,n=MC_DELTA*dz;k<density->Nz-1-MC_DELTA;k++,n+=(1+MC_DELTA)*dy)
	for (j=MC_DELTA,n+=MC_DELTA*dy;j<density->Ny-1-MC_DELTA;j++,n+=1+MC_DELTA)
	    for (i=MC_DELTA,n+=MC_DELTA;i<density->Nx-1-MC_DELTA;i++,n++)
	    {
		val = MC_Box[n];
		//if (!((!(flags&MCFLAG_CORRECT))||(!Object->Invalid[n])))
		//    printf ("deleted %d\n",n);

		if ((MC_NFaces[val])&&  ((!(flags&MCFLAG_CORRECT))||(!Object->Invalid[n])) )
		    //if (MC_NFaces[val])
		{
		    if (flags&MCFLAG_MARKFACES)
			for (l=nf/3;l<nf/3+MC_NFaces[val];l++)
			    Object->BoxIndex[l]=n;
		    
		    for (l=0;l<MC_NFaces[val]*3;l++)
		    {
			//here comes the boring part ...
			switch (triTable[val][l])
			{
			    //depending on the case, associate a face with three points ...
			    case 0:
				edge = edgeTable[val];
				if (edge&(1<<8))
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n]+6;
				    else
					Face[nf++]=NPoints[n]+3;
				else
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n]+3;
				    else
					Face[nf++]=NPoints[n];
				break;
			    case 1:
				edge = edgeTable[MC_Box[n+dx]];
				if (edge&(1<<8))
				    Face[nf++]=NPoints[n+dx]+3;
				else
				    Face[nf++]=NPoints[n+dx];
				
				break;
			    case 2:
				edge = edgeTable[MC_Box[n+dz]];
				if (edge&(1<<8))
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n+dz]+6;
				    else
					Face[nf++]=NPoints[n+dz]+3;   //Ugly, got to clean this ...

				else
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n+dz]+3;
				    else
					Face[nf++]=NPoints[n+dz];
				
				break;
			    case 3:
				edge = edgeTable[val];
				if (edge&(1<<8))
				    Face[nf++]=NPoints[n]+3;
				else
				    Face[nf++]=NPoints[n];
				
				
				break;
			    case 4:
				edge = edgeTable[MC_Box[n+dy]];
				if (edge&(1<<8))
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n+dy]+6;
				    else
					Face[nf++]=NPoints[n+dy]+3;
				else
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n+dy]+3;
				    else
					Face[nf++]=NPoints[n+dy];
				
				break;
			    case 5:
				edge = edgeTable[MC_Box[n+dx+dy]];
				if (edge&(1<<8))
				    Face[nf++]=NPoints[n+dx+dy]+3;
				else
				    Face[nf++]=NPoints[n+dx+dy];
				
				break;
			    case 6:
				edge = edgeTable[MC_Box[n+dy+dz]];
				if (edge&(1<<8))
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n+dy+dz]+6;
				    else
					Face[nf++]=NPoints[n+dy+dz]+3;
				else
				    if (edge&(1<<3))
					Face[nf++]=NPoints[n+dy+dz]+3;
				    else
					Face[nf++]=NPoints[n+dy+dz];
				
				break;
			    case 7:
				edge = edgeTable[MC_Box[n+dy]];
				if (edge&(1<<8))
				    Face[nf++]=NPoints[n+dy]+3;
				else
				    Face[nf++]=NPoints[n+dy];
				
				break;
			    case 8:
				Face[nf++]=NPoints[n];
				break;
			    case 9:
				Face[nf++]=NPoints[n+dx];
				break;
			    case 10:
				Face[nf++]=NPoints[n+dx+dz];
				break;
			    case 11:
				Face[nf++]=NPoints[n+dz];
				break;
			}
		    }
		    //Eventually compute the neighbours
		    if ((flags&MCFLAG_GETVERTEXTREE)&&(MC_NFaces[val]))
		    {
			//check all the faces we have just computed
			for (l=nf-MC_NFaces[val]*3;l<nf;l+=3)
			{
			    p1=Face[l]/3;
			    p2=Face[l+1]/3;
			    p3=Face[l+2]/3;
			    //printf ("(%d,%d,%d)",p1,p2,p3);
			    /*
			    Object->NNeighbours[p1]+=2;
			    Object->NNeighbours[p2]+=2;
			    Object->NNeighbours[p3]+=2;
			    if (Object->NNeighbours[p1]>=24)
				printf ("%d->%d (%f,%f,%f)\n",Object->NNeighbours[p1],p1
					,Object->Vertex[3*p1],Object->Vertex[3*p1+1],Object->Vertex[3*p1+2]);
			    if (Object->NNeighbours[p2]>=24)
				printf ("%d->%d (%f,%f,%f)\n",Object->NNeighbours[p2],p1
					,Object->Vertex[3*p1],Object->Vertex[3*p1+1],Object->Vertex[3*p1+2]);
			    if (Object->NNeighbours[p3]>=24)
				printf ("%d->%d (%f,%f,%f)\n",Object->NNeighbours[p3],p1
					,Object->Vertex[3*p1],Object->Vertex[3*p1+1],Object->Vertex[3*p1+2]);
			    
			    if (Object->NNeighbours[p1]>=26) printf("OOPS%d\n",Object->NNeighbours[p1]);fflush(0);
			    if (Object->NNeighbours[p2]>=26) printf("OOPS%d\n",Object->NNeighbours[p2]);fflush(0);
			    if (Object->NNeighbours[p3]>=26) printf("OOPS%d\n",Object->NNeighbours[p3]);fflush(0);
			    */

			    Object->VertexNeighbour[Object->NNeighbours[p1]++][p1]=p2;
			    Object->VertexNeighbour[Object->NNeighbours[p1]++][p1]=p3;
			    
			    Object->VertexNeighbour[Object->NNeighbours[p2]++][p2]=p3;
			    Object->VertexNeighbour[Object->NNeighbours[p2]++][p2]=p1;
			      
			    Object->VertexNeighbour[Object->NNeighbours[p3]++][p3]=p1;
			    Object->VertexNeighbour[Object->NNeighbours[p3]++][p3]=p2;
			    
			}
		    }
		}
	    }
/*
    for (i=0;i<Object->NVertices;i++) if (Object->NNeighbours[i]==0) 
	printf ("Point %d (%f %f %f)is alone !!!\n",i,Object->Vertex[3*i],Object->Vertex[3*i+1],Object->Vertex[3*i+2]);
*/
    Object->NFaces=nf/3;
    free (MC_Box);
    //printf ("p:%d,f:%d\n",Object->NVertices,Object->NFaces);
    /*
    for (l=0;l<Object->NVertices;l++)
	if (Object->NNeighbours[l]>=26)
	{
	    Object->NNeighbours[l]=0;
	    for (i=0;i<Object->NFaces;i++)
	    {
		if (Object->Face[3*i]/3==l) Object->NNeighbours[l]+=2;
		if (Object->Face[3*i+1]/3==l) Object->NNeighbours[l]+=2;
		if (Object->Face[3*i+2]/3==l) Object->NNeighbours[l]+=2;
	    }
	    printf ("%d\n",Object->NNeighbours[l]);fflush(0);
	}
    */
        
    //Finally compute the normals
    if (flags&MCFLAG_GETNORMALS)
    {
	Object->Normal=(float *) realloc(Object->Normal,3*Object->NVertices*sizeof(float));
	NPoints = (int *) realloc (NPoints,Object->NVertices*sizeof(int));

	Normal=Object->Normal;
	Point=Object->Vertex;
	Face=Object->Face;
	memset(Normal,0,3*Object->NVertices*sizeof(float));
	memset(NPoints,0,Object->NVertices*sizeof(int));
	for (i=0;i<Object->NFaces*3;i+=3)
	{
	    x=Point[Face[i]]-Point[Face[i+1]];
	    y=Point[Face[i]+1]-Point[Face[i+1]+1];
	    z=Point[Face[i]+2]-Point[Face[i+1]+2];
	    //norm=sqrt(x*x+y*y+z*z);
	    //x/=norm;y/=norm;z/=norm;
	    u=Point[Face[i]]-Point[Face[i+2]];
	    v=Point[Face[i]+1]-Point[Face[i+2]+1];
	    w=Point[Face[i]+2]-Point[Face[i+2]+2];
	    //norm=sqrt(u*u+v*v+w*w);
	    //u/=norm;v/=norm;w/=norm;
	    a=y*w-z*v;
	    b=u*z-x*w;
	    c=x*v-y*u;
	    norm=sqrt(a*a+b*b+c*c);
	    if (norm!=0)
	    {
		a/=norm;b/=norm;c/=norm;
		
		Normal[Face[i]]+=a;Normal[Face[i]+1]+=b;Normal[Face[i]+2]+=c;
		Normal[Face[i+1]]+=a;Normal[Face[i+1]+1]+=b;Normal[Face[i+1]+2]+=c;
		Normal[Face[i+2]]+=a;Normal[Face[i+2]+1]+=b;Normal[Face[i+2]+2]+=c;
	    
		//here NPoints stores the number of faces that contributed to the computation of a normal 
		NPoints[Face[i]/3]++;NPoints[Face[i+1]/3]++;NPoints[Face[i+2]/3]++;
	    
		//printf ("Points: %lg %lg %lg   %lg %lg %lg\n",x,y,z,u,v,w);
	    }
	}
	
	for (i=0;i<Object->NVertices*3;i+=3)
	{
	    //u=Normal[i];v=Normal[i+1];w=Normal[i+2];
	    //x=sqrt(u*u+v*v+w*w);
	    if (NPoints[i/3])
	    {
		Normal[i]/=(NPoints[i/3]);
		Normal[i+1]/=(NPoints[i/3]);
		Normal[i+2]/=(NPoints[i/3]);
	    }
	    //printf ("_%lg_ _%lg_ _%lg_  %d\n",Normal[i],Normal[i+1],Normal[i+2],NPoints[i/3]);
	}
	
	//printf("Nfaces:%d %d\n",np/3,Object->NFaces);
    }

    //Compute the norm of the gradiant, with sign
    if (flags&MCFLAG_GETGRAD)
    {
	float min=1.E13,max=-1.E13,dminmax;
	//int *VertexFlag;
	//VertexFlag = (int *) calloc(Object->NVertices,sizeof(int));

	OGrad = Object->grad;
	Object->gradNorm=(float *) realloc(Object->gradNorm,Object->NVertices*sizeof(float));
	for (i=0;i<Object->NVertices*3;i+=3)
	{
	    
	    //Sign of the gradiant (+ --> same as the normal)
	    /*
	    Object->gradNorm[i/3]=Normal[i]*OGrad[i]+
		Normal[i+1]*OGrad[i+1]+
		Normal[i+2]*OGrad[i+2];
	    Object->gradNorm[i/3]/=fabs(Object->gradNorm[i/3]);
	    */

	    Object->gradNorm[i/3]=OGrad[i]*OGrad[i]+
		OGrad[i+1]*OGrad[i+1]+
		OGrad[i+2]*OGrad[i+2];
	    
	}
	for (i=0;i<Object->NVertices;i++)
	{
	    if (min>Object->gradNorm[i]) min=Object->gradNorm[i];
	    if (max<Object->gradNorm[i]) max=Object->gradNorm[i];
	}
	dminmax= max-min;
	for (i=0;i<Object->NVertices;i++)
	{
	    Object->gradNorm[i]=(Object->gradNorm[i]-min)/dminmax;
	}
	//printf ("%lg %lg \n",min,max);
	//MUpdateMaxList(Object,VertexFlag,0);
	/*
	for (i=0,b=0,a=1000000;i<Object->NVertices;i++)
	{
	    if ((Point[3*i]<(density->Nx-3)*density->dx+density->x0)&&
		(Point[3*i+1]<(density->Ny-3)*density->dy+density->y0)&&
		(Point[3*i+2]<(density->Nz-3)*density->dz+density->z0)&&
		(Point[3*i]>(3)*density->dx+density->x0)&&
		(Point[3*i+1]>(3)*density->dy+density->y0)&&
		(Point[3*i+2]>(3)*density->dz+density->z0))
	    {
		if (Object->gradNorm[i]<a) a=Object->gradNorm[i];
		if (Object->gradNorm[i]>b) b=Object->gradNorm[i];
	    }
	}
	//b-=a;b=1;
	printf("a=%f b=%f\n",a,b);
	for (i=0;i<Object->NVertices;i++) Object->gradNorm[i]=(Object->gradNorm[i]-a)/b;
	*/
    }

    free (NPoints);

    return 0;
}

