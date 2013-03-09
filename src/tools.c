#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "field.h"
#include "mcubes.h"
#include "tools.h"

//Diagonalise 3x3 symetric matrix M, L are the eigenvalues (L[0]>=L[1]>=L[2]) and vec the eigenvectors
//returns 0 if OK, true if degenerate
//                   (M0 M1 M2)
//M = (M0 .. M5) --> (M1 M3 M4)
//                   (M2 M4 M5)
int Diagonalise3x3(double *M,double *L,double *vec)
{
    double a,b,c,u,v,w;
    double B,C,D;
    double p,q,r;
    double phi;
    double tmp;
    int i;
    
    a=M[0];u=M[1];v=M[2];
    b=M[3];c=M[5];w=M[4];

    //if ((a==b)&&(b==c)&&(c==0)) return -1;

    //Compute characteristic polynome
    B=-(a+b+c);
    C=a*b+a*c+b*c-u*u-v*v-w*w;
    D=u*u*c+v*v*b+w*w*a-a*b*c-2*u*v*w;

    //solve 3rd degree equation
    p=(3.*C-B*B)/9.;
    q=(2.*B*B*B/27.-B*C/3.+D)/2.;
    
    if (q*q+p*p*p>=0) 
    {
      L[0]=L[1]=L[2]=0;
      return -1;
    }
    
   
    r=((q<0)?-1:1)*sqrt(fabs(p));
    phi = acos(q/(r*r*r))/3.;
	
    r*=2;//B/=3;
    L[0]=-r*cos(phi)-B/3.;
    L[1]=r*cos(PI/3-phi)-B/3.;
    L[2]=r*cos(PI/3+phi)-B/3.;

    
    
    if ((L[2])>(L[1])) {tmp=L[1];L[1]=L[2];L[2]=tmp;}
    if ((L[1])>(L[0])) {tmp=L[1];L[1]=L[0];L[0]=tmp;}
    if ((L[2])>(L[1])) {tmp=L[1];L[1]=L[2];L[2]=tmp;}
    
    //Some special case ...

    //matrix is already diagonal
    if ((u==0)&&(v==0)&&(w==0))
    {
	vec[0]=0;vec[1]=0;vec[2]=0;
	//associate each axis (x,y and z) with it s eigenvalues
	for (i=0;i<3;i++)
	{
	    if (fabs(L[i])>1.E-8)
	    {
		vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;
		if (fabs((a-L[i])/L[i])<1.E-4) vec[3*i+0]=1;
		else if (fabs((b-L[i])/L[i])<1.E-4) vec[3*i+1]=1;
		else if (fabs((c-L[i])/L[i])<1.E-4) vec[3*i+2]=1;
		else printf ("OOPSSSSSS\n");
	    }
	    else {L[i]=0;vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;}
	}
    }
    
    //We already have 1 eigenvalue (matrix is block diagonal)
    else if (((u==0)&&(v==0))||((v==0)&&(w==0))||((u==0)&&(w==0)))
    {
      double tmpa,tmpb,tmpu;
	int tmpd;
	//Finds which axis are already eigenvector (x,y or z)
	vec[0]=0;vec[1]=0;vec[2]=0;
	for (i=0;i<3;i++)
	{
	    if (fabs(L[i])>1.E-8)
	    {
		vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;
		//Check if we already have this eigenvector
		if (fabs((a-L[i])/L[i])<1.E-4) vec[3*i+0]=1;
		else if (fabs((b-L[i])/L[i])<1.E-4) vec[3*i+1]=1;
		else if (fabs((c-L[i])/L[i])<1.E-4) vec[3*i+2]=1;
		//no ...
		//diagonalise 2x2 block ...
		else
		{
		    if ((u==0)&&(v==0)) {tmpa=b;tmpb=c;tmpu=w;tmpd=1;}
		    if ((v==0)&&(w==0)) {tmpa=a;tmpb=b;tmpu=u;tmpd=0;}
		    if ((u==0)&&(w==0)) {tmpa=a;tmpb=c;tmpu=v;tmpd=2;}
		    if ((tmpa==tmpb)&&(tmpb==0)) 
			B=C=0;
		    else
		    {
			B=(L[i]+tmpu-tmpa)/(L[i]+tmpu-tmpb);
			C=1./sqrt(1+B*B);
		    }
		    if (tmpd==2)
		    {
			vec[3*i]=C;
			vec[3*i+2]=B*C;
		    }
		    else
		    {
			vec[3*i+tmpd]=C;
			vec[3*i+tmpd+1]=B*C;
		    }
		}
	    }
	    else {L[i]=0;vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;}
	}
    }

    //General case ...
    else for (i=0;i<3;i++)
    {

	B=L[i]-b-u*u/(L[i]-a);
	B=(w+u*v/(L[i]-a))/B;
	
	C=u/(L[i]-a)*B+v/(L[i]-a);
	D=1./sqrt(B*B+C*C+1.);

	vec[3*i+2]=D;
	vec[3*i+1]=B*D;
	vec[3*i]=C*D;
    }
  
    //Now rearrange vectors to have a direct eigenbasis
    /*
    B=vec[0]*(vec[3+1]*vec[6+2]-vec[3+2]*vec[6+1])-
	vec[3]*(vec[0+1]*vec[6+2]-vec[0+2]*vec[6+1])+
	vec[6]*(vec[0+1]*vec[3+2]-vec[0+2]*vec[3+1]);

    if (B<0) {vec[0]=-vec[0];vec[1]=-vec[1];vec[2]=-vec[2];}
    */
    
    
    return 0;
}

//Returns in inter the intersection of
//line defined  by two points p1 and p2
//and plane defined by it s normal n and a point p3
int Line2plane_Inter(float *p1,float *p2, float *p3, float *n,float *inter)
{
    double A,B;

    A=n[0]*(p3[0]-p1[0])+n[1]*(p3[1]-p1[1])+n[2]*(p3[2]-p1[2]);
    B=n[0]*(p2[0]-p1[0])+n[1]*(p2[1]-p1[1])+n[2]*(p2[2]-p1[2]);
   
    if (B==0) {printf ("Faces are coplanar, I don't deal with this case ....\n");return 0;}
    //if (fabs(B)<0.5*sqrt()) return 0;
    A/=B;

    inter[0]=p1[0]+A*(p2[0]-p1[0]);
    inter[1]=p1[1]+A*(p2[1]-p1[1]);
    inter[2]=p1[2]+A*(p2[2]-p1[2]);

    return 1;
}

//(p1,q1,r1) are points coordinates that define face 1, (p2,q2,r2) define face 2
//resulting intersection segment is returned in seg[6] ->[x1,y1,z1,x2,y2,z2]
//function returns 0 if no intersection 
//returns 1 if there is an intersection
//returns -1 if triangles are coplanar
//A Faire: coplanaire --> renvoie un truc special (devrait etre OK)
int Find_2Faces_Inter(float *p_p1, float *p_q1, float *p_r1,float *p_p2, float *p_q2, float *p_r2,float *seg)
{
    double a,b,c,d,e,f,g,h,i;
    float *p1,*p2,*q1,*q2,*r1,*r2,*tmp;
    double s1,s2,s3;
    //float u,v,w;

    p1=p_p1;p2=p_p2;
    q1=p_q1;q2=p_q2;
    r1=p_r1;r2=p_r2;

    //checks wether p2 is above/below/on triangle 1 (s1>0 / s1<0 / s1=0)
    a=p1[0]-p2[0];b=p1[1]-p2[1];c=p1[2]-p2[2];
    d=q1[0]-p2[0];e=q1[1]-p2[1];f=q1[2]-p2[2];
    g=r1[0]-p2[0];h=r1[1]-p2[1];i=r1[2]-p2[2];

    s1=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);

    //checks wether q2 is above/below/on triangle 1
    a=p1[0]-q2[0];b=p1[1]-q2[1];c=p1[2]-q2[2];
    d=q1[0]-q2[0];e=q1[1]-q2[1];f=q1[2]-q2[2];
    g=r1[0]-q2[0];h=r1[1]-q2[1];i=r1[2]-q2[2];

    s2=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);

    //checks wether r2 is above/below/on triangle 1
    a=p1[0]-r2[0];b=p1[1]-r2[1];c=p1[2]-r2[2];
    d=q1[0]-r2[0];e=q1[1]-r2[1];f=q1[2]-r2[2];
    g=r1[0]-r2[0];h=r1[1]-r2[1];i=r1[2]-r2[2];

    s3=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);
    
    //triangles are coplanar (or an edge/vertex is on the other's plane) 
    if ((s1==0)||(s2==0)||(s3==0))
    {
	//triangles are coplanar
	if ((fabs(s1)<5E8)&&(fabs(s2)<5E8)&&(fabs(s3)<5E8)) return -1;
	//one vertex is on the other triangle plane
	
	if ((s1==0)&&(s2!=0)&&(s3!=0)) {/*printf ("Case no handled so far ...\n");*/return 0;}
	if ((s1!=0)&&(s2==0)&&(s3!=0)) {/*printf ("Case no handled so far ...\n");*/return 0;}
	if ((s1!=0)&&(s2!=0)&&(s3==0)) {/*printf ("Case no handled so far ...\n");*/return 0;}

	//two vertices are on the other triangle plane, that's fine ...
    }

    //Triangles don't intersect (all vertices of face 2 are above/below face 1) 
    if (((s1>0)&&(s2>0)&&(s3>0))||((s1<0)&&(s2<0)&&(s3<0))) return 0;

    //They intersect ... maybe
    
    //p2 is set to be the point alone on its side of face 1
    //point of face 1 are eventually rearranged for p2 to be above it
    //in the sense defined above
    if (SIGN(s1)==SIGN(s2)) 
    {
	tmp=p2;p2=r2;r2=tmp;
	if (s3<0) {tmp=q1;q1=r1;r1=tmp;}
    }
    else if (SIGN(s1)==SIGN(s3)) 
    {
	tmp=p2;p2=q2;q2=tmp;
	if (s2<0) {tmp=q1;q1=r1;r1=tmp;}
    }
    else if (s1<0) {tmp=q1;q1=r1;r1=tmp;}
    
    //checks wether p1 is above/below/on triangle 2 (s1>0 / s1<0 / s1=0)
    a=p2[0]-p1[0];b=p2[1]-p1[1];c=p2[2]-p1[2];
    d=q2[0]-p1[0];e=q2[1]-p1[1];f=q2[2]-p1[2];
    g=r2[0]-p1[0];h=r2[1]-p1[1];i=r2[2]-p1[2];
    
    s1=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);
    
    //checks wether q1 is above/below/on triangle 2
    a=p2[0]-q1[0];b=p2[1]-q1[1];c=p2[2]-q1[2];
    d=q2[0]-q1[0];e=q2[1]-q1[1];f=q2[2]-q1[2];
    g=r2[0]-q1[0];h=r2[1]-q1[1];i=r2[2]-q1[2];
    
    s2=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);
    
    //checks wether r1 is above/below/on triangle 2
    a=p2[0]-r1[0];b=p2[1]-r1[1];c=p2[2]-r1[2];
    d=q2[0]-r1[0];e=q2[1]-r1[1];f=q2[2]-r1[2];
    g=r2[0]-r1[0];h=r2[1]-r1[1];i=r2[2]-r1[2];
    
    s3=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);
    
    //some special case
    if ((s1==0)||(s2==0)||(s3==0))
    {
	//triangles are coplanar
	if ((s1==0)&&(s2==0)&&(s3==0)) return -1;
	//one vertex is on the other triangle plane
	if ((s1==0)&&(s2!=0)&&(s3!=0)) {/*printf ("Case no handled so far ...\n");*/return 0;}
	if ((s1!=0)&&(s2==0)&&(s3!=0)) {/*printf ("Case no handled so far ...\n");*/return 0;}
	if ((s1!=0)&&(s2!=0)&&(s3==0)) {/*printf ("Case no handled so far ...\n");*/return 0;}
	//two vertices are on the other triangle plane, that's fine ...
    }

    //Triangles don't intersect (all vertices of face 1 are above/below face 2) 
    else if (((s1>0)&&(s2>0)&&(s3>0))||((s1<0)&&(s2<0)&&(s3<0))) return 0;

    //p1 is set to be the point alone on its side of face 2
    //point of face 2 are eventually rearranged for p1 to be above it
    //in the sense defined above
    if (SIGN(s1)==SIGN(s2)) 
    {
	tmp=p1;p1=r1;r1=tmp;
	tmp=q1;q1=r1;r1=tmp;
	if (s3<0) {tmp=q2;q2=r2;r2=tmp;}
    }
    else if (SIGN(s1)==SIGN(s3)) 
    {
	tmp=p1;p1=q1;q1=tmp;
	tmp=q1;q1=r1;r1=tmp;
	if (s2<0) {tmp=q2;q2=r2;r2=tmp;}
    }
    else if (s1<0) {tmp=q2;q2=r2;r2=tmp;}
    
    //[p1,r1] and [p1,q1] intersect face 2 in i and j
    //[p2,r2] and [p2,q2] intersect face 1 in l and k
    //[i,j] and [k,l] are intersection of the triangles
    //with the intersection of the planes (L) they define
    //Because of the reorientation of the faces, the condition 
    //of intersection of the triangles is k<=j AND i<=l
    //on the properly oriented intersection L 

    a=p1[0]-p2[0];b=p1[1]-p2[1];c=p1[2]-p2[2];
    d=q1[0]-p2[0];e=q1[1]-p2[1];f=q1[2]-p2[2];
    g=r2[0]-p2[0];h=r2[1]-p2[1];i=r2[2]-p2[2];

    s1=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);

    if (s1>0)
    {
	a=p1[0]-p2[0];b=p1[1]-p2[1];c=p1[2]-p2[2];
	d=r1[0]-p2[0];e=r1[1]-p2[1];f=r1[2]-p2[2];
	g=r2[0]-p2[0];h=r2[1]-p2[1];i=r2[2]-p2[2];
	
	s2=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);

	//no intersection
	if (s2>0) return 0;
	
	a=p1[0]-p2[0];b=p1[1]-p2[1];c=p1[2]-p2[2];
	d=r1[0]-p2[0];e=r1[1]-p2[1];f=r1[2]-p2[2];
	g=q2[0]-p2[0];h=q2[1]-p2[1];i=q2[2]-p2[2];

	s3=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);
	
	
	//intersection is [i,l]
	if (s3>0) 
	{
	    float normal[3];

	    //printf ("i,l :\n");

	    a=q2[0]-p2[0];b=q2[1]-p2[1];c=q2[2]-p2[2];
	    d=r2[0]-p2[0];e=r2[1]-p2[1];f=r2[2]-p2[2];
	    normal[0]= b*f-c*e;
	    normal[1]= c*d-a*f;
	    normal[2]= a*e-b*d;
	    Line2plane_Inter(p1,r1,p2,normal,seg);
	    
	    a=q1[0]-p1[0];b=q1[1]-p1[1];c=q1[2]-p1[2];
	    d=r1[0]-p1[0];e=r1[1]-p1[1];f=r1[2]-p1[2];
	    normal[0]= b*f-c*e;
	    normal[1]= c*d-a*f;
	    normal[2]= a*e-b*d;
	    Line2plane_Inter(p2,r2,p1,normal,&seg[3]);
	}
	//intersection is [k,l]
	else 
	{
	    float normal[3];

	    //printf ("k,l :\n");

	    a=q1[0]-p1[0];b=q1[1]-p1[1];c=q1[2]-p1[2];
	    d=r1[0]-p1[0];e=r1[1]-p1[1];f=r1[2]-p1[2];
	    normal[0]= b*f-c*e;
	    normal[1]= c*d-a*f;
	    normal[2]= a*e-b*d;
	    
	    Line2plane_Inter(p2,q2,p1,normal,seg);
	    Line2plane_Inter(p2,r2,p1,normal,&seg[3]);
	}
    }
    else
    {
	a=p1[0]-p2[0];b=p1[1]-p2[1];c=p1[2]-p2[2];
	d=q1[0]-p2[0];e=q1[1]-p2[1];f=q1[2]-p2[2];
	g=q2[0]-p2[0];h=q2[1]-p2[1];i=q2[2]-p2[2];
	
	s2=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);

	//no intersection
	if (s2<0) return 0;

	a=p1[0]-p2[0];b=p1[1]-p2[1];c=p1[2]-p2[2];
	d=r1[0]-p2[0];e=r1[1]-p2[1];f=r1[2]-p2[2];
	g=q2[0]-p2[0];h=q2[1]-p2[1];i=q2[2]-p2[2];
	
	s3=a*(e*i-f*h) - b*(d*i-g*f) + c*(d*h-g*e);

	//intersection is [k,j]
	if (s3<0)
	{
	    float normal[3];
	    
	    //printf ("k,j :\n");

	    a=q1[0]-p1[0];b=q1[1]-p1[1];c=q1[2]-p1[2];
	    d=r1[0]-p1[0];e=r1[1]-p1[1];f=r1[2]-p1[2];
	    normal[0]= b*f-c*e;
	    normal[1]= c*d-a*f;
	    normal[2]= a*e-b*d;
	    
	    Line2plane_Inter(p2,q2,p1,normal,seg);

	    a=q2[0]-p2[0];b=q2[1]-p2[1];c=q2[2]-p2[2];
	    d=r2[0]-p2[0];e=r2[1]-p2[1];f=r2[2]-p2[2];
	    normal[0]= b*f-c*e;
	    normal[1]= c*d-a*f;
	    normal[2]= a*e-b*d;
	    Line2plane_Inter(p1,q1,p2,normal,&seg[3]);
	   
	}
	//intersection is [i,j]
	else 
	{
	    float normal[3];

	    //printf ("i,j :\n");

	    a=q2[0]-p2[0];b=q2[1]-p2[1];c=q2[2]-p2[2];
	    d=r2[0]-p2[0];e=r2[1]-p2[1];f=r2[2]-p2[2];
	    normal[0]= b*f-c*e;
	    normal[1]= c*d-a*f;
	    normal[2]= a*e-b*d;

	    Line2plane_Inter(p1,r1,p2,normal,seg);
	    Line2plane_Inter(p1,q1,p2,normal,&seg[3]);
	}
    }

    //printf ("OK\n");

    if ((seg[3]-seg[0])+(seg[4]-seg[1])+(seg[5]-seg[2]) == 0) return 0;

    return 1;
}

//returns in ipt the intersection of segment s1 and s2
// returns 0 if no intersect., -1 if colinear and 1 if they intersect
// algo is explained on http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
int SegIntersect(float *seg1,float *seg2,float *ipt)
{
    double a1,b1,c1,a2,b2,c2;
    double u1,v1,w1;
    double u2,v2,w2;
    double rx,ry,sy;
    double k,l;
    


    a1=seg1[3]-seg1[0];b1=seg1[3+1]-seg1[1];c1=seg1[3+2]-seg1[2];
    a2=seg2[3]-seg2[0];b2=seg2[3+1]-seg2[1];c2=seg2[3+2]-seg2[2];

    u1=b1*c2-c1*b2;v1=a2*c1-a1*c2;w1=a1*b2-b1*a2;

    //lines are not colinear ...
    if ((u1!=0)||(v1!=0)||(w1!=0))
    {
      if ((fabs(a1)>=fabs(b1))&&(fabs(a1)>=fabs(c1))) {u2=-(b1+c1)/a1;v2=1;w2=1;}
	else if ((fabs(b1)>=fabs(a1))&&(fabs(b1)>=fabs(c1))) {u2=1;v2=-(a1+c1)/b1;w2=1;}
	else {u2=1;v2=1;w2=-(a1+b1)/c1;}
	
	ry = sqrt(u2*u2+v2*v2+w2*w2);
	rx = sqrt(a1*a1+b1*b1+c1*c1);
	
	u2/=ry;v2/=ry;w2/=ry;
	ry=u2*(seg2[0]-seg1[0])+v2*(seg2[1]-seg1[1])+w2*(seg2[2]-seg1[2]);

	sy=ry+a2*u2+b2*v2+c2*w2;

	l=-(ry)/(sy-ry);
	
	if ((l<0)||(l>1)) return 0;
	u1=seg2[0]+l*a2;
	v1=seg2[1]+l*b2;
	w1=seg2[2]+l*c2;

	ry = rx*sqrt((u1-seg1[0])*(u1-seg1[0])+(v1-seg1[1])*(v1-seg1[1])+(w1-seg1[2])*(w1-seg1[2]));

	u2=(a1*(v1-seg1[1])-b1*(u1-seg1[0]))/ry; 
	v2=(a1*(w1-seg1[2])-c1*(u1-seg1[0]))/ry; 
	w2=(b1*(w1-seg1[2])-c1*(v1-seg1[1]))/ry;
	
	if ((fabs(u2)>2.E-2)||(fabs(v2)>2.E-2)||(fabs(w2)>2.E-2)) return 0;
	
	k=((u1-seg1[0])*a1+(v1-seg1[1])*b1+(w1-seg1[2])*c1)/(rx*rx);
	
	if ((k<=1)&&(k>=0))
	{
	    ipt[0]=u1;ipt[1]=v1;ipt[2]=w1;
	    return 1;
	}
	else return 0;

    }
    else 
	return -1;
}


//gives the density at a point of coordinates x,y,z
double GetPointDensity(density_grid *density,float x,float y,float z)
{
  double v;
  unsigned int nx,ny,nz;
  
/*
    i=(x-density->x0)/density->dx;
    j=(y-density->y0)/density->dy;
    k=(z-density->z0)/density->dz;
    
    dx=1;
    dy=density->Nx;
    dz=density->Ny*dy;
*/
    nx=density->Nx;
    ny=density->Ny;
    nz=density->Nz;
   
    //  l=i+j*dy+k*dz;

    InterpolateVector(density->grid,&v,1,nx,ny,nz,
		      (x-density->x0)/density->dx,
		      (y-density->y0)/density->dy,
		      (z-density->z0)/density->dz);

    return v;

    /*
    if (k==density->Nz-1) dz=-dz*(density->Nz-1);
    if (j==density->Ny-1) dy=-dy*(density->Ny-1);
    if (i==density->Nx-1) dx=-dx*(density->Nx-1);
	

    u=(x-(double)i*density->dx-density->x0)/density->dx;
    u1=u*density->grid[l+dx]+(1-u)*density->grid[l];
    u2=u*density->grid[l+dx+dy]+(1-u)*density->grid[l+dy];
    u=(y-(double)j*density->dy-density->y0)/density->dy;
    v1=u*u2+(1-u)*u1;
    
    u=(x-(double)i*density->dx-density->x0)/density->dx;
    u1=u*density->grid[l+dx+dz]+(1-u)*density->grid[l+dz];
    u2=u*density->grid[l+dx+dy+dz]+(1-u)*density->grid[l+dy+dz];
    u=(y-(double)j*density->dy-density->y0)/density->dy;
    v2=u*u2+(1-u)*u1;
    
    u=(z-(double)k*density->dz-density->z0)/density->dz;
    v=u*v2+(1-u)*v1;
    printf ("v=%f\n",v);
    

    return v;
    */
}


//returns int result the interpolated value of val 
//du, dv and dw are 0<..<1
//dx,dy and dz are  the increments to go to the next point along x,y or z;
int InterpolateVector(FLOAT *p_val,double *result,int size,int nx,int ny,int nz,double pp_ddu,double pp_ddv,double pp_ddw)
{
  int j;
  double du,dv,dw;
  double ddu,ddv,ddw;
  double p_ddu,p_ddv,p_ddw;
  int dx,dy,dz;
  int a,b,c;
  FLOAT *val;
  

  p_ddu = pp_ddu;
  p_ddv = pp_ddv;
  p_ddw = pp_ddw;

  if (p_ddu>=nx)
    while(p_ddu>=nx) p_ddu-=nx;
  if (p_ddv>=ny)
    while(p_ddv>=ny) p_ddv-=ny;
  if (p_ddw>=nz)
    while(p_ddw>=nz) p_ddw-=nz;

  if (p_ddu<0)
    while(p_ddu<0) p_ddu+=nx;
  if (p_ddv<0)
    while(p_ddv<0) p_ddv+=ny;
  if (p_ddw<0)
    while(p_ddw<0) p_ddw+=nz;
  
    
  a=(int)p_ddu;
  b=(int)p_ddv;
  c=(int)p_ddw;
    
  ddu = p_ddu-a;
  ddv = p_ddv-b;
  ddw = p_ddw-c;

  dx = size;
  dy = nx*dx;
  dz = dy*ny;

  val = &p_val[a*dx+b*dy+c*dz];

  if (a>=nx-1) dx *= -(nx-1);
  if (b>=ny-1) dy *= -(ny-1);
  if (c>=nz-1) dz *= -(nz-1);

  if (ddu<0) du=0;
  else if (ddu>1) du=1;
  else du = ddu;

  if (ddv<0) dv=0;
  else if (ddv>1) dv=1;
  else dv = ddv;

  if (ddw<0) dw=0;
  else if (ddw>1) dw=1;
  else dw = ddw;

  for (j=0;j<size;j++)
    result[j] = (1.-du)*(1.-dv)*(1.-dw)*val[j]+(du)*(1.-dv)*(1.-dw)*val[j+dx]+
      (1.-du)*(dv)*(1.-dw)*val[j+dy]+(1.-du)*(1.-dv)*(dw)*val[j+dz]+
      (du)*(dv)*(1.-dw)*val[j+dx+dy]+(1.-du)*(dv)*(dw)*val[j+dy+dz]+
      (du)*(1.-dv)*(dw)*val[j+dx+dz]+(du)*(dv)*(dw)*val[j+dx+dy+dz];

  return 1;
}

int ExtractDensityBlock(density_grid *density,density_grid **p_sub_density,int imin,int imax,int jmin,int jmax,int kmin,int kmax)
{
    density_grid *sub_density=*p_sub_density;

    int di,dj,dk;
    int is,js,ks;
    int i,j,k;

    int n,ns;
    
    di = 1+imax-imin;
    dj = 1+jmax-jmin;
    dk = 1+kmax-kmin;

    if (*p_sub_density != NULL)
    {
	free(sub_density->grid);
	free(sub_density->grad);
	free(sub_density->hessian);
	free(*p_sub_density);
	*p_sub_density=NULL;
    }

    sub_density = calloc(1,sizeof(density_grid));
    sub_density->grad=sub_density->hessian=NULL;
    sub_density->Nx=di;sub_density->Ny=dj;sub_density->Nz=dk;
    sub_density->NNodes = di*dj*dk;
    sub_density->x0 = density->x0+density->dx*imin;
    sub_density->y0 = density->y0+density->dy*jmin;
    sub_density->z0 = density->z0+density->dz*kmin;
    sub_density->dx = density->dx;
    sub_density->dy = density->dy;
    sub_density->dz = density->dz;
    
    sub_density->grid = calloc (di*dj*dk,sizeof(double));
    //printf ("comp ...");fflush(0);
    for (ks=0,ns=0;ks<dk;ks++)
	for (js=0;js<dj;js++)
	    for (is=0;is<di;is++,ns++)
	    {
		i=is+imin;
		j=js+jmin;
		k=ks+kmin;

		if (i<0) i+= density->Nx;
		else if (i>=density->Nx) i-= (density->Nx);
	
		if (j<0) j+= density->Ny;
		else if (j>=density->Ny) j-= (density->Ny);

		if (k<0) k+= density->Nz;
		else if (k>=density->Nz) k-= (density->Nz);

		n=i+j*density->Nx+k*(density->Nx*density->Ny);

		sub_density->grid[ns] = density->grid[n];

	    }
    *p_sub_density = sub_density;

    return 1;
}

//given 3 surfaces s1,s2 and s3, finds the N intersection points and return their coordinates in
//itrlist[3*N].
//returns the number of points found
int Find_3Surf_Inter_Points(MC_Object *s1,MC_Object *s2,MC_Object *s3,float **itrlist)
{
    int i,j,k;
    int a,b,c;
    int n,p1,p2,p3;
    int L1[10],L2[10],L3[10]; //List of the faces in the same box
    float seg1[6];
    float seg2[6];
    float seg3[6];
    float *P1,*P2,*P3;//pointer on points 
    int *F1,*F2,*F3;
    int i1,i2,i3;
    int count=0;

    i=j=k=0;

    if ((!s1->NFaces)||(!s2->NFaces)||(!s3->NFaces)) return 0;

    n=s1->BoxIndex[0];
    if (s2->BoxIndex[0]>n) n=s2->BoxIndex[0];
    if (s3->BoxIndex[0]>n) n=s3->BoxIndex[0];
 
    P1=s1->Vertex;F1=s1->Face;
    P2=s2->Vertex;F2=s2->Face;
    P3=s3->Vertex;F3=s3->Face;

    //Just checks if there exist a face in nth cell in every MC_Object
    while ((i<s1->NFaces)&&(j<s2->NFaces)&&(k<s3->NFaces))
    {
	for (;(i<s1->NFaces)&&(s1->BoxIndex[i]<n);i++);
	for (;(j<s2->NFaces)&&(s2->BoxIndex[j]<n);j++);
	for (;(k<s3->NFaces)&&(s3->BoxIndex[k]<n);k++);

	//Checks if there is a face from every MC_Object
	if ((i<s1->NFaces)&&(j<s2->NFaces)&&(k<s3->NFaces))
	if ((s1->BoxIndex[i]==n)&&(s2->BoxIndex[j]==n)&&(s3->BoxIndex[k]==n))
	{
	    for (p1=0;(i<s1->NFaces)&&(s1->BoxIndex[i]==n);i++,p1++) {s1->BoxIndex[i]=-s1->BoxIndex[i];L1[p1]=3*i;}
	    for (p2=0;(j<s2->NFaces)&&(s2->BoxIndex[j]==n);j++,p2++) {s2->BoxIndex[j]=-s2->BoxIndex[j];L2[p2]=3*j;}
	    for (p3=0;(k<s3->NFaces)&&(s3->BoxIndex[k]==n);k++,p3++) {s3->BoxIndex[k]=-s3->BoxIndex[k];L3[p3]=3*k;}
	    
	    //Now look for the intersection of these faces
	    for (c=0;c<p3;c++)
		for (b=0;b<p2;b++)
		{
		    i2=L2[b];i3=L3[c];
		    if (Find_2Faces_Inter(&P2[F2[i2]],&P2[F2[i2+1]],&P2[F2[i2+2]],&P3[F3[i3]],&P3[F3[i3+1]],&P3[F3[i3+2]],seg1)>0)
			for (a=0;a<p1;a++)
			{
			    i1=L1[a];
			    if ((Find_2Faces_Inter(&P2[F2[i2]],&P2[F2[i2+1]],&P2[F2[i2+2]],
						   &P1[F1[i1]],&P1[F1[i1+1]],&P1[F1[i1+2]],seg2)>0)&&
			    	(Find_2Faces_Inter(&P1[F1[i1]],&P1[F1[i1+1]],&P1[F1[i1+2]],
						   &P3[F3[i3]],&P3[F3[i3+1]],&P3[F3[i3+2]],seg3)>0))
			    {
				float segint1[3],segint2[3],segint3[3];
			
				if ((SegIntersect(seg1,seg2,segint1)>0)&&
				    (SegIntersect(seg2,seg3,segint2)>0)&&
				    (SegIntersect(seg1,seg3,segint3)>0))
				{
				    *itrlist = (float *) realloc (*itrlist,3*sizeof(float)*(count+1));
				    (*itrlist)[3*count]=segint1[0];
				    (*itrlist)[3*count+1]=segint1[1];
				    (*itrlist)[3*count+2]=segint1[2];
				    count++;
				    //printf ("Got One more (%d)...\n",count);
				    c=p3;b=p2;a=p1;
				}
			    }
			}
		}

	}

	if (i<s1->NFaces) n=s1->BoxIndex[i];
	if ((j<s2->NFaces)&&(s2->BoxIndex[j]>n)) n=s2->BoxIndex[j];
	if ((k<s3->NFaces)&&(s3->BoxIndex[k]>n)) n=s3->BoxIndex[k];
    }
    
    return count;
}
