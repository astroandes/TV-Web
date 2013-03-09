#ifndef FFT_DIFF_H
#define FFT_DIFF_H

#ifndef PI
#define PI 3.141592653589793238462643383279502884
#endif

#ifdef SIMPLEPRECISION
#define FLOAT float
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN fftwf_plan
#define FFTWNAME(x) fftwf_##x
#else
#define FLOAT double
#define FFTW_COMPLEX fftw_complex
#define FFTW_PLAN fftw_plan
#define FFTWNAME(x) fftw_##x
#endif

#define GRAV 6.6742E-11
#define RHO0 (0.3*10E-29)

//differentiation
#define FLAG_DX (1<<0)
#define FLAG_DY (1<<1)
#define FLAG_DZ (1<<2)
#define FLAG_D2X (1<<3)
#define FLAG_D2Y (1<<4)
#define FLAG_D2Z (1<<5)
#define FLAG_D3X (1<<6)
#define FLAG_D3Y (1<<7)
#define FLAG_D3Z (1<<8)
#define FLAG_DD (1<<9)  /*Laplacian*/

//Integration
#define FLAG_IX (1<<10)
#define FLAG_IY (1<<11)
#define FLAG_IZ (1<<12)
#define FLAG_I2X (1<<13)
#define FLAG_I2Y (1<<14)
#define FLAG_I2Z (1<<15)
#define FLAG_I3X (1<<16)
#define FLAG_I3Y (1<<17)
#define FLAG_I3Z (1<<18)
#define FLAG_II (1<<19)  /*Laplacian*/

int FFT_diff(FLOAT *,int,int,int,double,double,double,int,FLOAT **);

#endif
