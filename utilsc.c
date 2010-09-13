#include <stdio.h>
#include <stdlib.h>

// qsort wrapper based on example from
// http://www.cplusplus.com/reference/clibrary/cstdlib/qsort/

int meshindex_(float*);

int compare_floats(const void * a, const void * b) {
  return ( *(float*)a - *(float*)b );
}

int gen_index_6Dfloat(float *a) {
  double r, st, sp, rsp, psi;
  
  rsp = a[0]*a[0] + a[1]*a[1];
  r = sqrt( rsp + a[2]*a[2] );
  rsp = sqrt(rsp);
  st = rsp/r;
  //sp = a[1]/rsp;
  psi = atan2(a[0],a[1]);

  return (r*100. + (1.1+st)*1000000. + (3.2+psi)*100000000.);
}

int compare_6Dfloats(const void * a, const void * b) {

  return ( gen_index_6Dfloat((float*)a) - gen_index_6Dfloat((float*)b) );
  //return ( meshindex_((float*)a) - meshindex_((float*)b) );
}

void quicksort6dfloats_(float *vectors, int *n) {
  qsort(vectors, *n, 6*sizeof(float), compare_6Dfloats);
}

/*
int meshindex(float *a) {
  float small, rsp, rp, st, ct, sp, cp;
  float rf, tf, pf;
  int irl, itl, ipl;

  zro = 1.e-9;

  rsp = a(1)*a(1) + a(2)*a(2);
  rp = sqrt( rsp + a(3)*a(3) );
  rsp = sqrt(rsp);

  if (rsp > small) {
    sp = a(2)/rsp;
    cp = a(1)/rsp;
  } else {
    sp = 0.;
    cp = 1.;
  }
  ipl = interppsi_(sp,cp,pf);

  st = rsp/rp;
  ct = a(3)/rp;
  itl = interpth_(ct,tf);

  irl = 
}
      rsp=sqrt(rsp)
       then
         cp=x/rsp
         sp=y/rsp
      else
         cp=1.
         sp=0.
      endif

      ipl=interppsi(sp,cp,pf)
      
c theta sin/cos
      st=rsp/rp
      ct=z/rp
      ithl=interpth(ct,thf)

      irl=irpre(1+int((rp-r(1))*rfac))

*/
