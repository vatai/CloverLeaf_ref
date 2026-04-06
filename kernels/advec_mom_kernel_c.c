/*Crown Copyright 2012 AWE.
*
* This file is part of CloverLeaf.
*
* CloverLeaf is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 3 of the License, or (at your option)
* any later version.
*
* CloverLeaf is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* CloverLeaf. If not, see http://www.gnu.org/licenses/. */

/**
 *  @brief C momentum advection kernel
 *  @author Wayne Gaudin
 *  @details Performs a second order advective remap on the vertex momentum
 *  using van-Leer limiting and directional splitting.
 *  Note that although pre_vol is only set and not used in the update, please
 *  leave it in the method.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

#include <sys/time.h>
#define START_TIMER(i) gettimeofday(&__start, NULL);
#define END_TIMER(i) {gettimeofday(&__end, NULL); __timers[i] += (__end.tv_sec - __start.tv_sec) + (__end.tv_usec - __start.tv_usec) / 1000000.0;};
double __timers[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
struct timeval __start, __end;

void advec_mom_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                      double *vel1,
                      double *mass_flux_x,
                      double *vol_flux_x,
                      double *mass_flux_y,
                      double *vol_flux_y,
                      double *volume,
                      double *density1,
                      double *node_flux,
                      double *node_mass_post,
                      double *node_mass_pre,
                      double *mom_flux,
                      double *pre_vol,
                      double *post_vol,
                      double *celldx,
                      double *celldy,
                         int *whch_vl,
                         int *swp_nmbr,
                         int *drctn)

{
  // int x_min=*xmin;
  // int x_max=*xmax;
  // int y_min=*ymin;
  // int y_max=*ymax;
#define x_min 1
#define x_max 960
#define y_min 1
#define y_max 960
  int which_vel=*whch_vl;
  int sweep_number=*swp_nmbr;
  int direction=*drctn;

  int j,k,mom_sweep;
  int upwind,donor,downwind,dif;
  double sigma,wind,width;
  double vdiffuw,vdiffdw,auw,adw,limiter;

  double advec_vel_s;


  mom_sweep=direction+2*(sweep_number-1);
  
#pragma omp parallel
 {
  if(mom_sweep==1){
#pragma scop // 0
START_TIMER(0)
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                           -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                          +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                          -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 0
END_TIMER(0)
  } else if(mom_sweep==2){
#pragma scop // 1
START_TIMER(1)
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                                -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                               +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                               -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 1
END_TIMER(1)
  } else if(mom_sweep==3){
#pragma scop // 2
START_TIMER(2)
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                          +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                          -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 2
END_TIMER(2)
  } else if(mom_sweep==4){
#pragma scop // 3
START_TIMER(3)
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                          +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                          -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 3
END_TIMER(3)
  }

  if(direction==1) {
#pragma scop // 4
START_TIMER(4)
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma omp simd
      for (j=x_min-2;j<=x_max+2;j++) {
        node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                            *(mass_flux_x[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j+1,k-1,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma omp simd
      for (j=x_min-1;j<=x_max+2;j++) {
        node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                                 *(density1[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma omp simd
      for (j=x_min-1;j<=x_max+2;j++) {
        node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                          -node_flux[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]+node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 4
END_TIMER(4)
#pragma omp for private(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,j,advec_vel_s)
    for (k=y_min;k<=y_max+1;k++) {
      for (j=x_min-1;j<=x_max+1;j++) {
        if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
          upwind=j+2;
          donor=j+1;
          downwind=j;
          dif=donor;
        }
        else{
          upwind=j-1;
          donor=j;
          downwind=j+1;
          dif=upwind;
        }
        sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]);
        width=celldx[FTNREF1D(j,x_min-2)];
        vdiffuw=vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(upwind,k  ,x_max+5,x_min-2,y_min-2)];
        vdiffdw=vel1[FTNREF2D(downwind,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)];
        limiter=0.0;
        if(vdiffuw*vdiffdw>0.0){
          auw=fabs(vdiffuw);
          adw=fabs(vdiffdw);
          wind=1.0;
          if(vdiffdw<=0.0) wind=-1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx[FTNREF1D(dif,x_min-2)])/6.0,MIN(auw,adw));
        }
        advec_vel_s=vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
        mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel_s
                                                           *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }

#pragma scop // 5
START_TIMER(5)
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma omp simd
      for (j=x_min;j<=x_max+1;j++) {
        vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        *node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        +mom_flux[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]
                                                        -mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                                        /node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 5
END_TIMER(5)
  }
  else if(direction==2){
#pragma scop // 6
START_TIMER(6)
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min;j<=x_max+1;j++) {
        node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                            *(mass_flux_y[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j-1,k+1,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min-1;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min;j<=x_max+1;j++) {
        node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                                 *(density1[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min-1;k<=y_max+2;k++) {
#pragma omp simd
      for (j=x_min;j<=x_max+1;j++) {
        node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                          -node_flux[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]+node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 6
END_TIMER(6)
#pragma omp for private(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,j,advec_vel_s)
    for (k=y_min-1;k<=y_max+1;k++) {
      for (j=x_min;j<=x_max+1;j++) {
        if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
          upwind=k+2;
          donor=k+1;
          downwind=k;
          dif=donor;
        }
        else{
          upwind=k-1;
          donor=k;
          downwind=k+1;
          dif=upwind;
        }
        sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]);
        width=celldy[FTNREF1D(k,y_min-2)];
        vdiffuw=vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,upwind,x_max+5,x_min-2,y_min-2)];
        vdiffdw=vel1[FTNREF2D(j  ,downwind ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)];
        limiter=0.0;
        if(vdiffuw*vdiffdw>0.0){
          auw=fabs(vdiffuw);
          adw=fabs(vdiffdw);
          wind=1.0;
          if(vdiffdw<=0.0) wind=-1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy[FTNREF1D(dif,y_min-2)])/6.0,MIN(auw,adw));
        }
        advec_vel_s=vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
        mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel_s
                                                           *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma scop // 7
START_TIMER(7)
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma omp simd
      for (j=x_min;j<=x_max+1;j++) {
        vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        *node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        +mom_flux[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                        -mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                                        /node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma endscop // 7
END_TIMER(7)

  }

 }

 printf("@@@@ |");
 for(int __i; __i < sizeof(__timers)/sizeof(*__timers); ++__i)
	 printf("%i %lf |", __i, __timers[__i]);
 printf("\n");
}

