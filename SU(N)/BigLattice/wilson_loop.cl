/******************************************************************************
 * @file     wilson_loop.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Measurement of the Wilson loop on lattice devided into parts
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2017 Vadim Demchik, Natalia Kolomoyets
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *    Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *
 *    Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 *****************************************************************************/
//BigLattice!!!
#ifndef WILSON_LOOP_CL
#define WILSON_LOOP_CL

#include "complex.h"
#include "model.cl"
#include "misc.cl"
#if SUN == 2
#include "su2cl.cl"
#include "su2_matrix_memory.cl"
#include "su2_measurements_cl.cl"
#endif
#if SUN == 3
#include "su3cl.cl"
#include "su3_matrix_memory.cl"
#include "su3_measurements_cl.cl"
#endif

                    HGPU_INLINE_PREFIX_VOID void
lattice_zero3(gpu_su_3* matrix)
{
    (*matrix).uv1 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
    (*matrix).uv2 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
    (*matrix).uv3 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_Lt_store_3(__global hgpu_float4 * lattice_table, gpu_su_3* m, uint gindex){
            lattice_table[gindex +  0 * ROWSIZE_LT] = (*m).uv1;
            lattice_table[gindex +  1 * ROWSIZE_LT] = (*m).uv2;
            lattice_table[gindex +  2 * ROWSIZE_LT] = (*m).uv3;
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_Lr2_store_3(__global hgpu_float4 * lattice_table, gpu_su_3* m, uint gindex){
            lattice_table[gindex +  0 * ROWSIZE_LR2] = (*m).uv1;
            lattice_table[gindex +  1 * ROWSIZE_LR2] = (*m).uv2;
            lattice_table[gindex +  2 * ROWSIZE_LR2] = (*m).uv3;
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_L_store_3(__global hgpu_float4 * lattice_table, gpu_su_3* m, uint gindex, int rowsize){
            lattice_table[gindex +  0 * rowsize] = (*m).uv1;
            lattice_table[gindex +  1 * rowsize] = (*m).uv2;
            lattice_table[gindex +  2 * rowsize] = (*m).uv3;
}

                    HGPU_INLINE_PREFIX gpu_su_3
lattice_Lt_get_3(__global hgpu_float4 * lattice_table, uint gindex){
    gpu_su_3 matrix;
            matrix.uv1 = lattice_table[gindex +  0 * ROWSIZE_LT];
            matrix.uv2 = lattice_table[gindex +  1 * ROWSIZE_LT];
            matrix.uv3 = lattice_table[gindex +  2 * ROWSIZE_LT];
    return matrix;
}

                    HGPU_INLINE_PREFIX gpu_su_3
lattice_Lr2_get_3(__global hgpu_float4 * lattice_table, uint gindex){
    gpu_su_3 matrix;
            matrix.uv1 = lattice_table[gindex +  0 * ROWSIZE_LR2];
            matrix.uv2 = lattice_table[gindex +  1 * ROWSIZE_LR2];
            matrix.uv3 = lattice_table[gindex +  2 * ROWSIZE_LR2];
    return matrix;
}

                    HGPU_INLINE_PREFIX gpu_su_3
lattice_L_get_3(__global hgpu_float4 * lattice_table, uint gindex, int rowsize){
    gpu_su_3 matrix;
            matrix.uv1 = lattice_table[gindex +  0 * rowsize];
            matrix.uv2 = lattice_table[gindex +  1 * rowsize];
            matrix.uv3 = lattice_table[gindex +  2 * rowsize];
    return matrix;
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_zero2(gpu_su_2* matrix)
{
    (*matrix).uv1 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_L_store_2(__global hgpu_float4 * lattice_table, gpu_su_2* m, uint gindex){
            lattice_table[gindex] = (*m).uv1;
}

                    HGPU_INLINE_PREFIX gpu_su_2
lattice_L_get_2(__global hgpu_float4 * lattice_table, uint gindex){
    gpu_su_2 matrix;
            matrix.uv1 = lattice_table[gindex];
    return matrix;
}

                                       __kernel void
lattice_measurement_Lt(__global hgpu_float4  * lattice_table,
                       __global hgpu_float4  * Lt,
                       __global hgpu_float   * lattice_parameters){
#if SUN == 2
    gpu_su_2 gLt, gLt2, gLt3;
    su2_twist twist;
        twist.phi   = lattice_parameters[1];
#endif
#if SUN == 3
    gpu_su_3 gLt, gLt2, gLt3;
    su3_twist twist;
        twist.phi   = lattice_parameters[1];
        twist.omega = lattice_parameters[2];
#endif
    coords_4 coord, coord1, coord2;
    uint gid2;
    int j;

    const uint gid = GID;
    if(GID < (N1 + 1) * N2 * N3 * N4){
    lattice_gid_to_coords(&gid,&coord);
#if SUN == 2
    gLt = lattice_table_2(lattice_table, &coord, GID, T, &twist);
#endif
#if SUN == 3
    gLt = lattice_table_3(lattice_table, &coord, GID, T, &twist);
#endif

    coord2 = coord;
    for (j = 1; j < WLT; j++){
        coord1 = coord2;
        lattice_neighbours_gid(&coord1, &coord2, &gid2, T);
#if SUN == 2
        gLt2 = lattice_table_2(lattice_table, &coord2, gid2, T, &twist);
        gLt3 = matrix_times2(&gLt, &gLt2);
        gLt = gLt3;
#endif
#if SUN == 3
        gLt2 = lattice_table_3(lattice_table, &coord2, gid2, T, &twist);
        gLt3 = matrix_times3(&gLt, &gLt2);
        gLt = gLt3;
#endif
    }
#if SUN == 2
    lattice_L_store_2(Lt, &gLt, gid);
#endif
#if SUN == 3
#if PREC == 2
    lattice_Lt_store_3_double(Lt, &gLt, gid);
#else
    lattice_Lt_store_3(Lt, &gLt, gid);
#endif
#endif
    }
}

                                       __kernel void
lattice_measurement_Lr(__global hgpu_float4  * lattice_table,
                       __global hgpu_float4  * Lr1,
                       __global hgpu_float4  * Lr2,
                       __global hgpu_float   * lattice_parameters)
{
    coords_4 coord, coord1, coord2, coord0;
    int M, j, i;
    uint gid2, gid0;
#if SUN == 2
    gpu_su_2 gLr1, gLr3, gLr4, gLr5;
    gpu_su_2 gLr2[WLN];
    su2_twist twist;
        twist.phi   = lattice_parameters[1];
#endif
#if SUN == 3
    gpu_su_3 gLr1, gLr3, gLr4, gLr5;
    gpu_su_3 gLr2[WLN];
    su3_twist twist;
        twist.phi   = lattice_parameters[1];
        twist.omega = lattice_parameters[2];
#endif

    const uint gid = GID;

    if(GID < SITES){
        lattice_gid_to_coords(&gid,&coord);
        M = min(WLR, (int)(N1 - coord.x));

#if SUN == 2
        gLr1 = lattice_table_2(lattice_table, &coord, GID, X, &twist);
#endif
#if SUN == 3
        gLr1 = lattice_table_3(lattice_table, &coord, GID, X, &twist);
#endif
        coord2 = coord;
        for(j = 1; j < M; j++){
            coord1 = coord2;
            lattice_neighbours_gid(&coord1, &coord2, &gid2, X);
#if SUN == 2
            gLr4 = lattice_table_2(lattice_table, &coord2, gid2, X, &twist);
            gLr3 = matrix_times2(&gLr1, &gLr4);
            gLr1 = gLr3;
#endif
#if SUN == 3
            gLr4 = lattice_table_3(lattice_table, &coord2, gid2, X, &twist);
            gLr3 = matrix_times3(&gLr1, &gLr4);
            gLr1 = gLr3;
#endif
        }

        coord0.x = 0;
        coord0.y = coord.y;     coord0.z = coord.z;     coord0.t = coord.t;
        lattice_coords_to_gid(&gid0, &coord0);
        coord2 = coord0;
        for(i = 0; i < WLN; i++){
#if SUN == 2
            gLr5 = lattice_table_2(lattice_table, &coord0, gid0, X, &twist);
#endif
#if SUN == 3
            gLr5 = lattice_table_3(lattice_table, &coord0, gid0, X, &twist);
#endif
            for(j = 1; j < i; j++){
                coord1 = coord2;
                lattice_neighbours_gid(&coord1, &coord2, &gid2, X);
#if SUN == 2
                gLr4 = lattice_table_2(lattice_table, &coord2, gid2, X, &twist);
                gLr3 = matrix_times2(&gLr5, &gLr4);
                gLr5 = gLr3;
#endif
#if SUN == 3
                gLr4 = lattice_table_3(lattice_table, &coord2, gid2, X, &twist);
                gLr3 = matrix_times3(&gLr5, &gLr4);
                gLr5 = gLr3;
#endif
            }
            gLr2[i] = gLr5;
        }

#if SUN == 2
        lattice_L_store_2(Lr1, &gLr1, gid);
        for(i = 0; i < WLN; i++)
            lattice_L_store_2(Lr2, &gLr2[i], gid + i * 1 * ROWSIZE_LT);
#endif
#if SUN == 3
        lattice_Lt_store_3(Lr1, &gLr1, gid);
        for(i = 0; i < WLN; i++)
            lattice_Lt_store_3(Lr2, &gLr2[i], gid + i * 3 * ROWSIZE_LT);
#endif
    }
}

                                       __kernel void
lattice_measurement_Lr1(__global hgpu_float4  * lattice_table,
                       __global hgpu_float4  * Lr1,
                       __global hgpu_float   * lattice_parameters)
{
    coords_4 coord, coord1, coord2;
    int M, j;
    uint gid2;
#if SUN == 2
    gpu_su_2 gLr1, gLr3, gLr4;
    su2_twist twist;
        twist.phi   = lattice_parameters[1];
#endif
#if SUN == 3
    gpu_su_3 gLr1, gLr3, gLr4;
    su3_twist twist;
        twist.phi   = lattice_parameters[1];
        twist.omega = lattice_parameters[2];
#endif

    const uint gid = GID;

    if(GID < SITES){
        lattice_gid_to_coords(&gid,&coord);
        M = min(WLR, (int)(N1 - coord.x));

#if SUN == 2
        gLr1 = lattice_table_2(lattice_table, &coord, GID, X, &twist);
#endif
#if SUN == 3
        gLr1 = lattice_table_3(lattice_table, &coord, GID, X, &twist);
#endif
        coord2 = coord;
        for(j = 1; j < M; j++){
            coord1 = coord2;
            lattice_neighbours_gid(&coord1, &coord2, &gid2, X);
#if SUN == 2
            gLr4 = lattice_table_2(lattice_table, &coord2, gid2, X, &twist);
            gLr3 = matrix_times2(&gLr1, &gLr4);
#endif
#if SUN == 3
            gLr4 = lattice_table_3(lattice_table, &coord2, gid2, X, &twist);
            gLr3 = matrix_times3(&gLr1, &gLr4);
#endif
            gLr1 = gLr3;
        }

#if SUN == 2
        lattice_L_store_2(Lr1, &gLr1, gid);
#endif
#if SUN == 3
        lattice_L_store_3(Lr1, &gLr1, gid, ROWSIZE_WLX);
#endif
    }
}

                                       __kernel void
lattice_measurement_Lr2(__global hgpu_float4  * lattice_table,
                       __global hgpu_float4  * Lr2,
                       __global hgpu_float   * lattice_parameters){
    coords_4 coord1, coord2, coord0;
    int j, i;
    uint gid2, gid0;
#if SUN == 2
    gpu_su_2 gLr3, gLr4, gLr5;
    su2_twist twist;
        twist.phi   = lattice_parameters[1];
#endif
#if SUN == 3
    gpu_su_3 gLr3, gLr4, gLr5;
    su3_twist twist;
        twist.phi   = lattice_parameters[1];
        twist.omega = lattice_parameters[2];
#endif

    gid0 = GID;

    if(GID < N2N3N4){
        lattice_gid_to_coords(&gid0,&coord0);
        
        for(i = 0; i < WLN; i++){
#if SUN == 2
            gLr5 = lattice_table_2(lattice_table, &coord0, gid0, X, &twist);
#endif
#if SUN == 3
            gLr5 = lattice_table_3(lattice_table, &coord0, gid0, X, &twist);
#endif
            coord2 = coord0;
            for(j = 1; j < i + 1; j++){
                coord1 = coord2;
                lattice_neighbours_gid(&coord1, &coord2, &gid2, X);
#if SUN == 2
                gLr4 = lattice_table_2(lattice_table, &coord2, gid2, X, &twist);
                gLr3 = matrix_times2(&gLr5, &gLr4);
#endif
#if SUN == 3
                gLr4 = lattice_table_3(lattice_table, &coord2, gid2, X, &twist);
                gLr3 = matrix_times3(&gLr5, &gLr4);
#endif
                gLr5 = gLr3;
            }
#if SUN == 2
            lattice_L_store_2(Lr2, &gLr5, gid0 + i * 1 * ROWSIZE_LR2);
#endif
#if SUN == 3
            lattice_Lr2_store_3(Lr2, &gLr5, gid0 + i * 3 * ROWSIZE_LR2);
#endif
    }
    }
}

                                       __kernel void
lattice_measurement_WLx0(__global hgpu_float4  * Lt,
                       __global hgpu_float4  * Lr1,
                       __global hgpu_float4  * WLx,
                       __global hgpu_float   * lattice_parameters){
#if SUN == 2
    gpu_su_2 wlx1, wlx2, wlx3, wlx4;
#endif
#if SUN == 3
    gpu_su_3 wlx1, wlx2, wlx3, wlx4;
#endif
    coords_4 coord;
    uint gid1;
    uint gid = GID;
    
    if (GID < SITES){
        lattice_gid_to_coords(&gid,&coord);
        coord.t = (coord.t + WLT) % N4;
        lattice_coords_to_gid(&gid1, &coord);

#if SUN == 2
        wlx2 = lattice_L_get_2(Lr1, gid1);
        wlx4 = lattice_L_get_2(Lt, gid);

        wlx1 = matrix_times2(&wlx4, &wlx2);
        wlx4 = matrix_hermitian_gpu_su_2(&wlx1);

        wlx2 = lattice_L_get_2(Lr1, gid);
        wlx3 = matrix_times2(&wlx4, &wlx2);

        lattice_L_store_2(WLx, &wlx3, gid);
#endif
#if SUN == 3
        wlx2 = lattice_L_get_3(Lr1, gid1, ROWSIZE_WLX);
        wlx4 = lattice_Lt_get_3(Lt, gid);

        wlx1 = matrix_times3(&wlx4, &wlx2);
        wlx4 = matrix_hermitian_gpu_su_3(&wlx1);

        wlx2 = lattice_L_get_3(Lr1, gid, ROWSIZE_WLX);
        wlx3 = matrix_times3(&wlx4, &wlx2);

        lattice_L_store_3(WLx, &wlx3, gid, ROWSIZE_WLX);
#endif
    }
}

                                       __kernel void
lattice_measurement_WLx1(__global hgpu_float4  * WLx,
                         __global hgpu_float4  * Lt,
                         __global int          * m){
    unsigned int gid = GID;
    unsigned int gid1;
    coords_4 coord, coord1;
#if SUN == 2
    gpu_su_2 wl, wl1, gLt;
#endif
#if SUN == 3
    gpu_su_3 wl, wl1, gLt;
#endif

    if(GID < SITES)
    {
        lattice_gid_to_coords(&gid, &coord);
#if SUN == 2
        wl1 = lattice_L_get_2(WLx, gid);
#endif
#if SUN == 3
        wl1 = lattice_L_get_3(WLx, gid, ROWSIZE_WLX);
#endif

        coord1 = coord;
        coord1.x +=WLR;
        lattice_coords_to_gid(&gid1, &coord1);

        wl = wl1;
        if(m[coord.x] == K){
#if SUN == 2
            gLt = lattice_L_get_2(Lt, gid1);
            wl = matrix_times2(&wl1, &gLt);
#endif
#if SUN == 3
            gLt = lattice_Lt_get_3(Lt, gid1);
            wl = matrix_times3(&wl1, &gLt);
#endif
        }

#if SUN == 2
        lattice_L_store_2(WLx, &wl, gid);
#endif
#if SUN == 3
        lattice_L_store_3(WLx, &wl, gid, ROWSIZE_WLX);
#endif
    }
}

                                       __kernel void
lattice_measurement_WLx2(__global hgpu_float4  * WLx,
                         __global hgpu_float4  * Lr2,
                         __global int          * m,
                         __global hgpu_float   * lattice_parameters)
{
    int Nx1 = lattice_parameters[3];
    int Nxx = (int)lattice_parameters[4];
    int rowsize = (int)lattice_parameters[5];

    uint gid = GID;
    coords_4 coord, coord0;

    uint gid1, gid2;

    int l1, l;
#if SUN == 2
    gpu_su_2 lr21, lr22, wlx1;
#endif
#if SUN == 3
    gpu_su_3 lr21, lr22, wlx1;
#endif

    if(GID < SITES)
    {
#if SUN == 2
//        lattice_zero2(&wlx1);//DELETE?
        wlx1 = lattice_L_get_2(WLx, gid);
#endif
#if SUN == 3
        //lattice_zero3(&wlx1);//DELETE?
        wlx1 = lattice_L_get_3(WLx, gid, ROWSIZE_WLX);
#endif

        lattice_gid_to_coords(&gid, &coord);
        coord0 = coord;
        coord0.x = 0;

        l1 = WLR - Nxx + coord.x;
        l = (l1 < Nx1) ? l1 - 1 : Nx1 - 1;

        lattice_coords_to_gid(&gid2, &coord0);
        coord0.t = (coord0.t + WLT) % N4;
        lattice_coords_to_gid(&gid1, &coord0);
#if SUN == 2
        gid1 += l * 1 * rowsize;
        gid2 += l * 1 * rowsize;
#endif
#if SUN == 3
        gid1 += l * 3 * rowsize;
        gid2 += l * 3 * rowsize;
#endif

        if(m[coord.x] > K && l >= 0)
        {
#if SUN == 2
            lr22 = lattice_L_get_2(Lr2, gid1);
            lr21 = matrix_hermitian_gpu_su_2(&lr22);

            lr22 = matrix_times2(&lr21, &wlx1);

            lr21 = lattice_L_get_2(Lr2, gid2);

            wlx1 = matrix_times2(&lr22, &lr21);
        }

        lattice_L_store_2(WLx, &wlx1, gid);
#endif
#if SUN == 3
            lr22 = lattice_L_get_3(Lr2, gid1, rowsize);
            lr21 = matrix_hermitian_gpu_su_3(&lr22);

            lr22 = matrix_times3(&lr21, &wlx1);

            lr21 = lattice_L_get_3(Lr2, gid2, rowsize);

            wlx1 = matrix_times3(&lr22, &lr21);
        }

        lattice_L_store_3(WLx, &wlx1, gid, ROWSIZE_WLX);
#endif
    }
}

                                       __kernel void
lattice_measurement_WLx3(__global hgpu_float4  * WLx,
                         __global hgpu_float4  * Lt,
                         __global int          * m,
                         __global hgpu_float   * lattice_parameters)
{
    int KK      = lattice_parameters[3];
    int Nxx     = lattice_parameters[4];
    int rowsize = lattice_parameters[5];

    uint gid = GID;
    coords_4 coord;

    int l1;
    uint gid1;
#if SUN == 2
    gpu_su_2 gLt, wlx1, wlx2;
#endif
#if SUN == 3
    gpu_su_3 gLt, wlx1, wlx2;
#endif

    if(GID < SITES)
    {
        lattice_gid_to_coords(&gid, &coord);
        l1 = WLR - Nxx + coord.x;

#if SUN == 2        
        wlx2 = lattice_L_get_2(WLx, gid);
#endif
#if SUN == 3        
        wlx2 = lattice_L_get_3(WLx, gid, ROWSIZE_WLX);
#endif

        if(m[coord.x] == KK){
            wlx1 = wlx2;

            coord.x = l1;
            lattice_coords_to_gid(&gid1, &coord);
#if SUN == 2
            gLt = lattice_L_get_2(Lt, gid1);

            wlx2 = matrix_times2(&wlx1, &gLt);
        }

        lattice_L_store_2(WLx, &wlx2, gid);
#endif
#if SUN == 3
            gLt = lattice_L_get_3(Lt, gid1, rowsize);

            wlx2 = matrix_times3(&wlx1, &gLt);
        }

        lattice_L_store_3(WLx, &wlx2, gid, ROWSIZE_WLX);
#endif
    }
}

                                       __kernel void
lattice_measurement_WLx4(__global hgpu_float4  * WLx,
                         __global hgpu_double * lattice_measurement0,
                         __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    gpu_su_2 wlX;
    double_su_2 wlx;
#endif
#if SUN == 3
    gpu_su_3 wlX;
    double_su_3 wlx;
#endif
    hgpu_double wl = 0.0;
    int gid;
    gid = GID;

    if(GID < SITES){
#if SUN == 2
        wlX = lattice_L_get_2(WLx, gid);

        wlx = lattice_reconstruct2_double(&wlX);        
        wl = matrix_retrace2_double(&wlx);
#endif
#if SUN == 3
        wlX = lattice_L_get_3(WLx, gid, ROWSIZE_WLX);
        
        wlx = lattice_reconstruct3_double(&wlX);
        wl = matrix_retrace3_double(&wlx);
#endif
    }

    lattice_measurement0[GID] = wl;
}

                                        __kernel void
lattice_measurement_wilson(__global hgpu_float4  * lattice_table,
                           __global hgpu_double  * lattice_measurement0,
                           __global hgpu_double2 * lattice_measurement,
                           __global hgpu_float   * lattice_parameters,
                           __local hgpu_double2  * lattice_lds)
{
#if SUN == 2
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    uint gdi = GID;
    hgpu_double wilson_loop;
    coords_4 coord,  coord2;
    coords_4 coordX, coordX2;
    coords_4 coord_1,coord_2,coord_4;
    uint gdi_1,gdi_2,gdi_4,gdiT;
    gpu_su_2 u1;
    double_su_2 m1, m2, m3, m4, m5, w1;

    su2_twist twist;
        twist.phi   = lattice_parameters[1];
    hgpu_float wilson_R,wilson_T,tmp_R,tmp_T;
        wilson_R    = lattice_parameters[3];
        wilson_T    = lattice_parameters[4];
    wilson_loop = 0.0;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        // _______________________________ bottom link
        lattice_gid_to_coords(&gdi,&coord_1);
        gdi_1 = gdi;
        coord_2 = coord_1;
        u1 = lattice_table_2(lattice_table,&coord_2, gdi,T,&twist);         // [p,T]
        m1 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_2;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_2,&gdiT,T);
            u1 = lattice_table_2(lattice_table,&coord_2,gdiT,T,&twist);     // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m1,&m5);
            m1 = w1;

            tmp_T += 1.0;
        }
        lattice_neighbours_gid(&coord_2,&coord,&gdi_2,T);
        coord_2 = coord;

        // _______________________________ left and right X links
        wilson_loop = lattice_measurement0[gdi];

        // _______________________________ left and right Y links
        // coord_1 - (p,Y)
        // coord_2 - (p+T,Y)
        u1 = lattice_table_2(lattice_table,&coord_1,gdi_1,Y,&twist);        // [p,Y]
        m4 = lattice_reconstruct2_double(&u1);
        u1 = lattice_table_2(lattice_table,&coord_2,gdi_2,Y,&twist);        // [p+T,Y]
        m2 = lattice_reconstruct2_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Y);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Y);

            u1 = lattice_table_2(lattice_table,&coord, gdi, Y,&twist);      // [p,Y]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_2(lattice_table,&coord2,gdiT,Y,&twist);      // [p+T,Y]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Y);

        // _______________________________ top link
        u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Y,T]
        m3 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette2_double(&m1,&m2,&m3,&m4); 

        // _______________________________ left and right Z links
        // coord_1 - (p,Z)
        // coord_2 - (p+T,Z)
        u1 = lattice_table_2(lattice_table,&coord_1,gdi_1,Z,&twist);        // [p,Z]
        m4 = lattice_reconstruct2_double(&u1);
        u1 = lattice_table_2(lattice_table,&coord_2,gdi_2,Z,&twist);        // [p+T,Z]
        m2 = lattice_reconstruct2_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Z);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Z);

            u1 = lattice_table_2(lattice_table,&coord, gdi, Z,&twist);      // [p,Z]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_2(lattice_table,&coord2,gdiT,Z,&twist);      // [p+T,Z]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Z);

        // _______________________________ top link
        u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Z,T]
        m3 = lattice_reconstruct2_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_2(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct2_double(&u1);
            w1 = matrix_times_su2_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette2_double(&m1,&m2,&m3,&m4); 

        out = wilson_loop;
    }
    // first reduction
    reduce_first_step_val_double(lattice_lds,&out, &out2);
    if((TID == 0)) lattice_measurement[BID].x = out2;
#endif

#if SUN == 3
    hgpu_double out  = 0.0;
    hgpu_double out2 = 0.0;
    uint gdi = GID;
    hgpu_double wilson_loop;
    coords_4 coord,  coord2;
    coords_4 coordX, coordX2;
    coords_4 coord_1,coord_2,coord_4;
    uint gdi_1,gdi_2,gdi_4,gdiT;
    gpu_su_3 u1;
    double_su_3 m1, m2, m3, m4, m5, w1;

    su3_twist twist;
        twist.phi   = lattice_parameters[1];
        twist.omega = lattice_parameters[2];
    hgpu_float wilson_R,wilson_T,tmp_R,tmp_T;
        wilson_R    = lattice_parameters[3];
        wilson_T    = lattice_parameters[4];
    wilson_loop = 0.0;

    lattice_lds[TID] = (hgpu_double2) 0.0;
    if (GID<SITES) {
        // _______________________________ bottom link
        lattice_gid_to_coords(&gdi,&coord_1);
        gdi_1 = gdi;
        coord_2 = coord_1;
        u1 = lattice_table_3(lattice_table,&coord_2, gdi,T,&twist);         // [p,T]
        m1 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_2;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_2,&gdiT,T);
            u1 = lattice_table_3(lattice_table,&coord_2,gdiT,T,&twist);     // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m1,&m5);
            m1 = w1;

            tmp_T += 1.0;
        }
        lattice_neighbours_gid(&coord_2,&coord,&gdi_2,T);
        coord_2 = coord;

        // _______________________________ left and right X links
        wilson_loop = lattice_measurement0[gdi];

        // _______________________________ left and right Y links
        u1 = lattice_table_3(lattice_table,&coord_1,gdi_1,Y,&twist);        // [p,Y]
        m4 = lattice_reconstruct3_double(&u1);
        u1 = lattice_table_3(lattice_table,&coord_2,gdi_2,Y,&twist);        // [p+T,Y]
        m2 = lattice_reconstruct3_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Y);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Y);

            u1 = lattice_table_3(lattice_table,&coord, gdi, Y,&twist);      // [p,Y]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_3(lattice_table,&coord2,gdiT,Y,&twist);      // [p+T,Y]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Y);

        // _______________________________ top link
        u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Y,T]
        m3 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette3_double(&m1,&m2,&m3,&m4); 

        // _______________________________ left and right Z links
        u1 = lattice_table_3(lattice_table,&coord_1,gdi_1,Z,&twist);        // [p,Z]
        m4 = lattice_reconstruct3_double(&u1);
        u1 = lattice_table_3(lattice_table,&coord_2,gdi_2,Z,&twist);        // [p+T,Z]
        m2 = lattice_reconstruct3_double(&u1);
        tmp_R  = 1.1;
        coord  = coord_1;
        coord2 = coord_2;

        while (tmp_R<wilson_R){
            coordX  = coord;
            coordX2 = coord2;
            // prepare neighbours
            lattice_neighbours_gid(&coordX,  &coord,  &gdi, Z);
            lattice_neighbours_gid(&coordX2, &coord2, &gdiT,Z);

            u1 = lattice_table_3(lattice_table,&coord, gdi, Z,&twist);      // [p,Z]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m4,&m5);
            m4 = w1;

            u1 = lattice_table_3(lattice_table,&coord2,gdiT,Z,&twist);      // [p+T,Z]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m2,&m5);
            m2 = w1;

            tmp_R += 1.0;
        }
        lattice_neighbours_gid(&coord, &coord_4, &gdi_4,Z);

        // _______________________________ top link
        u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);        // [p+Z,T]
        m3 = lattice_reconstruct3_double(&u1);
        tmp_T = 1.1;

        while (tmp_T<wilson_T){
            coord = coord_4;
            // prepare neighbours
            lattice_neighbours_gid(&coord,&coord_4, &gdi_4,T);
            u1 = lattice_table_3(lattice_table,&coord_4,gdi_4,T,&twist);    // [p,T]
            m5 = lattice_reconstruct3_double(&u1);
            w1 = matrix_times_su3_double(&m3,&m5);
            m3 = w1;

            tmp_T += 1.0;
        }
        wilson_loop += lattice_retrace_plaquette3_double(&m1,&m2,&m3,&m4); 

        out = wilson_loop;
    }
    // first reduction
    reduce_first_step_val_double(lattice_lds,&out, &out2);
    if((TID == 0)) lattice_measurement[BID].x = out2;
#endif
}

                                        __kernel void
reduce_wilson_double2(__global hgpu_double2 * lattice_measurement,
                      __global hgpu_double  * lattice_wilson_loop,
                      __local hgpu_double2  * lattice_lds,
                      uint size,
                      uint index)
{
    reduce_final_step_double2(lattice_lds,lattice_measurement,size);
    if (GID==0) lattice_wilson_loop[index] = lattice_lds[0].x;
}






#endif
