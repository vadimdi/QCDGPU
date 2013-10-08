/******************************************************************************
 * @file     boundary.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Definition of boundaries of sublattices in divided lattice
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013, Vadim Demchik, Natalia Kolomoyets
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

#ifndef BOUNDARY_CL
#define BOUNDARY_CL

#include "complex.h"
#include "model.cl"
#include "su3_matrix_memory.cl"

                                        __kernel void
lattice_get_boundary_low(__global hgpu_float4 * lattice_table,
                         __global hgpu_float4 * lattice_boundary)
{
#if SUN == 3
    gpu_su_3 matrix;
    uint gindex = GID;

    if (GID<N2N3N4) {
        matrix = lattice_table_notwist_3(lattice_table,gindex,X);
        lattice_store_3(lattice_boundary,&matrix,GID,X);
        matrix = lattice_table_notwist_3(lattice_table,gindex,Y);
        lattice_store_3(lattice_boundary,&matrix,GID,Y);
        matrix = lattice_table_notwist_3(lattice_table,gindex,Z);
        lattice_store_3(lattice_boundary,&matrix,GID,Z);
        matrix = lattice_table_notwist_3(lattice_table,gindex,T);
        lattice_store_3(lattice_boundary,&matrix,GID,T);
    }

#endif
}

                                        __kernel void
lattice_get_boundary_high(__global hgpu_float4 * lattice_table,
                          __global hgpu_float4 * lattice_boundary)
{
#if SUN == 3
    gpu_su_3 matrix;
    uint gindex = GID + abs(N1-3) * N2N3N4;

    if (GID<N2N3N4) {
        matrix = lattice_table_notwist_3(lattice_table,gindex,X);
        lattice_store_3(lattice_boundary,&matrix,GID,X);
        matrix = lattice_table_notwist_3(lattice_table,gindex,Y);
        lattice_store_3(lattice_boundary,&matrix,GID,Y);
        matrix = lattice_table_notwist_3(lattice_table,gindex,Z);
        lattice_store_3(lattice_boundary,&matrix,GID,Z);
        matrix = lattice_table_notwist_3(lattice_table,gindex,T);
        lattice_store_3(lattice_boundary,&matrix,GID,T);
    }

#endif
}

                                        __kernel void
lattice_put_boundary_low(__global hgpu_float4 * lattice_table,
                         __global hgpu_float4 * lattice_boundary)
{
#if SUN == 3
    gpu_su_3 matrix;
    uint gindex = GID;

    if (GID<N2N3N4) {
        matrix = lattice_table_notwist_3(lattice_boundary,GID,X);
        lattice_store_3(lattice_table,&matrix,gindex,X);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Y);
        lattice_store_3(lattice_table,&matrix,gindex,Y);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Z);
        lattice_store_3(lattice_table,&matrix,gindex,Z);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,T);
        lattice_store_3(lattice_table,&matrix,gindex,T);
    }

#endif
}

                                        __kernel void
lattice_put_boundary_high(__global hgpu_float4 * lattice_table,
                          __global hgpu_float4 * lattice_boundary)
{
#if SUN == 3
    gpu_su_3 matrix;
    uint gindex = GID + abs(N1-3) * N2N3N4;

    if (GID<N2N3N4) {
        matrix = lattice_table_notwist_3(lattice_boundary,GID,X);
        lattice_store_3(lattice_table,&matrix,gindex,X);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Y);
        lattice_store_3(lattice_table,&matrix,gindex,Y);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Z);
        lattice_store_3(lattice_table,&matrix,gindex,Z);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,T);
        lattice_store_3(lattice_table,&matrix,gindex,T);
    }

#endif
}

                                        __kernel void
lattice_put_boundary_next(__global hgpu_float4 * lattice_table,
                          __global hgpu_float4 * lattice_boundary)
{
#if SUN == 3
    gpu_su_3 matrix;
    uint gindex = GID + abs(N1-2) * N2N3N4;

    if (GID<N2N3N4) {
        matrix = lattice_table_notwist_3(lattice_boundary,GID,X);
        lattice_store_3(lattice_table,&matrix,gindex,X);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Y);
        lattice_store_3(lattice_table,&matrix,gindex,Y);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Z);
        lattice_store_3(lattice_table,&matrix,gindex,Z);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,T);
        lattice_store_3(lattice_table,&matrix,gindex,T);
    }
#endif
}

                                        __kernel void
lattice_put_boundary_previous(__global hgpu_float4 * lattice_table,
                              __global hgpu_float4 * lattice_boundary)
{
#if SUN == 3
    gpu_su_3 matrix;
    uint gindex = GID + abs(N1-1) * N2N3N4;

    if (GID<N2N3N4) {
        matrix = lattice_table_notwist_3(lattice_boundary,GID,X);
        lattice_store_3(lattice_table,&matrix,gindex,X);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Y);
        lattice_store_3(lattice_table,&matrix,gindex,Y);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,Z);
        lattice_store_3(lattice_table,&matrix,gindex,Z);
        matrix = lattice_table_notwist_3(lattice_boundary,GID,T);
        lattice_store_3(lattice_table,&matrix,gindex,T);
    }
#endif
}



#endif
