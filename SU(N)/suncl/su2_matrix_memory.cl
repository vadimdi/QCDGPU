/******************************************************************************
 * @file     su2_matrix_memory.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Matrix memory organization for the SU(2) gauge group
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
#ifndef SU2_MATRIX_MEMORY_CL
#define SU2_MATRIX_MEMORY_CL

                    __attribute__((always_inline)) __private gpu_su_2
lattice_table_2(__global hgpu_float4 * lattice_table,const coords_4 * coord,uint gindex,const uint dir,const su2_twist * twist)
{
    gpu_su_2 m;
    switch (dir){
        case X: m.uv1 = lattice_table[gindex +  0 * ROWSIZE]; break;
        case Y: m.uv1 = lattice_table[gindex +  1 * ROWSIZE];

#ifdef  TBC
//twist here if coord.z=N3
    if ((*coord).z == (N3-1)){

        hgpu_float tmp_x;
        hgpu_float tmp_z;

        hgpu_float cosphi, sinphi;
        sinphi = (hgpu_float) sin((*twist).phi);
        cosphi = (hgpu_float) cos((*twist).phi);

        tmp_x = m.uv1.x * cosphi - m.uv1.z * sinphi;
        tmp_z = m.uv1.z * cosphi + m.uv1.x * sinphi;

        m.uv1.x = tmp_x;
        m.uv1.z = tmp_z;
    }
#endif
            break;
        case Z: m.uv1 = lattice_table[gindex +  2 * ROWSIZE]; break;
        case T: m.uv1 = lattice_table[gindex +  3 * ROWSIZE]; break;
        default: break;
    }

    return m;
}                                                                                                                                                

                    __attribute__((always_inline)) __private gpu_su_2
lattice_table_notwist_2(__global hgpu_float4 * lattice_table,uint gindex,const uint dir)
{
    gpu_su_2 m;
    switch (dir){
        case X:
            m.uv1 = lattice_table[gindex +  0 * ROWSIZE];
            break;
        case Y:
            m.uv1 = lattice_table[gindex +  1 * ROWSIZE];
            break;
        case Z:
            m.uv1 = lattice_table[gindex +  2 * ROWSIZE];
            break;
        case T:
            m.uv1 = lattice_table[gindex +  3 * ROWSIZE];
            break;
        default:
            break;
    }

    return m;
}                                                                                                                                                


                    __attribute__((always_inline)) void
lattice_store_2(__global hgpu_float4 * lattice_table,gpu_su_2* m,uint gindex,const uint dir){
    switch (dir){
        case 0:
            lattice_table[gindex +  0 * ROWSIZE] = (*m).uv1;
            break;
        case 1:
            lattice_table[gindex +  1 * ROWSIZE] = (*m).uv1;
            break;
        case 2:
            lattice_table[gindex +  2 * ROWSIZE] = (*m).uv1;
            break;
        case 3:
            lattice_table[gindex +  3 * ROWSIZE] = (*m).uv1;
            break;
        default:
            break;
    }
}  




#endif
