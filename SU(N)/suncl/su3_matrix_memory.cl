/******************************************************************************
 * @file     su3_matrix_memory.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Matrix memory organization for the SU(3) gauge group
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
#ifndef SU3_MATRIX_MEMORY_CL
#define SU3_MATRIX_MEMORY_CL

#include "su3cl.cl"

                    HGPU_INLINE_PREFIX gpu_su_3
lattice_table_3(__global hgpu_float4 * lattice_table,const coords_4 * coord,uint gindex,const uint dir,const su3_twist * twist)
{
    gpu_su_3 m;
    switch (dir){
        case X:
            m.uv1 = lattice_table[gindex +  0 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  4 * ROWSIZE];
            m.uv3 = lattice_table[gindex +  8 * ROWSIZE];
            break;
        case Y:
            m.uv1 = lattice_table[gindex +  1 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  5 * ROWSIZE];
            m.uv3 = lattice_table[gindex +  9 * ROWSIZE];

#ifdef  TBC
//twist here if coord.x = N1 - 1; H = (0, 0, Hz)
#ifdef BIGLAT
    if (((*coord).x + LEFT_SITES/N2N3N4 == 0)||((*coord).x + LEFT_SITES/N2N3N4 == FULL_SITES/N2N3N4)){
#else
    if ((*coord).x == (N1-1)){
#endif
       hgpu_float4 m1,m2,m3,m4,m5,m6;
       hgpu_float phi_p_omega_2 = ((*twist).phi + (*twist).omega)/2;
       hgpu_float omega_m_phi_2 = ((*twist).omega - (*twist).phi)/2;

       hgpu_float4 a1,a2,a3,a4;

       a1.x = (hgpu_float) cos(phi_p_omega_2);
       a1.y = a1.x;
       a1.z = a1.x;
       a1.w = (hgpu_float) cos(omega_m_phi_2);

       a2.x = (hgpu_float) sin(phi_p_omega_2);
       a2.y = a2.x;
       a2.z = a2.x;
       a2.w = (hgpu_float) sin(omega_m_phi_2);

       a3 = (hgpu_float4) cos(omega_m_phi_2);

       a4.x = (hgpu_float) -sin(omega_m_phi_2);
       a4.y =  a4.x;
       a4.z = -a4.x;
       a4.w = -a4.x;

       m1 = m.uv1 * a1;
       m2 = m.uv2 * a2;
       m3 = m.uv2 * a1;
       m4 = m.uv1 * a2;
       m5 = m.uv3 * a3;
       m6 = m.uv3.zwxy * a4;

       m.uv1 = m1 - m2;
       m.uv2 = m3 + m4;
       m.uv3 = m5 + m6;
    }
#endif
            break;
        case Z:
            m.uv1 = lattice_table[gindex +  2 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  6 * ROWSIZE];
            m.uv3 = lattice_table[gindex + 10 * ROWSIZE];
            break;
        case T:
            m.uv1 = lattice_table[gindex +  3 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  7 * ROWSIZE];
            m.uv3 = lattice_table[gindex + 11 * ROWSIZE];
            break;
        default:
            break;
    }

    return m;
}                                                                                                                                                

                    HGPU_INLINE_PREFIX gpu_su_3
lattice_table_notwist_3(__global hgpu_float4 * lattice_table,uint gindex,const uint dir)
{
    gpu_su_3 m;
    switch (dir){
        case X:
            m.uv1 = lattice_table[gindex +  0 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  4 * ROWSIZE];
            m.uv3 = lattice_table[gindex +  8 * ROWSIZE];
            break;
        case Y:
            m.uv1 = lattice_table[gindex +  1 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  5 * ROWSIZE];
            m.uv3 = lattice_table[gindex +  9 * ROWSIZE];
            break;
        case Z:
            m.uv1 = lattice_table[gindex +  2 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  6 * ROWSIZE];
            m.uv3 = lattice_table[gindex + 10 * ROWSIZE];
            break;
        case T:
            m.uv1 = lattice_table[gindex +  3 * ROWSIZE];
            m.uv2 = lattice_table[gindex +  7 * ROWSIZE];
            m.uv3 = lattice_table[gindex + 11 * ROWSIZE];
            break;
        default:
            break;
    }

    return m;
}                                                                                                                                                


                    HGPU_INLINE_PREFIX_VOID void
lattice_store_3(__global hgpu_float4 * lattice_table,gpu_su_3* m,uint gindex,const uint dir){
    switch (dir){
        case 0:
            lattice_table[gindex +  0 * ROWSIZE] = (*m).uv1;
            lattice_table[gindex +  4 * ROWSIZE] = (*m).uv2;
            lattice_table[gindex +  8 * ROWSIZE] = (*m).uv3;
            break;
        case 1:
            lattice_table[gindex +  1 * ROWSIZE] = (*m).uv1;
            lattice_table[gindex +  5 * ROWSIZE] = (*m).uv2;
            lattice_table[gindex +  9 * ROWSIZE] = (*m).uv3;
            break;
        case 2:
            lattice_table[gindex +  2 * ROWSIZE] = (*m).uv1;
            lattice_table[gindex +  6 * ROWSIZE] = (*m).uv2;
            lattice_table[gindex + 10 * ROWSIZE] = (*m).uv3;
            break;
        case 3:
            lattice_table[gindex +  3 * ROWSIZE] = (*m).uv1;
            lattice_table[gindex +  7 * ROWSIZE] = (*m).uv2;
            lattice_table[gindex + 11 * ROWSIZE] = (*m).uv3;
            break;
        default:
            break;
    }
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_store_3_rowsize(__global hgpu_float4 * lattice_table, gpu_su_3* m, uint gindex, const uint dir, int rowsize){
    switch (dir){
        case 0:
            lattice_table[gindex + 0 * rowsize] = (*m).uv1;
            lattice_table[gindex + 4 * rowsize] = (*m).uv2;
            lattice_table[gindex + 8 * rowsize] = (*m).uv3;
            break;
        case 1:
            lattice_table[gindex + 1 * rowsize] = (*m).uv1;
            lattice_table[gindex + 5 * rowsize] = (*m).uv2;
            lattice_table[gindex + 9 * rowsize] = (*m).uv3;
            break;
        case 2:
            lattice_table[gindex + 2 * rowsize] = (*m).uv1;
            lattice_table[gindex + 6 * rowsize] = (*m).uv2;
            lattice_table[gindex + 10 * rowsize] = (*m).uv3;
            break;
        case 3:
            lattice_table[gindex + 3 * rowsize] = (*m).uv1;
            lattice_table[gindex + 7 * rowsize] = (*m).uv2;
            lattice_table[gindex + 11 * rowsize] = (*m).uv3;
            break;
        default:
            break;
    }
}




#endif
