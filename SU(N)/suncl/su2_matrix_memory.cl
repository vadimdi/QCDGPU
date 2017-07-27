/******************************************************************************
 * @file     su2_matrix_memory.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Matrix memory organization for the SU(2) gauge group
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2016 Vadim Demchik, Natalia Kolomoyets
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
matrix_times2(gpu_su_2* u,gpu_su_2* v)
{
    gpu_su_2 a;

        a.uv1.x = -(*u).uv1.w * (*v).uv1.w + (*u).uv1.x * (*v).uv1.x - (*u).uv1.y * (*v).uv1.y - (*u).uv1.z * (*v).uv1.z;
        a.uv1.z =  (*u).uv1.y * (*v).uv1.w + (*u).uv1.z * (*v).uv1.x - (*u).uv1.w * (*v).uv1.y + (*u).uv1.x * (*v).uv1.z;

        a.uv1.y = -(*u).uv1.z * (*v).uv1.w + (*u).uv1.y * (*v).uv1.x + (*u).uv1.x * (*v).uv1.y + (*u).uv1.w * (*v).uv1.z;
        a.uv1.w =  (*u).uv1.x * (*v).uv1.w + (*u).uv1.w * (*v).uv1.x + (*u).uv1.z * (*v).uv1.y - (*u).uv1.y * (*v).uv1.z;

    return a;
}

                    __attribute__((always_inline)) __private gpu_su_2
lattice_table_2(__global hgpu_float4 * lattice_table,const coords_4 * coord,uint gindex,const uint dir,const su2_twist * twist)
{
    gpu_su_2 m;
    switch (dir){
        case X: m.uv1 = lattice_table[gindex +  0 * ROWSIZE]; break;
        case Y: m.uv1 = lattice_table[gindex +  1 * ROWSIZE];

#ifdef  TBC
//twist here if coord.x=N1; H = (0, 0, Hz)
    	gpu_su_2 Omega;	
	gpu_su_2 tmp;
	
#ifdef BIGLAT
    if (((*coord).x + LEFT_SITES/N2N3N4 == 0)||((*coord).x + LEFT_SITES/N2N3N4 == FULL_SITES/N2N3N4)){
#else
    if ((*coord).x == (N1-1)){
#endif

	hgpu_float cosphi, sinphi;
	sinphi = (hgpu_float) sin((*twist).phi/2);
	cosphi = (hgpu_float) cos((*twist).phi/2);
		
	Omega.uv1.x = cosphi;
	Omega.uv1.z = sinphi;
	
	Omega.uv1.y = 0.0;
	Omega.uv1.w = 0.0;

	tmp = matrix_times2(&Omega, &m);

	m.uv1 = tmp.uv1;
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
