/******************************************************************************
 * @file     su3cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines general procedures for lattice update (SU(2) gauge theory)
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
#ifndef SU2CL_CL
#define SU2CL_CL

#include "sun_common.cl"

                    __attribute__((always_inline)) __private su_2
matrix_times_su2(su_2* u,su_2* v)
{
    su_2 tmp;

    tmp.u1.re = -(*u).u1.im * (*v).u1.im + (*u).u1.re * (*v).u1.re - (*u).u2.im * (*v).v1.im + (*u).u2.re * (*v).v1.re;
    tmp.u1.im =  (*u).u1.re * (*v).u1.im + (*u).u1.im * (*v).u1.re + (*u).u2.re * (*v).v1.im + (*u).u2.im * (*v).v1.re;

    tmp.u2.re = -(*u).u1.im * (*v).u2.im + (*u).u1.re * (*v).u2.re - (*u).u2.im * (*v).v2.im + (*u).u2.re * (*v).v2.re;
    tmp.u2.im =  (*u).u1.re * (*v).u2.im + (*u).u1.im * (*v).u2.re + (*u).u2.re * (*v).v2.im + (*u).u2.im * (*v).v2.re; 

    tmp.v1.re = -(*u).v1.im * (*v).u1.im + (*u).v1.re * (*v).u1.re - (*u).v2.im * (*v).v1.im + (*u).v2.re * (*v).v1.re;
    tmp.v1.im =  (*u).v1.re * (*v).u1.im + (*u).v1.im * (*v).u1.re + (*u).v2.re * (*v).v1.im + (*u).v2.im * (*v).v1.re;

    tmp.v2.re = -(*u).v1.im * (*v).u2.im + (*u).v1.re * (*v).u2.re - (*u).v2.im * (*v).v2.im + (*u).v2.re * (*v).v2.re;
    tmp.v2.im =  (*u).v1.re * (*v).u2.im + (*u).v1.im * (*v).u2.re + (*u).v2.re * (*v).v2.im + (*u).v2.im * (*v).v2.re;

    return tmp;
}

                     __attribute__((always_inline)) __private su_2
lattice_reconstruct2(gpu_su_2* m)
{
    su_2 result;

    result.u1.re = (*m).uv1.x;
    result.u1.im = (*m).uv1.z;
    result.u2.re = (*m).uv1.y;
    result.u2.im = (*m).uv1.w;

    result.v1.re = -result.u2.re;
    result.v1.im =  result.u2.im;
    result.v2.re =  result.u1.re;
    result.v2.im = -result.u1.im;

    return result;
}

#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                   
                                                                                                                                                                   
