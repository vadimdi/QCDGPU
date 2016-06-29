/******************************************************************************
 * @file     update_su3.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Implementation of SU(3) heat-bath algorithm
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

#include "update_su3.h"

void update_link(su_3 *U, su_3 stap, hgpu_double beta, int nhit, int gid, int gid_start, int fsites, int final, PRNG_CL::PRNG *prngCPU){
    su_3 U0 = *U;
    su_3 x0, Vg;
    su_2 r0;
    
    x0 = U0 * stap;
    
        r0.u1.re = x0.u1.re;                                                         //   /  u1  u2  0  \   //
        r0.u1.im = x0.u1.im;                                                       //   |   v1  v2   0  |   //
        r0.u2.re = x0.u2.re;                                                         //   \  0    0    1  /   //
        r0.u2.im = x0.u2.im;                                                       //                            //
        r0.v1.re = x0.v1.re;                                                         /////////////////////
        r0.v1.im = x0.v1.im;
        r0.v2.re = x0.v2.re;
        r0.v2.im = x0.v2.im;
    update_link(&r0, r0, beta, nhit, gid, gid_start, fsites, 0, prngCPU);
    
    Vg.u1.re = r0.u1.re;            Vg.u2.re = r0.u2.re;            Vg.u3.re = 0.0;
    Vg.u1.im = r0.u1.im;          Vg.u2.im = r0.u2.im;           Vg.u3.im = 0.0;
    
    Vg.v1.re = r0.v1.re;            Vg.v2.re = r0.v2.re;            Vg.v3.re = 0.0;
    Vg.v1.im = r0.v1.im;          Vg.v2.im = r0.v2.im;           Vg.v3.im = 0.0;
    
    Vg.w1.re = 0.0;                  Vg.w2.re = 0.0;                   Vg.w3.re = 1.0;
    Vg.w1.im = 0.0;                 Vg.w2.im = 0.0;                  Vg.w3.im = 0.0;
    
    U0 = Vg * U0;
    x0 = U0 * stap;
    
        r0.u1.re = x0.u1.re;                                                         //   /  u1  0  u2  \   //
        r0.u1.im = x0.u1.im;                                                       //   |  0     1  0     |   //
        r0.u2.re = x0.u3.re;                                                         //   \  v1  0  v2  /   //
        r0.u2.im = x0.u3.im;                                                       //                            //
        r0.v1.re = x0.w1.re;                                                        /////////////////////
        r0.v1.im = x0.w1.im;
        r0.v2.re = x0.w3.re;
        r0.v2.im = x0.w3.im;
    update_link(&r0, r0, beta, nhit, gid, gid_start, fsites, 0, prngCPU);
    
    Vg.u1.re = r0.u1.re;            Vg.u2.re = 0.0;                    Vg.u3.re = r0.u2.re;
    Vg.u1.im = r0.u1.im;          Vg.u2.im = 0.0;                   Vg.u3.im = r0.u2.im;
    
    Vg.v1.re = 0.0;                   Vg.v2.re = 1.0;                   Vg.v3.re = 0.0;
    Vg.v1.im = 0.0;                  Vg.v2.im = 0.0;                  Vg.v3.im = 0.0;
    
    Vg.w1.re = r0.v1.re;           Vg.w2.re = 0.0;                   Vg.w3.re = r0.v2.re;
    Vg.w1.im = r0.v1.im;         Vg.w2.im = 0.0;                  Vg.w3.im = r0.v2.im;
    
    U0 = Vg * U0;
    x0 = U0 * stap;
    
        r0.u1.re = x0.v2.re;                                                         //   / 1   0   0    \   //
        r0.u1.im = x0.v2.im;                                                       //   |  0   u1  u2  |   //
        r0.u2.re = x0.v3.re;                                                         //   \ 0   v1  v2  /   //
        r0.u2.im = x0.v3.im;                                                       //                            //
        r0.v1.re = x0.w2.re;                                                        /////////////////////
        r0.v1.im = x0.w2.im;
        r0.v2.re = x0.w3.re;
        r0.v2.im = x0.w3.im;
    update_link(&r0, r0, beta, nhit, gid, gid_start, fsites, 0, prngCPU);
    
    Vg.u1.re = 1.0;                   Vg.u2.re = 0.0;                    Vg.u3.re = 0.0;
    Vg.u1.im = 0.0;                  Vg.u2.im = 0.0;                   Vg.u3.im = 0.0;
    
    Vg.v1.re = 0.0;                   Vg.v2.re = r0.u1.re;            Vg.v3.re = r0.u2.re;
    Vg.v1.im = 0.0;                  Vg.v2.im = r0.u1.im;           Vg.v3.im = r0.u2.im;
    
    Vg.w1.re = 0.0       ;           Vg.w2.re = r0.v1.re;            Vg.w3.re = r0.v2.re;
    Vg.w1.im = 0.0;                 Vg.w2.im = r0.v1.im;           Vg.w3.im = r0.v2.im;
    
    U0 = Vg * U0;
    GramSchmidt(&U0);
    
    *U = U0;
}
