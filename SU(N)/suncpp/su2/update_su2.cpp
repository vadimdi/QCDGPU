/******************************************************************************
 * @file     update_su2.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Implementation of SU(2) heat-bath algorithm
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

#include "update_su2.h"

void update_link(su_2 *a, su_2 stap, hgpu_double beta, int nhit, int gid, int gid_start, int fsites, int final, PRNG_CL::PRNG *prngCPU){
    su_2 aH, c;
    bool flag = false;
    float rnd[4];
    //double rnd[4];
    
    hgpu_double Mx, My, Mz, Mw;
    hgpu_double det, bdet, cosrnd;
    hgpu_double delta = 0.0;
    hgpu_double costh,sinth,cosal,sinal,phi,sinphi,cosphi;
    int i = 0;
    
    Mx = (stap.u1.re + stap.v2.re);
    My = (stap.u1.im - stap.v2.im);
    Mz = (stap.u2.re - stap.v1.re);
    Mw = (stap.u2.im + stap.v1.im);
    
    det = sqrt(Mx * Mx + My * My + Mz * Mz + Mw * Mw);
   
    bdet = beta * det;

    aH.u1.re =  Mx / det;
    aH.u1.im = -My / det;
    aH.u2.re = -Mz / det;
    aH.u2.im = -Mw / det;
    aH.v2.re =  aH.u1.re;
    aH.v2.im = -aH.u1.im;
    aH.v1.re = -aH.u2.re;
    aH.v1.im =  aH.u2.im;
    
    while ((i < nhit) && (flag == false)){
        if (gid_start){
            rnd[0] = fabs(sin((0.005*(1.0 + nhit)+270.0/fsites)*gid));
            rnd[1] = fabs(cos((0.005*(1.0 + nhit)+ 60.0/fsites)*gid));
            rnd[2] = fabs(sin((0.005*(1.0 + nhit)-150.0/fsites)*gid));
            rnd[3] = fabs(cos((0.005*(1.0 + nhit)-380.0/fsites)*gid));
        } else {
            /*rnd[0] = (float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd[1] = (float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd[2] = (float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            rnd[3] = (float) (((double) rand()) / ((double) RAND_MAX + 1.0));*/
            prngCPU->produce_CPU(rnd, 4);
        }
        cosrnd = cos(PI2 * rnd[1]);
        if(!gid_start){
            delta = -(log(1.0 - rnd[0]) + cosrnd * cosrnd * log(1.0 - rnd[2])) / bdet;
            if ((rnd[3] * rnd[3])<=(1.0-0.5*delta)) flag=true;
        } else {
            delta = cosrnd * cosrnd / bdet;  
            flag = true;
        }
        i++;
    }
    
    if (flag) {
        if (gid_start){
            rnd[0] = fabs(cos((0.08-270.0/fsites)*gid));
            rnd[1] = fabs(sin((0.08-60.0/fsites)*gid));
        } else {
            //rnd[0] = (float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            //rnd[1] = (float) (((double) rand()) / ((double) RAND_MAX + 1.0));
            prngCPU->produce_CPU(rnd, 2);
        }
            cosal = 1.0 - delta;
            costh = 2.0 * rnd[0] - 1.0;
            sinth = sqrt(1.0 - costh * costh);
            sinal = sqrt(1.0 - cosal * cosal);
            phi   = PI2 * rnd[1];

            sinphi = sin(phi);
            cosphi = cos(phi);
            
            c.u1.re =  cosal;
            c.u1.im =  sinal * costh;
            c.u2.re =  sinal * sinth * sinphi;
            c.u2.im =  sinal * sinth * cosphi;
            c.v2.re =  c.u1.re;
            c.v2.im = -c.u1.im;
            c.v1.re = -c.u2.re;
            c.v1.im =  c.u2.im;
            
            *a = c * aH;
    } 
    
    if ((!final) && (!flag))
        lattice_unity(a);
}
