/******************************************************************************
 * @file     su2_update_cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Contains functions for lattice update (SU(2) gauge theory)
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013, 2014 Vadim Demchik, Natalia Kolomoyets
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
#ifndef SU2UPDATECL_CL
#define SU2UPDATECL_CL

// _________________ single precision ________________________
                     __attribute__((always_inline)) void
lattice_unity_2(su_2* matrix)
{
    (*matrix).u1.re = 1.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 1.0;
    (*matrix).v2.im = 0.0;
}

                     __attribute__((always_inline)) void
lattice_zero_2(su_2* matrix)
{
    (*matrix).u1.re = 0.0;
    (*matrix).u1.im = 0.0;

    (*matrix).u2.re = 0.0;
    (*matrix).u2.im = 0.0;

    (*matrix).v1.re = 0.0;
    (*matrix).v1.im = 0.0;

    (*matrix).v2.re = 0.0;
    (*matrix).v2.im = 0.0;
}

                    __attribute__((always_inline)) void
lattice_unity2(gpu_su_2* matrix)
{
   (*matrix).uv1 = (hgpu_float4) (1.0, 0.0, 0.0, 0.0);
}

                    __attribute__((always_inline)) __private gpu_su_2
matrix_hermitian2(gpu_su_2* u)
{
    __private gpu_su_2 tmp;
    __private hgpu_float4 u1;
    
    u1 = (*u).uv1;
    tmp.uv1 = (hgpu_float4)(u1.x, -u1.y, -u1.z, -u1.w);

    return tmp;
}

                    __attribute__((always_inline)) __private su_2
matrix_add2(su_2* u,su_2* v)
{
    su_2 tmp;

    tmp.u1.re = (*u).u1.re + (*v).u1.re;
    tmp.u1.im = (*u).u1.im + (*v).u1.im;
       tmp.u2.re = (*u).u2.re + (*v).u2.re;
       tmp.u2.im = (*u).u2.im + (*v).u2.im;

    tmp.v1.re = (*u).v1.re + (*v).v1.re;
    tmp.v1.im = (*u).v1.im + (*v).v1.im;
       tmp.v2.re = (*u).v2.re + (*v).v2.re;
       tmp.v2.im = (*u).v2.im + (*v).v2.im;

    return tmp;
}

                    __attribute__((always_inline)) void
lattice_random2(gpu_su_2* matrix,__global const hgpu_prng_float4 * prns,uint gidprn)
{
    /* M. Di Pierro */
    gpu_su_2 m1;
    __private hgpu_float alpha,phi;
    __private hgpu_float sinth,costh;
    __private hgpu_float sinph,cosph;
    __private hgpu_float sinal;
    __private hgpu_float a0,a1,a2,a3;
    __private hgpu_float t1;

    __private hgpu_float4 rnd;
    rnd.x = (hgpu_float) prns[gidprn].x;
    rnd.y = (hgpu_float) prns[gidprn].y;
    rnd.z = (hgpu_float) prns[gidprn].z;
    rnd.w = (hgpu_float) prns[gidprn].w;

    lattice_unity2(&m1);

        alpha    = PI  * rnd.x;
        phi      = PI2 * rnd.y;
        costh    = 2.0f * rnd.z - 1.0f;

        sinal   = sin(alpha);
        a0      = cos(alpha);

        sinph   = sin(phi);
        cosph   = cos(phi);

        sinth     = sqrt(1.0f-costh*costh);
        t1        = sinal * sinth;
        a1        = t1 * cosph;
        a2        = t1 * sinph;
        a3        = sinal * costh;

     (*matrix).uv1 = (hgpu_float4)(a0, a2, a3, a1);
}

                    __attribute__((always_inline)) void
lattice_su2_Normalize(gpu_su_2* matrix)
{
    hgpu_float det;

    det = sqrt((*matrix).uv1.x * (*matrix).uv1.x + (*matrix).uv1.y * (*matrix).uv1.y + (*matrix).uv1.z * (*matrix).uv1.z + (*matrix).uv1.w * (*matrix).uv1.w);

    (*matrix).uv1 /= det;
}

                    __attribute__((always_inline)) __private su_2
lattice_staple_hermitian2(gpu_su_2* u1, gpu_su_2* u2, gpu_su_2* u3)
{                                             
    gpu_su_2 m1, m2, m3;                      ///////////////////////////
    su_2 result;                              //                       //
                                              //           u2          //
    m1 = matrix_hermitian2(u2);               //      _____/_____      //
    m2 = matrix_times2(u1,&m1);               //     |     \     |     //
    m1 = matrix_hermitian2(u3);               //     |           |     //
    m3 = matrix_times2(&m2,&m1);              // u3 \|/         /|\ u1 //
                                              //     |           |     //
    result = lattice_reconstruct2(&m3);       //     |           |     //
                                              //                       //
    return result;                            ///////////////////////////
}

                    __attribute__((always_inline)) __private su_2
lattice_staple_backward2(gpu_su_2* u1, gpu_su_2* u2, gpu_su_2* u3)
{                                             ///////////////////////////
    gpu_su_2 m1, m2, m3;                      //                       //
    su_2 result;                              //           u2          //
                                              //      _____\_____      //
    m1 = matrix_hermitian2(u1);               //     |     /     |     //
    m2 = matrix_times2(&m1,u2);               //     |           |     //
    m3 = matrix_times2(&m2,u3);               // u3 /|\         \|/ u1 //
                                              //     |           |     //
    result = lattice_reconstruct2(&m3);       //     |           |     //
                                              //                       //
    return result;                            ///////////////////////////
}

                    __attribute__((always_inline)) __private su_2
lattice_staple_hermitian_backward2(gpu_su_2* u1, gpu_su_2* u2, gpu_su_2* u3)
{
    gpu_su_2 m1, m2, m3;                      ///////////////////////////
    su_2 result;                              //                       //
                                              //     |           |     //
    m1 = matrix_hermitian2(u1);               //     |           |     //
    m2 = matrix_hermitian2(u2);               // u3 /|\         \|/ u1 //
    m3 = matrix_times2(&m1,&m2);              //     |           |     //
    m1 = matrix_times2(&m3,u3);               //     |_____/_____|     //
                                              //           \           //
    result = lattice_reconstruct2(&m1);       //           u2          //
                                              //                       //
    return result;                            ///////////////////////////
}

                    __attribute__((always_inline)) __private su_2
lattice_staple_2(__global hgpu_float4 * lattice_table, uint gindex,const uint dir,const su2_twist * twist)
{
        coords_4 coord,coordX,coordY,coordZ,coordT;
        coords_4 coord10,coord11,coord12,coord13,coord14,coord15;

        uint gdiX,   gdiY,   gdiZ,   gdiT;
        uint gdiXm,  gdiYm,  gdiZm,  gdiTm;

            uint gdiYmX, gdiZmX, gdiTmX;
            uint gdiXmY, gdiZmY, gdiTmY;
            uint gdiYmZ, gdiXmZ, gdiTmZ;
            uint gdiXmT, gdiYmT, gdiZmT;


        gpu_su_2 m1,m2,m3;
        su_2 staple, staple1, staple2;

        lattice_gid_to_coords(&gindex,&coord);

        // prepare neighbours
        lattice_neighbours_gid(&coord,&coordX,&gdiX,X);
        lattice_neighbours_gid(&coord,&coordY,&gdiY,Y);
        lattice_neighbours_gid(&coord,&coordZ,&gdiZ,Z);
        lattice_neighbours_gid(&coord,&coordT,&gdiT,T);

    switch (dir){
        case X:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiYm,Y);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiZm,Z);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiTm,T);

            lattice_neighbours_gid(&coord10,&coord13,&gdiYmX,X);
            lattice_neighbours_gid(&coord11,&coord14,&gdiZmX,X);
            lattice_neighbours_gid(&coord12,&coord15,&gdiTmX,X);

                 m1 = lattice_table_2(lattice_table,&coordX,gdiX,Y,twist);    // [p+X,Y]
                 m2 = lattice_table_2(lattice_table,&coordY,gdiY,X,twist);    // [p+Y,X]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,Y,twist);  // [p,Y]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  x-y: [p+X,Y] -[p+Y,X]*-[p,Y]*

                 m1 = lattice_table_2(lattice_table,&coord13,gdiYmX,Y,twist); // [p-Y+X,Y]
                 m2 = lattice_table_2(lattice_table,&coord10,gdiYm,X,twist);  // [p-Y,X]
                 m3 = lattice_table_2(lattice_table,&coord10,gdiYm,Y,twist);  // [p-Y,Y]
            staple2 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -x-y: [p-Y+X,Y]*-[p-Y,X]*-[p-Y,Y]

            staple = matrix_add2(&staple1,&staple2); // [x+y]+[x-y]

                 m1 = lattice_table_2(lattice_table,&coordX,gdiX,Z,twist);    // [p+X,Z]
                 m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,twist);    // [p+Z,X]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,Z,twist);  // [p,Z]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  x-z: [p+X,Z] -[p+Z,X]*-[p,Z]*

            staple2 = matrix_add2(&staple,&staple1); // [x+y]+[x-y]+[x+z]

                 m1 = lattice_table_2(lattice_table,&coord14,gdiZmX,Z,twist); // [p-Z+X,Z]
                 m2 = lattice_table_2(lattice_table,&coord11,gdiZm,X,twist);  // [p-Z,X]
                 m3 = lattice_table_2(lattice_table,&coord11,gdiZm,Z,twist);  // [p-Z,Z]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -x-z: [p-Z+X,Z]*-[p-Z,X]*-[p-Z,Z]

            staple = matrix_add2(&staple1,&staple2); // [x+y]+[x-y]+[x+z]+[x-z]

                 m1 = lattice_table_2(lattice_table,&coordX,gdiX,T,twist);    // [p+X,T]
                 m2 = lattice_table_2(lattice_table,&coordT,gdiT,X,twist);    // [p+T,X]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,T,twist);  // [p,T]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  x-t: [p+X,T] -[p+T,X]*-[p,T]*

            staple2 = matrix_add2(&staple,&staple1); // [x+y]+[x-y]+[x+z]+[x-z]+[x+t]

                 m1 = lattice_table_2(lattice_table,&coord15,gdiTmX,T,twist); // [p-T+X,T]
                 m2 = lattice_table_2(lattice_table,&coord12,gdiTm,X,twist);  // [p-T,X]
                 m3 = lattice_table_2(lattice_table,&coord12,gdiTm,T,twist);  // [p-T,T]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -x-t: [p-T+X,T]*-[p-T,X]*-[p-T,T]
            
            break;
        case Y:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiXm,X);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiZm,Z);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiTm,T);

            lattice_neighbours_gid(&coord10,&coord13,&gdiXmY,Y);
            lattice_neighbours_gid(&coord11,&coord14,&gdiZmY,Y);
            lattice_neighbours_gid(&coord12,&coord15,&gdiTmY,Y);
   
                 m1 = lattice_table_2(lattice_table,&coordY,gdiY,X,twist);    // [p+Y,X]
                 m2 = lattice_table_2(lattice_table,&coordX,gdiX,Y,twist);    // [p+X,Y]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,X,twist);  // [p,X]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  y-x: [p+Y,X] -[p+X,Y]*-[p,X]*

                 m1 = lattice_table_2(lattice_table,&coord13,gdiXmY,X,twist); // [p-X+Y,X]
                 m2 = lattice_table_2(lattice_table,&coord10,gdiXm,Y,twist);  // [p-X,Y]
                 m3 = lattice_table_2(lattice_table,&coord10,gdiXm,X,twist);  // [p-X,X]
            staple2 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -y-x: [p-X+Y,X]*-[p-X,Y]*-[p-X,X]

            staple = matrix_add2(&staple1,&staple2); // [y+x]+[y-x]

                 m1 = lattice_table_2(lattice_table,&coordY,gdiY,Z,twist);    // [p+Y,Z]
                 m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,twist);    // [p+Z,Y]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,Z,twist);  // [p,Z]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  y-z: [p+Y,Z] -[p+Z,Y]*-[p,Z]*

            staple2 = matrix_add2(&staple,&staple1); // [y+x]+[y-x]+[y+z]

                 m1 = lattice_table_2(lattice_table,&coord14,gdiZmY,Z,twist); // [p-Z+Y,Z]
                 m2 = lattice_table_2(lattice_table,&coord11,gdiZm,Y,twist);  // [p-Z,Y]
                 m3 = lattice_table_2(lattice_table,&coord11,gdiZm,Z,twist);  // [p-Z,Z]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -y-z: [p-Z+Y,Z]*-[p-Z,Y]*-[p-Z,Z] 

            staple = matrix_add2(&staple1,&staple2); // [y+x]+[y-x]+[y+z]+[y-z]

                 m1 = lattice_table_2(lattice_table,&coordY,gdiY,T,twist);    // [p+Y,T]
                 m2 = lattice_table_2(lattice_table,&coordT,gdiT,Y,twist);    // [p+T,Y]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,T,twist);  // [p,T]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  y-t: [p+Y,T] -[p+T,Y]*-[p,T]*   

            staple2 = matrix_add2(&staple,&staple1); // [y+x]+[y-x]+[y+z]+[y-z]+[y+t]

                 m1 = lattice_table_2(lattice_table,&coord15,gdiTmY,T,twist); // [p-T+Y,T]
                 m2 = lattice_table_2(lattice_table,&coord12,gdiTm,Y,twist);  // [p-T,Y]
                 m3 = lattice_table_2(lattice_table,&coord12,gdiTm,T,twist);  // [p-T,T]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -y-t: [p-T+Y,T]*-[p-T,Y]*-[p-T,T]   

            break;
        case Z:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiXm,X);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiYm,Y);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiTm,T);

            lattice_neighbours_gid(&coord10,&coord13,&gdiXmZ,Z);
            lattice_neighbours_gid(&coord11,&coord14,&gdiYmZ,Z);
            lattice_neighbours_gid(&coord12,&coord15,&gdiTmZ,Z);
    
                 m1 = lattice_table_2(lattice_table,&coordZ,gdiZ,Y,twist);    // [p+Z,Y]
                 m2 = lattice_table_2(lattice_table,&coordY,gdiY,Z,twist);    // [p+Y,Z]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,Y,twist);  // [p,Y]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  z-y: [p+Z,Y] -[p+Y,Z]*-[p,Y]*

                 m1 = lattice_table_2(lattice_table,&coord14,gdiYmZ,Y,twist); // [p-Y+Z,Y]
                 m2 = lattice_table_2(lattice_table,&coord11,gdiYm,Z,twist);  // [p+Y,Z]
                 m3 = lattice_table_2(lattice_table,&coord11,gdiYm,Y,twist);  // [p-Y,Y]
            staple2 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -z-y: [p-Y+Z,Y]*-[p-Y,Z]*-[p-Y,Y]

            staple = matrix_add2(&staple1,&staple2); // [z+y]+[z-y]

                 m1 = lattice_table_2(lattice_table,&coordZ,gdiZ,X,twist);    // [p+Z,X]
                 m2 = lattice_table_2(lattice_table,&coordX,gdiX,Z,twist);    // [p+X,Z]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,X,twist);  // [p,X]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  z-x: [p+Z,X] -[p+X,Z]*-[p,X]*

            staple2 = matrix_add2(&staple,&staple1); // [z+y]+[z-y]+[z+x]

                 m1 = lattice_table_2(lattice_table,&coord13,gdiXmZ,X,twist); // [p-X+Z,X]
                 m2 = lattice_table_2(lattice_table,&coord10,gdiXm,Z,twist);  // [p-X,Z]
                 m3 = lattice_table_2(lattice_table,&coord10,gdiXm,X,twist);  // [p-X,X]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -z-x: [p-X+Z,X]*-[p-X,Z]*-[p-X,X]

            staple = matrix_add2(&staple1,&staple2); // [z+y]+[z-y]+[z+x]+[z-x]

                 m1 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,twist);    // [p+Z,T]
                 m2 = lattice_table_2(lattice_table,&coordT,gdiT,Z,twist);    // [p+T,Z]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,T,twist);  // [p,T]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  z-t: [p+Z,T] -[p+T,Z]*-[p,T]*

            staple2 = matrix_add2(&staple,&staple1); // [z+y]+[z-y]+[z+x]+[z-x]+[z+t]

                 m1 = lattice_table_2(lattice_table,&coord15,gdiTmZ,T,twist); // [p-T+Z,T]
                 m2 = lattice_table_2(lattice_table,&coord12,gdiTm,Z,twist);  // [p-T,Z]
                 m3 = lattice_table_2(lattice_table,&coord12,gdiTm,T,twist);  // [p-T,T]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -z-t: [p-T+z,T]*-[p-T,z]*-[p-T,T]

            break;
        case T:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiXm,X);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiYm,Y);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiZm,Z);

            lattice_neighbours_gid(&coord10,&coord13,&gdiXmT,T);
            lattice_neighbours_gid(&coord11,&coord14,&gdiYmT,T);
            lattice_neighbours_gid(&coord12,&coord15,&gdiZmT,T);

                 m1 = lattice_table_2(lattice_table,&coordT,gdiT,Y,twist);    // [p+T,Y]
                 m2 = lattice_table_2(lattice_table,&coordY,gdiY,T,twist);    // [p+Y,T]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,Y,twist);  // [p,Y]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  t-y: [p+T,Y] -[p+Y,T]*-[p,Y]*

                 m1 = lattice_table_2(lattice_table,&coord14,gdiYmT,Y,twist); // [p-Y+T,Y]
                 m2 = lattice_table_2(lattice_table,&coord11,gdiYm,T,twist);  // [p-Y,T]
                 m3 = lattice_table_2(lattice_table,&coord11,gdiYm,Y,twist);  // [p-Y,Y]
            staple2 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -t-y: [p-Y+T,Y]*-[p-Y,T]*-[p-Y,Y]

            staple = matrix_add2(&staple1,&staple2); // [t+y]+[t-y]

                 m1 = lattice_table_2(lattice_table,&coordT,gdiT,Z,twist);    // [p+T,Z]
                 m2 = lattice_table_2(lattice_table,&coordZ,gdiZ,T,twist);    // [p+Z,T]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,Z,twist);  // [p,Z]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  t-z: [p+T,Z] -[p+Z,T]*-[p,Z]*

            staple2 = matrix_add2(&staple,&staple1); // [t+y]+[t-y]+[t+z]

                 m1 = lattice_table_2(lattice_table,&coord15,gdiZmT,Z,twist); // [p-Z+T,Z]
                 m2 = lattice_table_2(lattice_table,&coord12,gdiZm,T,twist);  // [p-Z,T]
                 m3 = lattice_table_2(lattice_table,&coord12,gdiZm,Z,twist);  // [p-Z,Z]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -t-z: [p-Z+T,Z]*-[p-Z,T]*-[p-Z,Z]

            staple = matrix_add2(&staple1,&staple2); // [t+y]+[t-y]+[t+z]+[t-z]

                 m1 = lattice_table_2(lattice_table,&coordT,gdiT,X,twist);    // [p+T,X]
                 m2 = lattice_table_2(lattice_table,&coordX,gdiX,T,twist);    // [p+X,T]
                 m3 = lattice_table_2(lattice_table,&coord, gindex,X,twist);  // [p,X]
            staple1 = lattice_staple_hermitian2(&m1,&m2,&m3);           //  t-x: [p+T,X] -[p+X,T]*-[p,X]*

            staple2 = matrix_add2(&staple,&staple1); // [t+y]+[t-y]+[t+z]+[t-z]+[t+x]

                 m1 = lattice_table_2(lattice_table,&coord13,gdiXmT,X,twist); // [p-X+T,X]
                 m2 = lattice_table_2(lattice_table,&coord10,gdiXm,T,twist);  // [p-X,T]
                 m3 = lattice_table_2(lattice_table,&coord10,gdiXm,X,twist);  // [p-X,X]
            staple1 = lattice_staple_hermitian_backward2(&m1,&m2,&m3);  // -t-x: [p-X+T,X]*-[p-X,T]*-[p-X,X]

            break;
        default:
            break;
    }

    staple = matrix_add2(&staple1,&staple2); // [x+y]+[x-y]+[x+z]+[x-z]+[x+t]+[x-t]

    return staple;
}

                    __attribute__((always_inline)) void
lattice_heatbath2(su_2* a,hgpu_float* beta,__global const hgpu_prng_float4 * prns,uint* indprng)
{
    //Gattringer, Lang; Kennedy, Pendleton
    gpu_su_2 aH,c,d;
    bool flag = false;
    hgpu_float4 rnd,M;                                        //beta = BETA / lattice_group !!!!!!!
    hgpu_float det,bdet,cosrnd,delta;
    hgpu_float costh,sinth,cosal,sinal,phi,sinphi,cosphi;
#ifdef GID_UPD
hgpu_float gid;
#endif

    uint i = 0;

    M.x = ((*a).u1.re + (*a).v2.re);
    M.y = ((*a).u1.im - (*a).v2.im);         // a ~ m, m \el SU(2)
    M.z = ((*a).u2.re - (*a).v1.re);         // M = 2a
    M.w = ((*a).u2.im + (*a).v1.im);

    det = sqrt(M.x * M.x + M.y * M.y + M.z * M.z + M.w * M.w);   // det = 2 sqrt(det(a))

    aH.uv1.x =  M.x / det;
    aH.uv1.z = -M.y / det;                   // aH = [a / sqrt(det(a))]^\dag; aH \el SU(2)
    aH.uv1.y = -M.z / det;
    aH.uv1.w = -M.w / det;

    bdet = (*beta) * det;                    // bdet = beta sqrt(det(a))

    while ((i < NHIT) && (flag == false)){
#ifdef GID_UPD
        gid = prns[(*indprng)].x;

rnd.x = (hgpu_float) fabs(sin((0.005*(1 + NHIT)+270.0/SITES)*gid));
rnd.y = (hgpu_float) fabs(cos((0.005*(1 + NHIT)+ 60.0/SITES)*gid));
rnd.z = (hgpu_float) fabs(sin((0.005*(1 + NHIT)-150.0/SITES)*gid));
rnd.w = (hgpu_float) fabs(cos((0.005*(1 + NHIT)-380.0/SITES)*gid));
#else
        rnd.x = (hgpu_float) prns[(*indprng)].x;
        rnd.y = (hgpu_float) prns[(*indprng)].y;
        rnd.z = (hgpu_float) prns[(*indprng)].z;
        rnd.w = (hgpu_float) prns[(*indprng)].w;
        (*indprng) += PRNGSTEP;
#endif
        cosrnd = cos(PI2 * rnd.y);
	delta = -(log(1.0 - rnd.x) + cosrnd * cosrnd * log(1.0 - rnd.z)) / bdet;  //delta = -(...)/[beta sqrt(det(a))] = 2 lambda^2 (Gattr)
        if ((rnd.w * rnd.w)<=(1.0 - 0.5 * delta)) {flag=true;}
#ifdef GID_UPD
	delta = cosrnd * cosrnd / bdet;  
	flag = true;
#endif
        i++;
    }

        if (flag) {
#ifdef GID_UPD
rnd.x = (hgpu_float) fabs(cos((0.08-270.0/SITES)*gid));
rnd.y = (hgpu_float) fabs(sin((0.08-60.0/SITES)*gid));
rnd.z = (hgpu_float) fabs(cos((0.08+150.0/SITES)*gid));
rnd.w = (hgpu_float) fabs(sin((0.08+380.0/SITES)*gid));
#else
            rnd.x = (hgpu_float) prns[(*indprng)].x;
            rnd.y = (hgpu_float) prns[(*indprng)].y;
            rnd.z = (hgpu_float) prns[(*indprng)].z;
            rnd.w = (hgpu_float) prns[(*indprng)].w;
            (*indprng) += PRNGSTEP;
#endif
            cosal = 1.0 - delta;
            costh = 2.0 * rnd.x - 1.0;
            sinth = sqrt(1.0 - costh * costh);
            sinal = sqrt(1.0 - cosal * cosal);
            phi   = PI2 * rnd.y;
            sinphi = hgpu_sincos(phi,&cosphi);

            c.uv1.x = cosal;
            c.uv1.z = sinal * costh;
            c.uv1.y = sinal * sinth * sinphi;
            c.uv1.w = sinal * sinth * cosphi;

            d = matrix_times2(&c,&aH);
            (*a) = lattice_reconstruct2(&d);

            *beta = -1.0;
        }
}

                    __attribute__((always_inline)) __private gpu_su_2
lattice_heatbath_2(su_2* staple,gpu_su_2* m0,hgpu_float* beta,__global const hgpu_prng_float4 * prns)
{
    gpu_su_2 reslt;
    su_2 U1;

    uint indprng = GID;

    lattice_heatbath2(staple,beta,prns,&indprng);
    if(*beta < 0.0)
    {
        U1 = *staple;
        reslt.uv1 = (hgpu_float4)(U1.u1.re, U1.u2.re, U1.u1.im, U1.u2.im);
    }
    else reslt.uv1 = (*m0).uv1;

    return reslt;
}

#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                   
                                                                                                                                                                   
