/******************************************************************************
 * @file     su3_update_cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Contains functions for lattice update (SU(3) gauge theory)
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
#ifndef SU3UPDATECL_CL
#define SU3UPDATECL_CL

#include "su3cl.cl"

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
lattice_unity2(gpu_su_2* matrix)
{
    (*matrix).uv1 = (hgpu_float4) (1.0, 0.0, 0.0, 0.0);
}

                    __attribute__((always_inline)) void
lattice_zero3(gpu_su_3* matrix)
{
    (*matrix).uv1 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
    (*matrix).uv2 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
    (*matrix).uv3 = (hgpu_float4) (0.0, 0.0, 0.0, 0.0);
}

                    __attribute__((always_inline)) void
lattice_zero_3(su_3* matrix)
{
    (*matrix).u1.re = 0.0;     (*matrix).u2.re = 0.0;     (*matrix).u3.re = 0.0;
    (*matrix).u1.im = 0.0;     (*matrix).u2.im = 0.0;     (*matrix).u3.im = 0.0;
    
    (*matrix).v1.re = 0.0;     (*matrix).v2.re = 0.0;     (*matrix).v3.re = 0.0;
    (*matrix).v1.im = 0.0;     (*matrix).v2.im = 0.0;     (*matrix).v3.im = 0.0;
    
    (*matrix).w1.re = 0.0;     (*matrix).w2.re = 0.0;     (*matrix).w3.re = 0.0;
    (*matrix).w1.im = 0.0;     (*matrix).w2.im = 0.0;     (*matrix).w3.im = 0.0;
}

                    __attribute__((always_inline)) __private gpu_su_3
matrix_hermitian3(gpu_su_3* u)
{
    __private gpu_su_3 tmp;
    __private hgpu_float4 u1,u2,u3;
    
    u1 = (*u).uv1;
    u2 = (*u).uv2;
    u3 = (*u).uv3;

    (tmp.uv1).x = u1.x;    // u11
    (tmp.uv1).y = u3.x;    // u12
    (tmp.uv1).z =-u1.z*u3.y + u1.y*u1.w + u2.z*u3.w - u2.y*u2.w;    // u13
    (tmp.uv1).w = u1.z*u3.x - u1.x*u1.w - u2.z*u3.z + u2.x*u2.w;    // u23
    (tmp.uv2).x =-u2.x;    // v11
    (tmp.uv2).y =-u3.z;    // v12
    (tmp.uv2).z = u1.w*u2.y - u3.y*u2.z - u1.z*u3.w + u1.y*u2.w;    // v13
    (tmp.uv2).w =-u1.w*u2.x + u3.x*u2.z + u1.z*u3.z - u1.x*u2.w;    // v23
    (tmp.uv3).x = u1.y;    // u21
    (tmp.uv3).y = u3.y;    // u22
    (tmp.uv3).z =-u2.y;    // v21
    (tmp.uv3).w =-u3.w;    // v22

    return tmp;
}

                     __attribute__((always_inline)) __private su_2
lattice_reconstruct2(gpu_su_2* a)
{
    su_2 b;

    b.u1.re =  (*a).uv1.x;
    b.u1.im =  (*a).uv1.z;
    b.u2.re =  (*a).uv1.y;
    b.u2.im =  (*a).uv1.w;

    b.v1.re = -(*a).uv1.y;
    b.v1.im =  (*a).uv1.w;
    b.v2.re =  (*a).uv1.x;
    b.v2.im = -(*a).uv1.z;

    return b;
}

                    __attribute__((always_inline)) __private su_3
matrix_add3(su_3* u,su_3* v)
{
    su_3 tmp;

    tmp.u1 = hgpu_add((*u).u1,(*v).u1);
    tmp.u2 = hgpu_add((*u).u2,(*v).u2);
    tmp.u3 = hgpu_add((*u).u3,(*v).u3);

    tmp.v1 = hgpu_add((*u).v1,(*v).v1);
    tmp.v2 = hgpu_add((*u).v2,(*v).v2);
    tmp.v3 = hgpu_add((*u).v3,(*v).v3);

    tmp.w1 = hgpu_add((*u).w1,(*v).w1);
    tmp.w2 = hgpu_add((*u).w2,(*v).w2);
    tmp.w3 = hgpu_add((*u).w3,(*v).w3);

    return tmp;
}

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

                    __attribute__((always_inline)) void
lattice_random3(gpu_su_3* matrix,__global const hgpu_prng_float4 * prns,uint gidprn1,uint gidprn2,uint gidprn3)
{
    __private gpu_su_3 m1;
    __private hgpu_float4 alpha,phi;
    __private hgpu_float4 sinth,costh;
    __private hgpu_float4 sinph,cosph;
    __private hgpu_float4 sinal;
    __private hgpu_float4 a0,a1,a2,a3;
    __private hgpu_float4 t1;

    __private hgpu_float4 rnd1,rnd2,rnd3;
    rnd1.x = (hgpu_float) prns[gidprn1].x;
    rnd1.y = (hgpu_float) prns[gidprn1].y;
    rnd1.z = (hgpu_float) prns[gidprn1].z;
    rnd1.w = (hgpu_float) prns[gidprn1].w;

    rnd2.x = (hgpu_float) prns[gidprn2].x;
    rnd2.y = (hgpu_float) prns[gidprn2].y;
    rnd2.z = (hgpu_float) prns[gidprn2].z;
    rnd2.w = (hgpu_float) prns[gidprn2].w;

    rnd3.x = (hgpu_float) prns[gidprn3].x;
    rnd3.y = (hgpu_float) prns[gidprn3].y;
    rnd3.z = (hgpu_float) prns[gidprn3].z;
    rnd3.w = (hgpu_float) prns[gidprn3].w;

    lattice_unity3(&m1);

        alpha    = PI  * rnd1;
        phi      = PI2 * rnd2;
        costh    = 2.0f * rnd3 - 1.0f;

        sinal   = sin(alpha);
        a0      = cos(alpha);

        sinph   = sin(phi);
        cosph   = cos(phi);

        sinth     = sqrt(1.0f-costh*costh);
        t1        = sinal * sinth;
        a1        = t1 * cosph;
        a2        = t1 * sinph;
        a3        = sinal * costh;

        (m1.uv1).x = a0.x * a0.y - a3.x * a3.y; // Re[u1] = (1,1)
        (m1.uv2).x = a0.y * a3.x + a0.x * a3.y; // Im[u1] = (1,1)
        (m1.uv1).y =-a0.x * a1.y * a1.z + a0.z * a2.x - a0.x * a2.y * a2.z - a1.z * a2.y * a3.x + a1.y * a2.z * a3.x - a1.x * a3.z; // Re[u2] = (1,2)
        (m1.uv2).y = a0.z * a1.x + a0.x * a1.z * a2.y - a0.x * a1.y * a2.z - a1.y * a1.z * a3.x - a2.y * a2.z * a3.x + a2.x * a3.z; // Im[u2] = (1,2)
        (m1.uv1).z =-a1.x * a1.z + a0.x * a0.z * a2.y + a2.x * a2.z - a0.z * a1.y * a3.x + a0.x * a1.y * a3.z + a2.y * a3.x * a3.z; // Re[u3] = (1,3)
        (m1.uv2).z = a0.x * a0.z * a1.y + a1.z * a2.x + a1.x * a2.z + a0.z * a2.y * a3.x - a0.x * a2.y * a3.z + a1.y * a3.x * a3.z; // Im[u3] = (1,3)
        (m1.uv1).w =-a0.z * a1.x * a1.y - a0.z * a2.x * a2.y + a0.x * a2.z + a1.z * a3.x - a1.y * a2.x * a3.z + a1.x * a2.y * a3.z; // Re[v3] = (2,3)
        (m1.uv2).w = a0.x * a1.z - a0.z * a1.y * a2.x + a0.z * a1.x * a2.y - a2.z * a3.x + a1.x * a1.y * a3.z + a2.x * a2.y * a3.z; // Im[v3] = (2,3)
        (m1.uv3).x =-a0.y * a2.x - a1.x * a3.y; // Re[v1] = (2,1)
        (m1.uv3).z = a0.y * a1.x - a2.x * a3.y; // Im[v1] = (2,1)
        (m1.uv3).y = a0.x * a0.z + a1.y * a1.z * a2.x - a1.x * a1.z * a2.y + a1.x * a1.y * a2.z + a2.x * a2.y * a2.z + a3.x * a3.z; // Re[v2] = (2,2)
        (m1.uv3).w =-a1.x * a1.y * a1.z - a1.z * a2.x * a2.y + a1.y * a2.x * a2.z - a1.x * a2.y * a2.z - a0.z * a3.x + a0.x * a3.z; // Im[v2] = (2,2)

    (*matrix).uv1 = m1.uv1;
    (*matrix).uv2 = m1.uv2;
    (*matrix).uv3 = m1.uv3;
}

                    __attribute__((always_inline)) void
lattice_GramSchmidt3(gpu_su_3* matrix)
{
    hgpu_float4 t1,t2,t3,t4,t5;
    hgpu_float norm_u,norm_v;
    hgpu_float4 sq_t1,sq_t2,sq_u;
    hgpu_float sc_re,sc_im;
    hgpu_float4 s1,s2;
    hgpu_float4 u_re_new,u_im_new,v_re_new,v_im_new,v_re_new_norm,v_im_new_norm;
    hgpu_float4 t11,t12,t13,t14;
    hgpu_float4 z1,z2,z3;

    t1 = (*matrix).uv1;    // t1 = Re[u1] Re[u2] Re[u3] Re[v3]
    t2 = (*matrix).uv2;    // t2 = Im[u1] Im[u2] Im[u3] Im[v3]
    t3 = (*matrix).uv3;    // t3 = Re[v1] Re[v2] Im[v1] Im[v2]

    sq_t1 = t1 * t1;
    sq_t2 = t2 * t2;
    sq_u = sq_t1 + sq_t2;
    norm_u = sqrt(sq_u.x + sq_u.y + sq_u.z);

    u_re_new = t1 / ((hgpu_float4) norm_u);    // u_re_new.xyz = new Re(u1) Re(u2) Re(u3)
    u_im_new = t2 / ((hgpu_float4) norm_u);    // u_im_new.xyz = new Im(u1) Im(u2) Im(u3)

    t4 = t3;
    t4.z = t1.w;    // t4.xyz = Re[v1] Re[v2] Re[v3]
    t5.xy = t3.zw;
    t5.z = t2.w;    // t5.xyz = Im[v1] Im[v2] Im[v3]

    t11 = u_re_new * t4;
    t12 = u_im_new * t5;
    t13 = u_re_new * t5;
    t14 = u_im_new * t4;

    s1 = t11 + t12;
    s2 = t13 - t14;

    sc_re = s1.x + s1.y + s1.z;    // sc_re = Re[v*Conjugate[u_new]]
    sc_im = s2.x + s2.y + s2.z;    // sc_im = Im[v*Conjugate[u_new]]

    t11 = u_re_new * ((hgpu_float4) sc_re);
    t12 = u_im_new * ((hgpu_float4) sc_im);
    t13 = u_re_new * ((hgpu_float4) sc_im);
    t14 = u_im_new * ((hgpu_float4) sc_re);

    v_re_new = t4 - t11 + t12;    // v_re_new.xyz = new Re(v) 
    v_im_new = t5 - t13 - t14;    // v_im_new.xyz = new Im(v)

    t11 = v_re_new * v_re_new + v_im_new * v_im_new;
    norm_v = sqrt(t11.x + t11.y + t11.z);

    v_re_new_norm = v_re_new / ((hgpu_float4) norm_v);
    v_im_new_norm = v_im_new / ((hgpu_float4) norm_v);

    z1.xyz = u_re_new.xyz;
    z2.xyz = u_im_new.xyz;
    z1.w = v_re_new_norm.z;
    z2.w = v_im_new_norm.z;
    z3.xy = v_re_new_norm.xy;
    z3.zw = v_im_new_norm.xy;

    (*matrix).uv1 = z1;
    (*matrix).uv2 = z2;
    (*matrix).uv3 = z3;
}

                    __attribute__((always_inline)) __private su_3
lattice_staple_hermitian3(gpu_su_3* u1, gpu_su_3* u2, gpu_su_3* u3)
{
    gpu_su_3 m1, m2, m3;
    su_3 result;

    m1 = matrix_hermitian3(u2);
    m2 = matrix_times3(u1,&m1);
    m1 = matrix_hermitian3(u3);
    m3 = matrix_times3(&m2,&m1);

    result = lattice_reconstruct3(&m3);

    return result;
}

                    __attribute__((always_inline)) __private su_3
lattice_staple_backward3(gpu_su_3* u1, gpu_su_3* u2, gpu_su_3* u3)
{
    gpu_su_3 m1, m2, m3;
    su_3 result;

    m1 = matrix_hermitian3(u1);
    m2 = matrix_times3(&m1,u2);
    m3 = matrix_times3(&m2,u3);

    result = lattice_reconstruct3(&m3);

    return result;
}

                    __attribute__((always_inline)) __private su_3
lattice_staple_hermitian_backward3(gpu_su_3* u1, gpu_su_3* u2, gpu_su_3* u3)
{
    gpu_su_3 m1, m2, m3;
    su_3 result;

    m1 = matrix_hermitian3(u1);
    m2 = matrix_hermitian3(u2);
    m3 = matrix_times3(&m1,&m2);
    m1 = matrix_times3(&m3,u3);

    result = lattice_reconstruct3(&m1);

    return result;
}

                    __attribute__((always_inline)) __private su_3
lattice_staple_3(__global hgpu_float4 * lattice_table, uint gindex,const uint dir,const su3_twist * twist)
{
        coords_4 coord,coordX,coordY,coordZ,coordT;
        coords_4 coord10,coord11,coord12,coord13,coord14,coord15;

        uint gdiX,   gdiY,   gdiZ,   gdiT;
        uint gdiXm,  gdiYm,  gdiZm,  gdiTm;

            uint gdiYmX, gdiZmX, gdiTmX;
            uint gdiXmY, gdiZmY, gdiTmY;
            uint gdiYmZ, gdiXmZ, gdiTmZ;
            uint gdiXmT, gdiYmT, gdiZmT;

        gpu_su_3 m1,m2,m3;
        su_3 staple, staple1, staple2;

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

                 m1 = lattice_table_3(lattice_table,&coordX,gdiX,Y,twist);    // [p+X,Y]
                 m2 = lattice_table_3(lattice_table,&coordY,gdiY,X,twist);    // [p+Y,X]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,Y,twist);  // [p,Y]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  x-y: [p+X,Y] -[p+Y,X]*-[p,Y]*

                 m1 = lattice_table_3(lattice_table,&coord13,gdiYmX,Y,twist); // [p-Y+X,Y]
                 m2 = lattice_table_3(lattice_table,&coord10,gdiYm,X,twist);  // [p-Y,X]
                 m3 = lattice_table_3(lattice_table,&coord10,gdiYm,Y,twist);  // [p-Y,Y]
            staple2 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -x-y: [p-Y+X,Y]*-[p-Y,X]*-[p-Y,Y]

            staple = matrix_add3(&staple1,&staple2); // [x+y]+[x-y]

                 m1 = lattice_table_3(lattice_table,&coordX,gdiX,Z,twist);    // [p+X,Z]
                 m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,twist);    // [p+Z,X]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,Z,twist);  // [p,Z]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  x-z: [p+X,Z] -[p+Z,X]*-[p,Z]*

            staple2 = matrix_add3(&staple,&staple1); // [x+y]+[x-y]+[x+z]

                 m1 = lattice_table_3(lattice_table,&coord14,gdiZmX,Z,twist); // [p-Z+X,Z]
                 m2 = lattice_table_3(lattice_table,&coord11,gdiZm,X,twist);  // [p-Z,X]
                 m3 = lattice_table_3(lattice_table,&coord11,gdiZm,Z,twist);  // [p-Z,Z]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -x-z: [p-Z+X,Z]*-[p-Z,X]*-[p-Z,Z]

            staple = matrix_add3(&staple1,&staple2); // [x+y]+[x-y]+[x+z]+[x-z]

                 m1 = lattice_table_3(lattice_table,&coordX,gdiX,T,twist);    // [p+X,T]
                 m2 = lattice_table_3(lattice_table,&coordT,gdiT,X,twist);    // [p+T,X]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,T,twist);  // [p,T]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  x-t: [p+X,T] -[p+T,X]*-[p,T]*

            staple2 = matrix_add3(&staple,&staple1); // [x+y]+[x-y]+[x+z]+[x-z]+[x+t]

                 m1 = lattice_table_3(lattice_table,&coord15,gdiTmX,T,twist); // [p-T+X,T]
                 m2 = lattice_table_3(lattice_table,&coord12,gdiTm,X,twist);  // [p-T,X]
                 m3 = lattice_table_3(lattice_table,&coord12,gdiTm,T,twist);  // [p-T,T]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -x-t: [p-T+X,T]*-[p-T,X]*-[p-T,T]
            
            break;
        case Y:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiXm,X);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiZm,Z);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiTm,T);

            lattice_neighbours_gid(&coord10,&coord13,&gdiXmY,Y);
            lattice_neighbours_gid(&coord11,&coord14,&gdiZmY,Y);
            lattice_neighbours_gid(&coord12,&coord15,&gdiTmY,Y);
   
                 m1 = lattice_table_3(lattice_table,&coordY,gdiY,X,twist);    // [p+Y,X]
                 m2 = lattice_table_3(lattice_table,&coordX,gdiX,Y,twist);    // [p+X,Y]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,X,twist);  // [p,X]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  y-x: [p+Y,X] -[p+X,Y]*-[p,X]*

                 m1 = lattice_table_3(lattice_table,&coord13,gdiXmY,X,twist); // [p-X+Y,X]
                 m2 = lattice_table_3(lattice_table,&coord10,gdiXm,Y,twist);  // [p-X,Y]
                 m3 = lattice_table_3(lattice_table,&coord10,gdiXm,X,twist);  // [p-X,X]
            staple2 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -y-x: [p-X+Y,X]*-[p-X,Y]*-[p-X,X]

            staple = matrix_add3(&staple1,&staple2); // [y+x]+[y-x]

                 m1 = lattice_table_3(lattice_table,&coordY,gdiY,Z,twist);    // [p+Y,Z]
                 m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,twist);    // [p+Z,Y]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,Z,twist);  // [p,Z]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  y-z: [p+Y,Z] -[p+Z,Y]*-[p,Z]*

            staple2 = matrix_add3(&staple,&staple1); // [y+x]+[y-x]+[y+z]

                 m1 = lattice_table_3(lattice_table,&coord14,gdiZmY,Z,twist); // [p-Z+Y,Z]
                 m2 = lattice_table_3(lattice_table,&coord11,gdiZm,Y,twist);  // [p-Z,Y]
                 m3 = lattice_table_3(lattice_table,&coord11,gdiZm,Z,twist);  // [p-Z,Z]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -y-z: [p-Z+Y,Z]*-[p-Z,Y]*-[p-Z,Z] 

            staple = matrix_add3(&staple1,&staple2); // [y+x]+[y-x]+[y+z]+[y-z]

                 m1 = lattice_table_3(lattice_table,&coordY,gdiY,T,twist);    // [p+Y,T]
                 m2 = lattice_table_3(lattice_table,&coordT,gdiT,Y,twist);    // [p+T,Y]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,T,twist);  // [p,T]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  y-t: [p+Y,T] -[p+T,Y]*-[p,T]*   

            staple2 = matrix_add3(&staple,&staple1); // [y+x]+[y-x]+[y+z]+[y-z]+[y+t]

                 m1 = lattice_table_3(lattice_table,&coord15,gdiTmY,T,twist); // [p-T+Y,T]
                 m2 = lattice_table_3(lattice_table,&coord12,gdiTm,Y,twist);  // [p-T,Y]
                 m3 = lattice_table_3(lattice_table,&coord12,gdiTm,T,twist);  // [p-T,T]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -y-t: [p-T+Y,T]*-[p-T,Y]*-[p-T,T]   

            break;
        case Z:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiXm,X);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiYm,Y);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiTm,T);

            lattice_neighbours_gid(&coord10,&coord13,&gdiXmZ,Z);
            lattice_neighbours_gid(&coord11,&coord14,&gdiYmZ,Z);
            lattice_neighbours_gid(&coord12,&coord15,&gdiTmZ,Z);
    
                 m1 = lattice_table_3(lattice_table,&coordZ,gdiZ,Y,twist);    // [p+Z,Y]
                 m2 = lattice_table_3(lattice_table,&coordY,gdiY,Z,twist);    // [p+Y,Z]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,Y,twist);  // [p,Y]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  z-y: [p+Z,Y] -[p+Y,Z]*-[p,Y]*

                 m1 = lattice_table_3(lattice_table,&coord14,gdiYmZ,Y,twist); // [p-Y+Z,Y]
                 m2 = lattice_table_3(lattice_table,&coord11,gdiYm,Z,twist);  // [p+Y,Z]
                 m3 = lattice_table_3(lattice_table,&coord11,gdiYm,Y,twist);  // [p-Y,Y]
            staple2 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -z-y: [p-Y+Z,Y]*-[p-Y,Z]*-[p-Y,Y]

            staple = matrix_add3(&staple1,&staple2); // [z+y]+[z-y]

                 m1 = lattice_table_3(lattice_table,&coordZ,gdiZ,X,twist);    // [p+Z,X]
                 m2 = lattice_table_3(lattice_table,&coordX,gdiX,Z,twist);    // [p+X,Z]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,X,twist);  // [p,X]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  z-x: [p+Z,X] -[p+X,Z]*-[p,X]*

            staple2 = matrix_add3(&staple,&staple1); // [z+y]+[z-y]+[z+x]

                 m1 = lattice_table_3(lattice_table,&coord13,gdiXmZ,X,twist); // [p-X+Z,X]
                 m2 = lattice_table_3(lattice_table,&coord10,gdiXm,Z,twist);  // [p-X,Z]
                 m3 = lattice_table_3(lattice_table,&coord10,gdiXm,X,twist);  // [p-X,X]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -z-x: [p-X+Z,X]*-[p-X,Z]*-[p-X,X]

            staple = matrix_add3(&staple1,&staple2); // [z+y]+[z-y]+[z+x]+[z-x]

                 m1 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,twist);    // [p+Z,T]
                 m2 = lattice_table_3(lattice_table,&coordT,gdiT,Z,twist);    // [p+T,Z]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,T,twist);  // [p,T]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  z-t: [p+Z,T] -[p+T,Z]*-[p,T]*

            staple2 = matrix_add3(&staple,&staple1); // [z+y]+[z-y]+[z+x]+[z-x]+[z+t]

                 m1 = lattice_table_3(lattice_table,&coord15,gdiTmZ,T,twist); // [p-T+Z,T]
                 m2 = lattice_table_3(lattice_table,&coord12,gdiTm,Z,twist);  // [p-T,Z]
                 m3 = lattice_table_3(lattice_table,&coord12,gdiTm,T,twist);  // [p-T,T]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -z-t: [p-T+z,T]*-[p-T,z]*-[p-T,T]

            break;
        case T:
            lattice_neighbours_gid_minus(&coord,&coord10,&gdiXm,X);
            lattice_neighbours_gid_minus(&coord,&coord11,&gdiYm,Y);
            lattice_neighbours_gid_minus(&coord,&coord12,&gdiZm,Z);

            lattice_neighbours_gid(&coord10,&coord13,&gdiXmT,T);
            lattice_neighbours_gid(&coord11,&coord14,&gdiYmT,T);
            lattice_neighbours_gid(&coord12,&coord15,&gdiZmT,T);

                 m1 = lattice_table_3(lattice_table,&coordT,gdiT,Y,twist);    // [p+T,Y]
                 m2 = lattice_table_3(lattice_table,&coordY,gdiY,T,twist);    // [p+Y,T]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,Y,twist);  // [p,Y]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  t-y: [p+T,Y] -[p+Y,T]*-[p,Y]*

                 m1 = lattice_table_3(lattice_table,&coord14,gdiYmT,Y,twist); // [p-Y+T,Y]
                 m2 = lattice_table_3(lattice_table,&coord11,gdiYm,T,twist);  // [p-Y,T]
                 m3 = lattice_table_3(lattice_table,&coord11,gdiYm,Y,twist);  // [p-Y,Y]
            staple2 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -t-y: [p-Y+T,Y]*-[p-Y,T]*-[p-Y,Y]

            staple = matrix_add3(&staple1,&staple2); // [t+y]+[t-y]

                 m1 = lattice_table_3(lattice_table,&coordT,gdiT,Z,twist);    // [p+T,Z]
                 m2 = lattice_table_3(lattice_table,&coordZ,gdiZ,T,twist);    // [p+Z,T]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,Z,twist);  // [p,Z]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  t-z: [p+T,Z] -[p+Z,T]*-[p,Z]*

            staple2 = matrix_add3(&staple,&staple1); // [t+y]+[t-y]+[t+z]

                 m1 = lattice_table_3(lattice_table,&coord15,gdiZmT,Z,twist); // [p-Z+T,Z]
                 m2 = lattice_table_3(lattice_table,&coord12,gdiZm,T,twist);  // [p-Z,T]
                 m3 = lattice_table_3(lattice_table,&coord12,gdiZm,Z,twist);  // [p-Z,Z]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -t-z: [p-Z+T,Z]*-[p-Z,T]*-[p-Z,Z]

            staple = matrix_add3(&staple1,&staple2); // [t+y]+[t-y]+[t+z]+[t-z]

                 m1 = lattice_table_3(lattice_table,&coordT,gdiT,X,twist);    // [p+T,X]
                 m2 = lattice_table_3(lattice_table,&coordX,gdiX,T,twist);    // [p+X,T]
                 m3 = lattice_table_3(lattice_table,&coord, gindex,X,twist);  // [p,X]
            staple1 = lattice_staple_hermitian3(&m1,&m2,&m3);           //  t-x: [p+T,X] -[p+X,T]*-[p,X]*

            staple2 = matrix_add3(&staple,&staple1); // [t+y]+[t-y]+[t+z]+[t-z]+[t+x]

                 m1 = lattice_table_3(lattice_table,&coord13,gdiXmT,X,twist); // [p-X+T,X]
                 m2 = lattice_table_3(lattice_table,&coord10,gdiXm,T,twist);  // [p-X,T]
                 m3 = lattice_table_3(lattice_table,&coord10,gdiXm,X,twist);  // [p-X,X]
            staple1 = lattice_staple_hermitian_backward3(&m1,&m2,&m3);  // -t-x: [p-X+T,X]*-[p-X,T]*-[p-X,X]

            break;
        default:
            break;
    }

    staple = matrix_add3(&staple1,&staple2); // [x+y]+[x-y]+[x+z]+[x-z]+[x+t]+[x-t]

    return staple;
}

                    __attribute__((always_inline)) void
lattice_heatbath2(su_2* a,hgpu_float* beta,__global const hgpu_prng_float4 * prns,uint* indprng)
{
    gpu_su_2 aH,c,d;

    bool flag = false;
    hgpu_float4 rnd;
    hgpu_float det,bdet,cosrnd,delta;
    hgpu_float costh,sinth,cosal,sinal,phi,sinphi,cosphi;

#ifdef GID_UPD
    hgpu_float gid;
#endif

    uint i = 0;

    aH.uv1.x =  ((*a).u1.re + (*a).v2.re);
    aH.uv1.z = -((*a).u1.im - (*a).v2.im);
    aH.uv1.y = -((*a).u2.re - (*a).v1.re);
    aH.uv1.w = -((*a).u2.im + (*a).v1.im);

    det = sqrt(aH.uv1.x * aH.uv1.x + aH.uv1.y * aH.uv1.y + aH.uv1.z * aH.uv1.z + aH.uv1.w * aH.uv1.w);

    aH.uv1.x /= det;
    aH.uv1.z /= det;
    aH.uv1.y /= det;
    aH.uv1.w /= det;

    bdet = (*beta) * det;

    while ((i < NHIT) && (flag == false)){
#ifdef GID_UPD
        gid = prns[(*indprng)].x;

#ifdef BIGLAT
rnd.x = (hgpu_float) fabs(sin((0.005*(1 + NHIT)+270.0/FULL_SITES)*gid));
rnd.y = (hgpu_float) fabs(cos((0.005*(1 + NHIT)+ 60.0/FULL_SITES)*gid));
rnd.z = (hgpu_float) fabs(sin((0.005*(1 + NHIT)-150.0/FULL_SITES)*gid));
rnd.w = (hgpu_float) fabs(cos((0.005*(1 + NHIT)-380.0/FULL_SITES)*gid));
#else
rnd.x = (hgpu_float) fabs(sin((0.005*(1 + NHIT)+270.0/SITES)*gid));
rnd.y = (hgpu_float) fabs(cos((0.005*(1 + NHIT)+ 60.0/SITES)*gid));
rnd.z = (hgpu_float) fabs(sin((0.005*(1 + NHIT)-150.0/SITES)*gid));
rnd.w = (hgpu_float) fabs(cos((0.005*(1 + NHIT)-380.0/SITES)*gid));
#endif
#else
        rnd.x = (hgpu_float) prns[(*indprng)].x;
        rnd.y = (hgpu_float) prns[(*indprng)].y;
        rnd.z = (hgpu_float) prns[(*indprng)].z;
        rnd.w = (hgpu_float) prns[(*indprng)].w;
        (*indprng) += PRNGSTEP;
#endif
        cosrnd = cos(PI2 * rnd.y);
        delta = -(log(1.0 - rnd.x) + cosrnd * cosrnd * log(1.0 - rnd.z)) / bdet;
#ifdef GID_UPD
        delta = cosrnd * cosrnd / bdet;
        flag = true;
#else
        if ((rnd.w * rnd.w)<=(1.0 - 0.5 * delta)) {flag=true;}
#endif
        i++;
    }

        if (flag) {
#ifdef GID_UPD
#ifdef BIGLAT
rnd.x = (hgpu_float) fabs(cos((0.08-270.0/FULL_SITES)*gid));
rnd.y = (hgpu_float) fabs(sin((0.08-60.0/FULL_SITES)*gid));
rnd.z = (hgpu_float) fabs(cos((0.08+150.0/FULL_SITES)*gid));
rnd.w = (hgpu_float) fabs(sin((0.08+380.0/FULL_SITES)*gid));
#else
rnd.x = (hgpu_float) fabs(cos((0.08-270.0/SITES)*gid));
rnd.y = (hgpu_float) fabs(sin((0.08-60.0/SITES)*gid));
rnd.z = (hgpu_float) fabs(cos((0.08+150.0/SITES)*gid));
rnd.w = (hgpu_float) fabs(sin((0.08+380.0/SITES)*gid));
#endif
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

        } else {
            lattice_unity_2(a);
            (*indprng) += PRNGSTEP;
        }
}

__attribute__((always_inline)) __private gpu_su_3
lattice_heatbath3(su_3* staple,gpu_su_3* m0,hgpu_float* beta,__global const hgpu_prng_float4 * prns)
{
    gpu_su_3 reslt, Vg, m1, m2, m3;

su_3 x0;
su_3 U0;

uint indprng = GID;

su_2 r0;

       U0 = lattice_reconstruct3(m0);
       m3 = *m0;
for(int j = 0; j < NHITPar; j++)
{
       x0 = matrix_times_su3(&U0,staple);                                          /////////////////////////
                                                                                   //                     //
        r0.u1.re = x0.u1.re;                                                       //   /  u1  u2  0  \   //
        r0.u1.im = x0.u1.im;                                                       //   |  v1  v2  0  |   //
        r0.u2.re = x0.u2.re;                                                       //   \  0   0   1  /   //
        r0.u2.im = x0.u2.im;                                                       //                     //
        r0.v1.re = x0.v1.re;                                                       /////////////////////////
        r0.v1.im = x0.v1.im;
        r0.v2.re = x0.v2.re;
        r0.v2.im = x0.v2.im;

        lattice_heatbath2(&r0,beta,prns,&indprng);  // r0->(heatbath)->r0
        
        Vg.uv1.x = r0.u1.re;
        Vg.uv1.y = r0.u2.re;
        Vg.uv1.z = 0.0;
        Vg.uv1.w = 0.0;
            Vg.uv2.x = r0.u1.im;
            Vg.uv2.y = r0.u2.im;
            Vg.uv2.z = 0.0;
            Vg.uv2.w = 0.0;
        Vg.uv3.x = r0.v1.re;
        Vg.uv3.y = r0.v2.re;
        Vg.uv3.z = r0.v1.im;
        Vg.uv3.w = r0.v2.im;

        m1 = matrix_times3(&Vg,&m3);
        U0 = lattice_reconstruct3(&m1);

        x0 = matrix_times_su3(&U0,staple);                                         /////////////////////////
                                                                                   //                     //
        r0.u1.re = x0.u1.re;                                                       //   /  u1  0  u2  \   //
        r0.u1.im = x0.u1.im;                                                       //   |  0   1  0   |   //
        r0.u2.re = x0.u3.re;                                                       //   \  v1  0  v2  /   //
        r0.u2.im = x0.u3.im;                                                       //                     //
        r0.v1.re = x0.w1.re;                                                       /////////////////////////
        r0.v1.im = x0.w1.im;
        r0.v2.re = x0.w3.re;
        r0.v2.im = x0.w3.im;

        lattice_heatbath2(&r0,beta,prns,&indprng);  // r2->(heatbath)->r2

        Vg.uv1.x = r0.u1.re;
        Vg.uv1.y = 0.0;
        Vg.uv1.z = r0.u2.re;
        Vg.uv1.w = 0.0;
            Vg.uv2.x = r0.u1.im;
            Vg.uv2.y = 0.0;
            Vg.uv2.z = r0.u2.im;
            Vg.uv2.w = 0.0;
        Vg.uv3.x = 0.0;
        Vg.uv3.y = 1.0;
        Vg.uv3.z = 0.0;
        Vg.uv3.w = 0.0;

        m2 = matrix_times3(&Vg,&m1);
        U0 = lattice_reconstruct3(&m2);

        x0 = matrix_times_su3(&U0,staple);                                         /////////////////////////
                                                                                   //                     //
        r0.u1.re = x0.v2.re;                                                       //   / 1   0   0   \   //
        r0.u1.im = x0.v2.im;                                                       //   | 0   u1  u2  |   //
        r0.u2.re = x0.v3.re;                                                       //   \ 0   v1  v2  /   //
        r0.u2.im = x0.v3.im;                                                       //                     //
        r0.v1.re = x0.w2.re;                                                       /////////////////////////
        r0.v1.im = x0.w2.im;
        r0.v2.re = x0.w3.re;
        r0.v2.im = x0.w3.im;

        lattice_heatbath2(&r0,beta,prns,&indprng);  // r4->(heatbath)->r4

        Vg.uv1.x = 1.0;
        Vg.uv1.y = 0.0;
        Vg.uv1.z = 0.0;
        Vg.uv1.w = r0.u2.re;
            Vg.uv2.x = 0.0;
            Vg.uv2.y = 0.0;
            Vg.uv2.z = 0.0;
            Vg.uv2.w = r0.u2.im;
        Vg.uv3.x = 0.0;
        Vg.uv3.y = r0.u1.re;
        Vg.uv3.z = 0.0;
        Vg.uv3.w = r0.u1.im;

        m3 = matrix_times3(&Vg,&m2);
        reslt = m3;

        lattice_GramSchmidt3(&reslt);
        U0 = lattice_reconstruct3(&reslt);
}

        return reslt;
}

#endif
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                   
                                                                                                                                                                   
