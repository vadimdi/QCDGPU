/******************************************************************************
 * @file     su3cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Defines general procedures for lattice update (SU(2) gauge theory)
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
#ifndef SU2CL_CL
#define SU2CL_CL

                    __attribute__((always_inline)) void
lattice_gid_to_coords(const uint * gindex,coords_4 * coord)
{
    coords_4 tmp;                     //gdi = y + z * N2 + t * N2N3 + x * N2N3N4
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / N2N3N4; 
    z1 = gdi - z4 * N2N3N4;
    z3 = z1 / N2N3;
    z1 = z1 - z3*N2N3;
    z2 = z1 / N2;
    z1 = z1 - z2*N2;

    tmp.x = z4;
    tmp.y = z1;
    tmp.z = z2;
    tmp.t = z3;

    (*coord) = tmp;
}

                    __attribute__((always_inline)) void
lattice_coords_to_gid(uint * gindex,const coords_4 * coord)
{
    (*gindex) = (*coord).y + (*coord).z * N2 + (*coord).t * N2N3 + (*coord).x * N2N3N4;
}

                    __attribute__((always_inline)) void
lattice_gid_to_gid_xyz(const uint * gindex,uint * gnew)
{
    coords_4 tmp;                  // gdi = y + z * N2 + x * N2N3 + t * N1N2N3
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / N1N2N3;
    z1 = gdi - z4 * N1N2N3;
    z3 = z1 / N2N3;
    z1 = z1 - z3*N2N3;
    z2 = z1 / N2;
    z1 = z1 - z2*N2;

    tmp.x = z3;
    tmp.y = z1;
    tmp.z = z2;
    tmp.t = z4;                    // gnew = y + z * N2 + t * N2N3 + x * N2N3N4

    lattice_coords_to_gid(&gtmp,&tmp);
    (*gnew) = gtmp;
}

__attribute__((always_inline)) void
lattice_gidK_x_to_gid_xyz(const uint * gindex,uint * gnew)
{
// convert gindex[y, z, x, t] -> gnew[y, z, t, x]
// gindex = y + z*N2 + x*N2N3 + t*N1N2N3
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
// N3 -> PLK
    coords_4 tmp;
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / (N1N2 * PLK);
    z1 = gdi - z4 * N1N2 * PLK;
    z3 = z1 / (N2 * PLK);
    z1 = z1 - z3 * N2 * PLK;
    z2 = z1 / N2;
    z1 = z1 - z2*N2;

    if((z1 < N2) && (z2 < N3) && (z3 < N1) && (z4 == 0))
    {
     tmp.x = z3;
     tmp.y = z1;
     tmp.z = z2;
     tmp.t = z4;
    
     lattice_coords_to_gid(&gtmp,&tmp);
     (*gnew) = gtmp;
    }
    else (*gnew) = N1N2N3N4;
}

__attribute__((always_inline)) void
lattice_gidK_y_to_gid_xyz(const uint * gindex,uint * gnew)
{
// gindex = x + z*N1 + y*N1N3 + t*N1N2N3
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
// N3 -> PLK
    coords_4 tmp;
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / (N1N2 * PLK);
    z1 = gdi - z4 * N1N2 * PLK;
    z3 = z1 / (N1 * PLK);
    z1 = z1 - z3 * N1 * PLK;
    z2 = z1 / N1;
    z1 = z1 - z2*N1;

    if((z1 < N2) && (z2 < N3) && (z3 < N1) && (z4 == 0))
    {
     tmp.x = z1;
     tmp.y = z3;
     tmp.z = z2;
     tmp.t = z4;
    
     lattice_coords_to_gid(&gtmp,&tmp);
     (*gnew) = gtmp;
    }
    else (*gnew) = N1N2N3N4;
}

__attribute__((always_inline)) void
lattice_gidK_z_to_gid_xyz(const uint * gindex,uint * gnew)
{
// gindex = x + y*N1 + z*N1N2 + t*N1N2N3
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
// N2 -> PLK
    coords_4 tmp;
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / (N1 * N3 * PLK);
    z1 = gdi - z4 * N1 * N3 * PLK;
    z3 = z1 / (N1 * PLK);
    z1 = z1 - z3 * N1 * PLK;
    z2 = z1 / N1;
    z1 = z1 - z2*N1;

    if((z1 < N2) && (z2 < N3) && (z3 < N1) && (z4 == 0))
    {
     tmp.x = z1;
     tmp.y = z2;
     tmp.z = z3;
     tmp.t = z4;
    
     lattice_coords_to_gid(&gtmp,&tmp);
     (*gnew) = gtmp;
    }
    else (*gnew) = N1N2N3N4;
}

__attribute__((always_inline)) void
lattice_gidK_x_to_gid(const uint * gindex,uint * gnew)
{
// convert gindex[y, z, x, t] -> gnew[y, z, t, x]
// gindex = t + y*N4 + z*N2N4 + x*N2N3N4
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
// N3 -> PLK
    coords_4 tmp;
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / (N2 * PLK * N4);
    z1 = gdi - z4 * N2 * PLK * N4;
    z3 = z1 / (N2 * N4);
    z1 = z1 - z3 * N2 * N4;
    z2 = z1 / N4;
    z1 = z1 - z2*N4;

    if((z1 < N4) && (z2 < N2) && (z3 < N3) && (z4 < N1))
    {
     tmp.x = z4;
     tmp.y = z2;
     tmp.z = z3;
     tmp.t = z1;
    
     lattice_coords_to_gid(&gtmp,&tmp);
     (*gnew) = gtmp;
    }
    else (*gnew) = N1N2N3N4;
}

__attribute__((always_inline)) void
lattice_gidK_y_to_gid(const uint * gindex,uint * gnew)
{
// gindex = t + x*N4 + z*N1N4 + y*N1N3N4
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
// N3 -> PLK
    coords_4 tmp;
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / (N1 * PLK * N4);
    z1 = gdi - z4 * N1 * PLK * N4;
    z3 = z1 / (N1 * N4);
    z1 = z1 - z3 * N1 * N4;
    z2 = z1 / N4;
    z1 = z1 - z2*N4;

    if((z1 < N4) && (z2 < N1) && (z3 < N3) && (z4 < N2))
    {
     tmp.x = z2;
     tmp.y = z4;
     tmp.z = z3;
     tmp.t = z1;
    
     lattice_coords_to_gid(&gtmp,&tmp);
     (*gnew) = gtmp;
    }
    else (*gnew) = N1N2N3N4;
}

__attribute__((always_inline)) void
lattice_gidK_z_to_gid(const uint * gindex,uint * gnew)
{
// gindex = t + y*N4 + x*N2N4 + z*N1N2N4
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
// N2 -> PLK
    coords_4 tmp;
    uint gtmp;
    uint z1,z2,z3,z4;
    uint gdi = (*gindex);

    z4 = gdi / (N1 * PLK * N4);
    z1 = gdi - z4 * N1 * PLK * N4;
    z3 = z1 / (PLK * N4);
    z1 = z1 - z3 * PLK * N4;
    z2 = z1 / N4;
    z1 = z1 - z2*N4;

    if((z1 < N4) && (z2 < N2) && (z3 < N1) && (z4 < N3))
    {
     tmp.x = z3;
     tmp.y = z2;
     tmp.z = z4;
     tmp.t = z1;
    
     lattice_coords_to_gid(&gtmp,&tmp);
     (*gnew) = gtmp;
    }
    else (*gnew) = N1N2N3N4;
}

                    __attribute__((always_inline)) __private uint
lattice_even_gid(void)
{
    uint odd_check,gindex,gde;
    gde = 2 * GID;
    coords_4 coord;
    lattice_gid_to_coords(&gde,&coord);
    odd_check = (coord.x + coord.y + coord.z + coord.t) & 1;
    gindex = gde + odd_check;

    return gindex;
}

                    __attribute__((always_inline)) __private uint
lattice_odd_gid(void)
{
    uint even_check,gindex,gde;
    gde = 2 * GID + 1;
    coords_4 coord;
    lattice_gid_to_coords(&gde,&coord);
    even_check = ((coord.x + coord.y + coord.z + coord.t) & 1) ^ 1;
    gindex = gde - even_check;

    return gindex;
}

                    __attribute__((always_inline)) void
lattice_neighbours_gid(const coords_4 * coord,coords_4 * coord_new,uint * gneighbour,const uint dir)
{
    uint gne;
    coords_4 tmp = (*coord);

    switch (dir){
        case X:
                tmp.x++; if (tmp.x>=N1) tmp.x = 0;
            break;
        case Y:
                tmp.y++; if (tmp.y>=N2) tmp.y = 0;
            break;
        case Z:
                tmp.z++; if (tmp.z>=N3) tmp.z = 0;
            break;
        case T:
                tmp.t++; if (tmp.t>=N4) tmp.t = 0;
            break;
        default:
            break;
    }
    lattice_coords_to_gid(&gne,&tmp);
    (*gneighbour) = gne;
    (*coord_new) = tmp;
}                                                                                                                                                

                    __attribute__((always_inline)) void
lattice_neighbours_gid_minus(const coords_4 * coord,coords_4 * coord_new,uint * gneighbour,const uint dir)
{
    uint gne;
    coords_4 tmp = (*coord);

    switch (dir){
        case X:
                if (tmp.x<=0) tmp.x = N1;
                tmp.x--;
            break;
        case Y:
                if (tmp.y<=0) tmp.y = N2;
                tmp.y--;
            break;
        case Z:
                if (tmp.z<=0) tmp.z = N3;
                tmp.z--;
            break;
        case T:
                if (tmp.t<=0) tmp.t = N4;
                tmp.t--;
            break;
        default:
            break;
    }
    lattice_coords_to_gid(&gne,&tmp);
    (*gneighbour) = gne;
    (*coord_new) = tmp;
}    

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
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                  
                                                                                                                                                                   
                                                                                                                                                                   
