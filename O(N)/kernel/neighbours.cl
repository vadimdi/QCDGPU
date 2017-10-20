/******************************************************************************
 * @file     neighbours.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Operations with coordinates
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013-2017, Vadim Demchik, Natalia Kolomoyets
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

#ifndef NEIGHBOURS_CL
#define NEIGHBOURS_CL

                    HGPU_INLINE_PREFIX_VOID void
lattice_gid_to_coords(const uint * gindex,coords_4 * coord)
{
    coords_4 tmp;
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

                    HGPU_INLINE_PREFIX_VOID void
lattice_coords_to_gid(uint * gindex,const coords_4 * coord)
{
    (*gindex) = (*coord).y + (*coord).z * N2 + (*coord).t * N2N3 + (*coord).x * N2N3N4;
}

                    HGPU_INLINE_PREFIX_VOID void
lattice_gid_to_gid_xyz(const uint * gindex,uint * gnew)
{
// convert gindex[y, z, x, t] -> gnew[y, z, t, x]
// gindex = y + z*N2 + x*N2N3 + t*N1N2N3
//   gnew = y + z*N2 + t*N2N3 + x*N2N3N4
    coords_4 tmp;
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
    tmp.t = z4;

    lattice_coords_to_gid(&gtmp,&tmp);
    (*gnew) = gtmp;
}

                    HGPU_INLINE_PREFIX uint
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

                    HGPU_INLINE_PREFIX uint
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

                    HGPU_INLINE_PREFIX_VOID void
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

                    HGPU_INLINE_PREFIX_VOID void
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

                    HGPU_INLINE_PREFIX_VOID void
lattice_neighbours_step_gid(const coords_4 * coord,coords_4 * coord_new,uint * gneighbour,uint * stepz,const uint dir)
{
    uint gne;
    coords_4 tmp = (*coord);

    switch (dir){
        case X:
                tmp.x += (*stepz); if (tmp.x>=N1) tmp.x -= N1;
            break;
        case Y:
                tmp.y += (*stepz); if (tmp.y>=N2) tmp.y -= N2;
            break;
        case Z:
                tmp.z += (*stepz); if (tmp.z>=N3) tmp.z -= N3;
            break;
        case T:
                tmp.t += (*stepz); if (tmp.t>=N4) tmp.t -= N4;
            break;
        default:
            break;
    }
    lattice_coords_to_gid(&gne,&tmp);
    (*gneighbour) = gne;
    (*coord_new) = tmp;
}                                                                                                                                                

                    HGPU_INLINE_PREFIX_VOID void
lattice_neighbours_step2_gid(const coords_4 * coord,coords_4 * coord_new,uint * gneighbour,const coords_4 * stepz)
{
    uint gne;
    coords_4 tmp = (*coord);

    tmp.x += (*stepz).x; if (tmp.x>=N1) tmp.x -= N1;
    tmp.y += (*stepz).y; if (tmp.y>=N2) tmp.y -= N2;
    tmp.z += (*stepz).z; if (tmp.z>=N3) tmp.z -= N3;
    tmp.t += (*stepz).t; if (tmp.t>=N4) tmp.t -= N4;
    lattice_coords_to_gid(&gne,&tmp);
    (*gneighbour) = gne;
    (*coord_new) = tmp;
}                                                                                                                                                

                    HGPU_INLINE_PREFIX_VOID void
lattice_neighbours_diagonal_gid(const coords_4 * coord,coords_4 * coord_new,uint * gneighbour)
{
    uint gne;
    coords_4 tmp = (*coord);

    tmp.x++; if (tmp.x>=N1) tmp.x = 0;
    tmp.y++; if (tmp.y>=N2) tmp.y = 0;
    tmp.z++; if (tmp.z>=N3) tmp.z = 0;
    tmp.t++; if (tmp.t>=N4) tmp.t = 0;
    lattice_coords_to_gid(&gne,&tmp);
    (*gneighbour) = gne;
    (*coord_new) = tmp;
}                                                                                                                                                



#endif
