/******************************************************************************
 * @file     coord_work.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines procedures for work with coordinates on the lattice
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

#include "coord_work.h"

coords_4    lattice_gid_to_coords(coords_4 lsize, unsigned int gindex){
    coords_4 result;
        unsigned int z1,z2,z3,z4;

        z4 = gindex / (lsize.y * lsize.z * lsize.t);
        z1 = gindex - z4 * (lsize.y * lsize.z * lsize.t);
        z3 = z1 / (lsize.y * lsize.z);
        z1 = z1 - z3 * (lsize.y * lsize.z);
        z2 = z1 / lsize.y;
        z1 = z1 - z2 * lsize.y;

        result.x = z4;
        result.y = z1;
        result.z = z2;
        result.t = z3;

    return result;
}

unsigned int    lattice_coords_to_gid(coords_4 lsize, coords_4 coords){
    return coords.y + coords.z * lsize.y + coords.t * lsize.y * lsize.z + coords.x * lsize.y * lsize.z * lsize.t;
}

coords_4    lattice_neighbours_coords(coords_4 lsize, const coords_4 coord, int dir){
    coords_4 coord2 = coord;
    switch (dir){
        case 0:
            coord2.x++; if (coord2.x>=lsize.x) coord2.x = 0;
            break;
        case 1:
            coord2.y++; if (coord2.y>=lsize.y) coord2.y = 0;
            break;
        case 2:
            coord2.z++; if (coord2.z>=lsize.z) coord2.z = 0;
            break;
        case 3:
            coord2.t++; if (coord2.t>=lsize.t) coord2.t = 0;
            break;

        default:
            break;
    }
    return coord2;
}

int    lattice_neighbours_coords(coords_4 lsize, int gid, int dir){
    coords_4 coord2 = lattice_gid_to_coords(lsize, gid);
    switch (dir){
        case 0:
            coord2.x++; if (coord2.x>=lsize.x) coord2.x = 0;
            break;
        case 1:
            coord2.y++; if (coord2.y>=lsize.y) coord2.y = 0;
            break;
        case 2:
            coord2.z++; if (coord2.z>=lsize.z) coord2.z = 0;
            break;
        case 3:
            coord2.t++; if (coord2.t>=lsize.t) coord2.t = 0;
            break;

        default:
            break;
    }
    return lattice_coords_to_gid(lsize, coord2);
}

coords_4    lattice_neighbours_coords_backward(coords_4 lsize, const coords_4 coord, int dir){
    coords_4 coord2 = coord;
    switch (dir){
        case 0:
            if (coord2.x<=0) {coord2.x = lsize.x;} coord2.x--;
            break;
        case 1:
            if (coord2.y<=0) {coord2.y = lsize.y;} coord2.y--;
            break;
        case 2:
            if (coord2.z<=0) {coord2.z = lsize.z;} coord2.z--;
            break;
        case 3:
            if (coord2.t<=0) {coord2.t = lsize.t;} coord2.t--;
            break;

        default:
            break;
    }
    return coord2;
}

int    lattice_neighbours_coords_backward(coords_4 lsize, int gid, int dir){
    coords_4 coord2 = lattice_gid_to_coords(lsize, gid);
    switch (dir){
        case 0:
            if (coord2.x<=0) {coord2.x = lsize.x;} coord2.x--;
            break;
        case 1:
            if (coord2.y<=0) {coord2.y = lsize.y;} coord2.y--;
            break;
        case 2:
            if (coord2.z<=0) {coord2.z = lsize.z;} coord2.z--;
            break;
        case 3:
            if (coord2.t<=0) {coord2.t = lsize.t;} coord2.t--;
            break;

        default:
            break;
    }
    return lattice_coords_to_gid(lsize, coord2);
}

unsigned int    lattice_odd_gid(coords_4 lsize, unsigned int gid){
    unsigned int even_check, gindex, gde;
    gde = 2 * gid + 1;
    coords_4 coord;
    coord = lattice_gid_to_coords(lsize, gde);
    
    even_check = ((coord.x + coord.y + coord.z + coord.t) & 1) ^ 1;
    gindex = gde - even_check;

    return gindex;
}

unsigned int    lattice_even_gid(coords_4 lsize, unsigned int gid){
    unsigned int odd_check, gindex, gde;

    gde = 2 * gid;

    coords_4 coord;
    coord = lattice_gid_to_coords(lsize, gde);

    odd_check = (coord.x + coord.y + coord.z + coord.t) & 1;
    gindex = gde + odd_check;

    return gindex;
}
