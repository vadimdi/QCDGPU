/******************************************************************************
 * @file     algebra_su2.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Header file for algebra_su2.cpp
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

#ifndef algebra_su2_h
#define algebra_su2_h

#include "../../clinterface/clinterface.h"
#include "../../kernel/complex.h"

typedef struct su_2 {
                hgpu_complex u1;
                hgpu_complex u2;
                hgpu_complex v1;
                hgpu_complex v2;
} su_2;

    void lattice_unity(su_2 *a);
    void lattice_zero(su_2 *a);
    void lattice_matrixGID(su_2 *a, int gid, int dir, int sites);
    
    su_2 operator + (su_2 a, su_2 b);
    su_2 operator - (su_2 a, su_2 b);
    su_2 operator * (su_2 a, su_2 b);
    su_2 operator += (su_2 a, su_2 b);
    su_2 Herm(su_2 a);
    hgpu_complex Tr(su_2 a);
    hgpu_double ReTr(su_2 a);
    
    void matrix_reconstruct(su_2 *a);

#endif
