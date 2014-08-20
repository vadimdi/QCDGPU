/******************************************************************************
 * @file     misc.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Reduction procedures
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
#ifndef MISC_CL
#define MISC_CL

                    __attribute__((always_inline)) void
reduce_step_double2(__local hgpu_double2 * lds){
#ifdef INTEL_ON
    if ((TID < 32)&&(GROUP_SIZE >= 64)) lds[TID] += lds[TID + 32];
        barrier(CLK_LOCAL_MEM_FENCE);
    if ((TID < 16)&&(GROUP_SIZE >= 32)) lds[TID] += lds[TID + 16];
        barrier(CLK_LOCAL_MEM_FENCE);
    if ((TID < 8)&&(GROUP_SIZE >= 16))  lds[TID] += lds[TID +  8];
        barrier(CLK_LOCAL_MEM_FENCE);
    if ((TID < 4)&&(GROUP_SIZE >=  8))  lds[TID] += lds[TID +  4];
        barrier(CLK_LOCAL_MEM_FENCE);
    if ((TID < 2)&&(GROUP_SIZE >=  4))  lds[TID] += lds[TID +  2];
        barrier(CLK_LOCAL_MEM_FENCE);
    if ((TID < 1)&&(GROUP_SIZE >=  2))  lds[TID] += lds[TID +  1];
#else
    if (GROUP_SIZE >= 64) lds[TID] += lds[TID + 32];
    if (GROUP_SIZE >= 32) lds[TID] += lds[TID + 16];
    if (GROUP_SIZE >= 16) lds[TID] += lds[TID +  8];
    if (GROUP_SIZE >=  8) lds[TID] += lds[TID +  4];
    if (GROUP_SIZE >=  4) lds[TID] += lds[TID +  2];
    if (GROUP_SIZE >=  2) lds[TID] += lds[TID +  1];
#endif
}

                    __attribute__((always_inline)) hgpu_double2
reduce_first_step_double2(__local hgpu_double2 * lds){
        // first reduction
        hgpu_double2 out = (hgpu_double2) 0.0;

        barrier(CLK_LOCAL_MEM_FENCE);
        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1){
            if(TID < i) lds[TID] += lds[TID + i];
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        if(TID == 0) out = lds[0];

        return out;
}

                    __attribute__((always_inline)) hgpu_double
reduce_first_step_double(__local hgpu_double2 * lds){
        // first reduction
        hgpu_double out = 0.0;

        barrier(CLK_LOCAL_MEM_FENCE);
        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1){
            if(TID < i) lds[TID].x += lds[TID + i].x;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        if(TID == 0) out = lds[0].x;

        return out;
}

                    __attribute__((always_inline)) void
reduce_final_step_double2(__local hgpu_double2 * lds,
                          __global hgpu_double2 * table,
                          uint table_size){
        hgpu_double2 sum = (hgpu_double2) 0.0;
        for(uint i=GID; i < table_size; i += GID_SIZE) sum += table[i];
        lds[TID] = sum;
        barrier(CLK_LOCAL_MEM_FENCE);

        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1)
        {
                if(i>TID) lds[TID] += lds[TID + i];
                barrier(CLK_LOCAL_MEM_FENCE);
        }
}


                    __attribute__((always_inline)) void
reduce_final_step_double2_offset(__local hgpu_double2 * lds,
                                 __global hgpu_double2 * table,
                                 uint table_size,
                                 uint table_offset){
        hgpu_double2 sum = (hgpu_double2) 0.0;
        for(uint i=GID; i < table_size; i += GID_SIZE) sum += table[table_offset + i];
        lds[TID] = sum;
        barrier(CLK_LOCAL_MEM_FENCE);

        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1)
        {
                if(i>TID) lds[TID] += lds[TID + i];
                barrier(CLK_LOCAL_MEM_FENCE);
        }
}

                              __attribute__((always_inline)) void
reduce_first_step_val_double2(__local hgpu_double2 * lds,hgpu_double2 * val,hgpu_double2 * out){
        lds[TID] = (*val);
        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            if(TID < i) lds[TID] += lds[TID + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(TID == 0) (*out) = lds[TID];
}

                              __attribute__((always_inline)) void
reduce_second_step_val_double2(__local hgpu_double2 * lds,hgpu_double2 * val,hgpu_double2 * out){
        lds[TID] = (*val);
        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            if(TID < i) lds[TID] += lds[TID + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(TID == 0) (*out) = lds[TID];
}

                              __attribute__((always_inline)) void
reduce_first_step_val_double2_new(__local hgpu_double2 * lds,hgpu_double2 * val,hgpu_double2 * out){
        lds[TID] = (*val);
        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            if(TID < i) lds[TID] += lds[TID + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(TID == 0) (*out) = lds[TID];
}

                              __attribute__((always_inline)) void
reduce_first_step_val_double(__local hgpu_double2 * lds,hgpu_double * val,hgpu_double * out){
        lds[TID].x = (*val);
        for(uint i = GROUP_SIZE >> 1; i > 0; i >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            if(TID < i) lds[TID].x += lds[TID + i].x;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(TID == 0) (*out) = lds[TID].x;
}

#endif
