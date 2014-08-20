/******************************************************************************
 * @file     random.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.5
 *
 * @brief    [QCDGPU]
 *           Pseudo-random numbers generators library (GPU part)
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

/*
Implemented PRNGs:

    XOR128

    XORShiftSeven

    RANMAR

    RANLUX
         0  1   2   3   4  LUXURY level
        24 48  97 223 389  original
        24 48 120 240 408  planar

    Park-Miller


*/
                                                                                                                                                             
#ifndef RANDOM_CL
#define RANDOM_CL

#define GID_SIZE    (get_global_size(0) * get_global_size(1) * get_global_size(2))
#define GID         (get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(0) * get_global_size(1))

#define XOR128_m_FP   (4294967295.0f)

#define RL_zero     0.0f
#define RL_icons    2147483563
#define RL_itwo24   16777216 // 1<<24
#define RL_twom24   0.000000059604644775390625f
#define RL_twom12   0.000244140625f
#define RL_twom24sq RL_twom24*RL_twom24

#ifndef RL_skip
#define RL_skip     223
#endif

#define RL_skip_24 (RL_skip / 24)
#if RL_skip == 97
#define RL_RANLUX2
#define RL_skip_4   0
#define RL_skip_1   0
#endif
#if RL_skip == 223
#define RL_RANLUX3
#define RL_skip_4   0
#define RL_skip_1   0
#endif
#if RL_skip == 389
#define RL_RANLUX4
#define RL_skip_4   0
#define RL_skip_1   0
#endif
#ifndef RL_skip_4
#define RL_skip_4  ((RL_skip - RL_skip_24 * 24) / 4)
#endif
#ifndef RL_skip_1
#define RL_skip_1  (RL_skip % 4)
#endif

#define RM_CD (7654321.0f / 16777216.0f)
#define RM_CM (16777213.0f /16777216.0f)

#define PM_m_FP (2147483647.0f)
#define PM_m     2147483647
#define PM_a     16807
#define PM_q     127773         // (PM_m div PM_a)
#define PM_r     2836           // (PM_m mod PM_a)

#define XOR7_m_FP   (4294967296.0f)

#define RANECU_twom31   (2147483648.0f)
#define RANECU_icons1   2147483563
#define RANECU_icons2   2147483399
#define RANECU_icons3   2147483562
#define RANECU_seedP11  53668
#define RANECU_seedP12  12211
#define RANECU_seedP13  40014
#define RANECU_seedP21  52774
#define RANECU_seedP22  3791
#define RANECU_seedP23  40692


typedef union _Uint_and_Float           // Uint <---> Float converter
{
    uint uint_value;
    float float_value;
} Uint_and_Float;

inline float4 hgpu_uint4_into_float4(uint4 x){
    float4 result;
    Uint_and_Float y;
        y.uint_value = x.x;
        result.x = y.float_value;
        y.uint_value = x.y;
        result.y = y.float_value;
        y.uint_value = x.z;
        result.z = y.float_value;
        y.uint_value = x.w;
        result.w = y.float_value;

    return result;
}
inline float4 hgpu_uint4_to_float4(uint4 x){
    float4 result;
    result.x = (float) x.x;
    result.y = (float) x.y;
    result.z = (float) x.z;
    result.w = (float) x.w;

    return result;
}

//________________________________________________________________________________________________________ XOR128 PRNG
__attribute__((always_inline)) uint4
xor128_step(uint4 seed)
{
    uint4 result;
    uint t = (seed.x^(seed.x<<11));

	result.x = seed.y;
	result.y = seed.z;
	result.z = seed.w;
	result.w = (seed.w^(seed.w>>19))^(t^(t>>8));

    return result;
}

__kernel void
xor128(__global uint4* seed_table, 
                     __global float4* randoms,
                     const uint N)
{
    uint giddst = GID;
    float4 result;
    float4 normal = (float4) XOR128_m_FP;
    uint4 seed = seed_table[GID];
    for (uint i = 0; i < N; i++) {
        seed = xor128_step(seed);
        result.x = (float) seed.x;
        seed = xor128_step(seed);
        result.y = (float) seed.x;
        seed = xor128_step(seed);
        result.z = (float) seed.x;
        seed = xor128_step(seed);
        result.w = (float) seed.x;
        randoms[giddst] = result / normal;
        giddst += GID_SIZE;
    }
    seed_table[GID] = seed;
}

//________________________________________________________________________________________________________ RANLUX PRNG
__attribute__((always_inline)) int4
rlseedint(int RL_jseed)
{
	int4 RL_output_int;
		int RL_k = RL_jseed / 53668;
		RL_jseed = 40014 * (RL_jseed - RL_k * 53668) - RL_k * 12211;
		if (RL_jseed < 0) {RL_jseed = RL_jseed + RL_icons;}
		RL_output_int.x = RL_jseed;

		RL_k = RL_jseed / 53668;
		RL_jseed = 40014 * (RL_jseed - RL_k * 53668) - RL_k * 12211;
		if (RL_jseed < 0) {RL_jseed = RL_jseed + RL_icons;}
		RL_output_int.y = RL_jseed;

		RL_k = RL_jseed / 53668;
		RL_jseed = 40014 * (RL_jseed - RL_k * 53668) - RL_k * 12211;
		if (RL_jseed < 0) {RL_jseed = RL_jseed + RL_icons;}
		RL_output_int.z = RL_jseed;

		RL_k = RL_jseed / 53668;
		RL_jseed = 40014 * (RL_jseed - RL_k * 53668) - RL_k * 12211;
		if (RL_jseed < 0) {RL_jseed = RL_jseed + RL_icons;}
		RL_output_int.w = RL_jseed;

    return RL_output_int;
}

__attribute__((always_inline)) float4
rlseedint2(int4 RL_temp)
{
	float4 RL_output_int;
	float4 RL_tmp = (float4) (RL_twom24,RL_twom24,RL_twom24,RL_twom24);
		RL_output_int.x = ((float) (RL_temp.x % RL_itwo24));
		RL_output_int.y = ((float) (RL_temp.y % RL_itwo24));
		RL_output_int.z = ((float) (RL_temp.z % RL_itwo24));
		RL_output_int.w = ((float) (RL_temp.w % RL_itwo24));
	RL_output_int = RL_output_int * RL_tmp;
    return RL_output_int;
}

__kernel void
rlseed(const __global uint * seeds, __global float4 * seedtable)
{
    int RL_jseed = (int) seeds[GID];
	int4		RL_temp1 = rlseedint(RL_jseed);
	seedtable[GID + 0 * GID_SIZE] = rlseedint2(RL_temp1);
				RL_temp1 = rlseedint(RL_temp1.w);
	seedtable[GID + 1 * GID_SIZE] = rlseedint2(RL_temp1);
				RL_temp1 = rlseedint(RL_temp1.w);
	seedtable[GID + 2 * GID_SIZE] = rlseedint2(RL_temp1);
				RL_temp1 = rlseedint(RL_temp1.w);
	seedtable[GID + 3 * GID_SIZE] = rlseedint2(RL_temp1);
				RL_temp1 = rlseedint(RL_temp1.w);
	seedtable[GID + 4 * GID_SIZE] = rlseedint2(RL_temp1);
				RL_temp1 = rlseedint(RL_temp1.w);
	float4		RL_temp2 = rlseedint2(RL_temp1);
	seedtable[GID + 5 * GID_SIZE] = RL_temp2;
	float4		RL_carin = (float4) (23.0f, 9.0f, 6.0f, 0.0f);
	if (RL_temp2.w == 0) {RL_carin.w = RL_twom24;}
	seedtable[GID + 6 * GID_SIZE] = RL_carin; // RL_i24, RL_j24, RL_in24, RL_carry
}

__attribute__((always_inline)) float4
rlproducefloat4(float4 * RL_seeds_i24, const float4 RL_seeds_j24, const float RL_seeds_j24P, float * RL_carry)
{
	float RL_c = *RL_carry;
	float4 RL_output_int = RL_seeds_j24 - (*RL_seeds_i24).wzyx;
	RL_output_int.x = RL_output_int.x - RL_c;
		if (RL_output_int.x < 0.0f) { RL_output_int.x = RL_output_int.x + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }
	RL_output_int.y = RL_output_int.y - RL_c;
		if (RL_output_int.y < 0.0f) { RL_output_int.y = RL_output_int.y + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }
	RL_output_int.z = RL_output_int.z - RL_c;
		if (RL_output_int.z < 0.0f) { RL_output_int.z = RL_output_int.z + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }
	RL_output_int.w = RL_output_int.w - RL_c;
		if (RL_output_int.w < 0.0f) { RL_output_int.w = RL_output_int.w + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }

	*RL_seeds_i24 = RL_output_int.wzyx;

    float4 RL_temporary = RL_seeds_j24.yzwx;
    RL_temporary.w = RL_seeds_j24P;

	float4 RL_output_int2 = select(RL_output_int, (RL_output_int + RL_twom24 * RL_temporary), RL_output_int < (float4) RL_twom12);
			RL_output_int = select(RL_output_int2, RL_twom24sq, RL_output_int2 == (float4) 0.0f);

	*RL_carry = RL_c;
	return RL_output_int;
}

__attribute__((always_inline)) void
rlproduceupd4(float4 * RL_seeds_i24, const float4 RL_seeds_j24, float * RL_carry)
{
	float RL_c = *RL_carry;
	float4 RL_output_int = RL_seeds_j24 - (*RL_seeds_i24).wzyx;
	RL_output_int.x = RL_output_int.x - RL_c;
		if (RL_output_int.x < 0.0f) { RL_output_int.x = RL_output_int.x + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }
	RL_output_int.y = RL_output_int.y - RL_c;
		if (RL_output_int.y < 0.0f) { RL_output_int.y = RL_output_int.y + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }
	RL_output_int.z = RL_output_int.z - RL_c;
		if (RL_output_int.z < 0.0f) { RL_output_int.z = RL_output_int.z + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }
	RL_output_int.w = RL_output_int.w - RL_c;
		if (RL_output_int.w < 0.0f) { RL_output_int.w = RL_output_int.w + 1.0f; RL_c = RL_twom24; } else { RL_c = 0.0f; }

	*RL_seeds_i24 = RL_output_int.wzyx;
	*RL_carry = RL_c;
}

__kernel void
rlproduce(__global float4 * seedtable,__global float4 * prns, const uint samples)
{
    uint giddst = GID;

	float4 RL_seeds_j24, uni, RL_seed_temp, RL_seed_temp2;
		float4 RL_seed0 = seedtable[GID + 0 * GID_SIZE];
		float4 RL_seed1 = seedtable[GID + 1 * GID_SIZE];
		float4 RL_seed2 = seedtable[GID + 2 * GID_SIZE];
		float4 RL_seed3 = seedtable[GID + 3 * GID_SIZE];
		float4 RL_seed4 = seedtable[GID + 4 * GID_SIZE];
		float4 RL_seed5 = seedtable[GID + 5 * GID_SIZE];
		float4 RL_carin = seedtable[GID + 6 * GID_SIZE];    // RL_i24, RL_j24, RL_in24, RL_carry
	float RL_carry = RL_carin.w;

	for (int i=0; i<samples; i++) {
		if (RL_carin.z >= 6.0f) {
			for (int t=0; t<RL_skip_24; t++) {
				RL_seeds_j24.xy = RL_seed2.yx;
				RL_seeds_j24.zw = RL_seed1.wz;
					rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed1.yx;
				RL_seeds_j24.zw = RL_seed0.wz;
					rlproduceupd4(&RL_seed4,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed0.yx;
				RL_seeds_j24.zw = RL_seed5.wz;
					rlproduceupd4(&RL_seed3,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed5.yx;
				RL_seeds_j24.zw = RL_seed4.wz;
					rlproduceupd4(&RL_seed2,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed4.yx;
				RL_seeds_j24.zw = RL_seed3.wz;
					rlproduceupd4(&RL_seed1,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed3.yx;
				RL_seeds_j24.zw = RL_seed2.wz;
					rlproduceupd4(&RL_seed0,RL_seeds_j24,&RL_carry);
			}
#if RL_skip_4 == 1
				RL_seeds_j24.xy = RL_seed2.yx;
				RL_seeds_j24.zw = RL_seed1.wz;
					rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
                RL_seed_temp = RL_seed5;    // 5->4, 4->3, 3->2, 2->1, 1->0, 0->5
                RL_seed5 = RL_seed4;
                RL_seed4 = RL_seed3;
                RL_seed3 = RL_seed2;
                RL_seed2 = RL_seed1;
                RL_seed1 = RL_seed0;
                RL_seed0 = RL_seed_temp;
#elif RL_skip_4 == 2
				RL_seeds_j24.xy = RL_seed2.yx;
				RL_seeds_j24.zw = RL_seed1.wz;
					rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed1.yx;
				RL_seeds_j24.zw = RL_seed0.wz;
					rlproduceupd4(&RL_seed4,RL_seeds_j24,&RL_carry);
                RL_seed_temp = RL_seed5;    // 5->3, 4->2, 3->1, 2->0, 1->5, 0->4
                RL_seed5 = RL_seed3;
                RL_seed3 = RL_seed1;
                RL_seed1 = RL_seed_temp;
                RL_seed_temp = RL_seed2;
                RL_seed2 = RL_seed0;
                RL_seed0 = RL_seed4;
                RL_seed4 = RL_seed_temp;
#elif RL_skip_4 == 3
				RL_seeds_j24.xy = RL_seed2.yx;
				RL_seeds_j24.zw = RL_seed1.wz;
					rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed1.yx;
				RL_seeds_j24.zw = RL_seed0.wz;
					rlproduceupd4(&RL_seed4,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed0.yx;
				RL_seeds_j24.zw = RL_seed5.wz;
					rlproduceupd4(&RL_seed3,RL_seeds_j24,&RL_carry);
                RL_seed_temp = RL_seed5;    // 5->2, 4->1, 3->0, 2->5, 1->4, 0->3
                RL_seed5 = RL_seed2;
                RL_seed2 = RL_seed_temp;
                RL_seed_temp = RL_seed4;
                RL_seed4 = RL_seed1;
                RL_seed1 = RL_seed_temp;
                RL_seed_temp = RL_seed3;
                RL_seed3 = RL_seed0;
                RL_seed0 = RL_seed_temp;
#elif RL_skip_4 == 4
				RL_seeds_j24.xy = RL_seed2.yx;
				RL_seeds_j24.zw = RL_seed1.wz;
					rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed1.yx;
				RL_seeds_j24.zw = RL_seed0.wz;
					rlproduceupd4(&RL_seed4,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed0.yx;
				RL_seeds_j24.zw = RL_seed5.wz;
					rlproduceupd4(&RL_seed3,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed5.yx;
				RL_seeds_j24.zw = RL_seed4.wz;
					rlproduceupd4(&RL_seed2,RL_seeds_j24,&RL_carry);
                RL_seed_temp = RL_seed5;    // 5->1, 4->0, 3->5, 2->4, 1->3, 0->2
                RL_seed5 = RL_seed1;
                RL_seed1 = RL_seed3;
                RL_seed3 = RL_seed_temp;
                RL_seed_temp = RL_seed4;
                RL_seed4 = RL_seed0;
                RL_seed0 = RL_seed2;
                RL_seed2 = RL_seed_temp;
#elif RL_skip_4 == 5
				RL_seeds_j24.xy = RL_seed2.yx;
				RL_seeds_j24.zw = RL_seed1.wz;
					rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed1.yx;
				RL_seeds_j24.zw = RL_seed0.wz;
					rlproduceupd4(&RL_seed4,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed0.yx;
				RL_seeds_j24.zw = RL_seed5.wz;
					rlproduceupd4(&RL_seed3,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed5.yx;
				RL_seeds_j24.zw = RL_seed4.wz;
					rlproduceupd4(&RL_seed2,RL_seeds_j24,&RL_carry);
				RL_seeds_j24.xy = RL_seed4.yx;
				RL_seeds_j24.zw = RL_seed3.wz;
					rlproduceupd4(&RL_seed1,RL_seeds_j24,&RL_carry);
                RL_seed_temp = RL_seed5;    // 5->0, 4->5, 3->4, 2->3, 1->2, 0->1
                RL_seed5 = RL_seed0;
                RL_seed0 = RL_seed1;
                RL_seed1 = RL_seed2;
                RL_seed2 = RL_seed3;
                RL_seed3 = RL_seed4;
                RL_seed4 = RL_seed_temp;
#endif
#if RL_skip_1 == 1
                float RL_output_int = RL_seed2.y - RL_seed5.w;
                RL_output_int = RL_output_int - RL_carry;
                    if (RL_output_int < 0.0f) { RL_output_int = RL_output_int + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_seed5.w = RL_output_int;

                RL_seed_temp = RL_seed5;
                RL_seed5.yzw = RL_seed5.xyz;
                RL_seed5.x = RL_seed4.w;
                RL_seed4.yzw = RL_seed4.xyz;
                RL_seed4.x = RL_seed3.w;
                RL_seed3.yzw = RL_seed3.xyz;
                RL_seed3.x = RL_seed2.w;
                RL_seed2.yzw = RL_seed2.xyz;
                RL_seed2.x = RL_seed1.w;
                RL_seed1.yzw = RL_seed1.xyz;
                RL_seed1.x = RL_seed0.w;
                RL_seed0.yzw = RL_seed0.xyz;
                RL_seed0.x = RL_seed_temp.w;
#elif RL_skip_1 == 2
                float4 RL_output_int;
                RL_output_int.xy = RL_seed2.yx - RL_seed5.wz;
                RL_output_int.x = RL_output_int.x - RL_carry;
                    if (RL_output_int.x < 0.0f) { RL_output_int.x = RL_output_int.x + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_output_int.y = RL_output_int.y - RL_carry;
                    if (RL_output_int.y < 0.0f) { RL_output_int.y = RL_output_int.y + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }

                RL_seed5.zw = RL_output_int.yx;

                RL_seed_temp = RL_seed5;
                RL_seed5.zw = RL_seed5.xy;
                RL_seed5.xy = RL_seed4.zw;
                RL_seed4.zw = RL_seed4.xy;
                RL_seed4.xy = RL_seed3.zw;
                RL_seed3.zw = RL_seed3.xy;
                RL_seed3.xy = RL_seed2.zw;
                RL_seed2.zw = RL_seed2.xy;
                RL_seed2.xy = RL_seed1.zw;
                RL_seed1.zw = RL_seed1.xy;
                RL_seed1.xy = RL_seed0.zw;
                RL_seed0.zw = RL_seed0.xy;
                RL_seed0.xy = RL_seed_temp.zw;
#elif RL_skip_1 == 3
                float4 RL_output_int;
                RL_output_int.xy = RL_seed2.yx - RL_seed5.wz;
                RL_output_int.z = RL_seed1.w - RL_seed5.y;
                RL_output_int.x = RL_output_int.x - RL_carry;
                    if (RL_output_int.x < 0.0f) { RL_output_int.x = RL_output_int.x + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_output_int.y = RL_output_int.y - RL_carry;
                    if (RL_output_int.y < 0.0f) { RL_output_int.y = RL_output_int.y + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_output_int.z = RL_output_int.z - RL_carry;
                    if (RL_output_int.z < 0.0f) { RL_output_int.z = RL_output_int.z + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }

                RL_seed5.yzw = RL_output_int.zyx;

                RL_seed_temp = RL_seed5;
                RL_seed5.w = RL_seed5.x;
                RL_seed5.xyz = RL_seed4.yzw;
                RL_seed4.w = RL_seed4.x;
                RL_seed4.xyz = RL_seed3.yzw;
                RL_seed3.w = RL_seed3.x;
                RL_seed3.xyz = RL_seed2.yzw;
                RL_seed2.w = RL_seed2.x;
                RL_seed2.xyz = RL_seed1.yzw;
                RL_seed1.w = RL_seed1.x;
                RL_seed1.xyz = RL_seed0.yzw;
                RL_seed0.w = RL_seed0.x;
                RL_seed0.xyz = RL_seed_temp.yzw;
#endif
#ifdef RL_RANLUX2
                float RL_output_int = RL_seed2.y - RL_seed5.w;
                RL_output_int = RL_output_int - RL_carry;
                    if (RL_output_int < 0.0f) { RL_output_int = RL_output_int + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_seed5.w = RL_output_int;

                RL_seed_temp = RL_seed5;
                RL_seed5.yzw = RL_seed5.xyz;
                RL_seed5.x = RL_seed4.w;
                RL_seed4.yzw = RL_seed4.xyz;
                RL_seed4.x = RL_seed3.w;
                RL_seed3.yzw = RL_seed3.xyz;
                RL_seed3.x = RL_seed2.w;
                RL_seed2.yzw = RL_seed2.xyz;
                RL_seed2.x = RL_seed1.w;
                RL_seed1.yzw = RL_seed1.xyz;
                RL_seed1.x = RL_seed0.w;
                RL_seed0.yzw = RL_seed0.xyz;
                RL_seed0.x = RL_seed_temp.w;
#endif
#ifdef RL_RANLUX3
                RL_seeds_j24.xy = RL_seed2.yx;
                RL_seeds_j24.zw = RL_seed1.wz;
                    rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
                float4 RL_output_int;
                RL_output_int.xy = RL_seed1.yx - RL_seed4.wz;
                RL_output_int.z = RL_seed0.w - RL_seed4.y;
                RL_output_int.x = RL_output_int.x - RL_carry;
                    if (RL_output_int.x < 0.0f) { RL_output_int.x = RL_output_int.x + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_output_int.y = RL_output_int.y - RL_carry;
                    if (RL_output_int.y < 0.0f) { RL_output_int.y = RL_output_int.y + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }
                RL_output_int.z = RL_output_int.z - RL_carry;
                    if (RL_output_int.z < 0.0f) { RL_output_int.z = RL_output_int.z + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }

                RL_seed_temp.xyz  = RL_output_int.zyx;
                RL_seed_temp.w    = RL_seed5.x;
                RL_seed_temp2.xyz = RL_seed5.yzw;
                RL_seed_temp2.w   = RL_seed0.x;

                RL_seed5.w   = RL_seed4.x;
                RL_seed5.xyz = RL_seed3.yzw;
                RL_seed4.w   = RL_seed3.x;
                RL_seed4.xyz = RL_seed2.yzw;
                RL_seed3.w   = RL_seed2.x;
                RL_seed3.xyz = RL_seed1.yzw;
                RL_seed2.w   = RL_seed1.x;
                RL_seed2.xyz = RL_seed0.yzw;
                RL_seed1     = RL_seed_temp2;
                RL_seed0     = RL_seed_temp;
#endif
#ifdef RL_RANLUX4
                RL_seeds_j24.xy = RL_seed2.yx;
                RL_seeds_j24.zw = RL_seed1.wz;
                    rlproduceupd4(&RL_seed5,RL_seeds_j24,&RL_carry);
                float RL_output_int = RL_seed1.y - RL_seed4.w;
                RL_output_int = RL_output_int - RL_carry;
                    if (RL_output_int < 0.0f) { RL_output_int = RL_output_int + 1.0f; RL_carry = RL_twom24; } else { RL_carry = 0.0f; }

                RL_seed_temp  = RL_seed5;
                RL_seed_temp2.x = RL_output_int;

                RL_seed5.yzw = RL_seed4.xyz;
                RL_seed5.x   = RL_seed3.w;
                RL_seed4.yzw = RL_seed3.xyz;
                RL_seed4.x   = RL_seed2.w;
                RL_seed3.yzw = RL_seed2.xyz;
                RL_seed3.x   = RL_seed1.w;
                RL_seed2.yzw = RL_seed1.xyz;
                RL_seed2.x   = RL_seed0.w;
                RL_seed1.yzw = RL_seed0.xyz;
                RL_seed1.x   = RL_seed_temp.w;
                RL_seed0.yzw = RL_seed_temp.xyz;
                RL_seed0.x   = RL_seed_temp2.x;
#endif
            RL_carin.z = 0.0f;
		}
        if (RL_carin.z < 1.0f) {
			RL_seeds_j24.xy = RL_seed2.yx;
			RL_seeds_j24.zw = RL_seed1.wz;
				uni = rlproducefloat4(&RL_seed5,RL_seeds_j24,RL_seed1.y,&RL_carry);
			RL_carin.z = 1.0f;
		} else if (RL_carin.z == 1.0f) {
			RL_seeds_j24.xy = RL_seed1.yx;
			RL_seeds_j24.zw = RL_seed0.wz;
				uni = rlproducefloat4(&RL_seed4,RL_seeds_j24,RL_seed0.y,&RL_carry);
			RL_carin.z = 2.0f;
		} else if (RL_carin.z == 2.0f) {
			RL_seeds_j24.xy = RL_seed0.yx;
			RL_seeds_j24.zw = RL_seed5.wz;
				uni = rlproducefloat4(&RL_seed3,RL_seeds_j24,RL_seed5.y,&RL_carry);
			RL_carin.z = 3.0f;
		} else if (RL_carin.z == 3.0f) {
			RL_seeds_j24.xy = RL_seed5.yx;
			RL_seeds_j24.zw = RL_seed4.wz;
				uni = rlproducefloat4(&RL_seed2,RL_seeds_j24,RL_seed4.y,&RL_carry);
			RL_carin.z = 4.0f;
		} else if (RL_carin.z == 4.0f) {
			RL_seeds_j24.xy = RL_seed4.yx;
			RL_seeds_j24.zw = RL_seed3.wz;
				uni = rlproducefloat4(&RL_seed1,RL_seeds_j24,RL_seed3.y,&RL_carry);
			RL_carin.z = 5.0f;
		} else if (RL_carin.z == 5.0f) {
			RL_seeds_j24.xy = RL_seed3.yx;
			RL_seeds_j24.zw = RL_seed2.wz;
				uni = rlproducefloat4(&RL_seed0,RL_seeds_j24,RL_seed2.y,&RL_carry);
			RL_carin.z = 6.0f;
		}
		prns[giddst] = uni;
		giddst += GID_SIZE;
	}

	seedtable[GID + 0 * GID_SIZE] = RL_seed0;
	seedtable[GID + 1 * GID_SIZE] = RL_seed1;
	seedtable[GID + 2 * GID_SIZE] = RL_seed2;
	seedtable[GID + 3 * GID_SIZE] = RL_seed3;
	seedtable[GID + 4 * GID_SIZE] = RL_seed4;
	seedtable[GID + 5 * GID_SIZE] = RL_seed5;

	RL_carin.w = RL_carry;
	seedtable[GID + 6 * GID_SIZE] = RL_carin;
}

//________________________________________________________________________________________________________ RANMAR PRNG
__kernel void
rmseed(const __global uint4 * seeds, __global float4 * seedtable)
{
    float4 seed;
    //
    uint4 RM_seed1 = seeds[GID];
    uint4 RM_seed2 = seeds[GID + GID_SIZE];
    //
    uint4 i = ((RM_seed1 / 177) % 177) + 2;
    uint4 j = (RM_seed1 % 177) + 2;
    uint4 k = ((RM_seed2 / 169) % 178) + 1;
    uint4 l = (RM_seed2 % 169);
    for (int n = 0; n < 97; n++)
    {
        uint4 s = (uint4) 0;
        uint  t = 8388608;
        for (int m = 0; m < 24;  m++)
        {
            uint4 u = (i * j) % 179;
            u = (u * k) % 179;
            i = j;
            j = k;
            k = u;
            l = (53 * l + 1) % 169;
            uint4 lus = (l * u) % 64;

            if ( lus.x >= 32) s.x += t;
            if ( lus.y >= 32) s.y += t;
            if ( lus.z >= 32) s.z += t;
            if ( lus.w >= 32) s.w += t;
            t = t >> 1;
        }
        seed.x = ((float) s.x);
        seed.y = ((float) s.y);
        seed.z = ((float) s.z);
        seed.w = ((float) s.w);
        seedtable[GID + GID_SIZE * n] = seed / ((float4) 16777216.0f);
	}
    Uint_and_Float indx_I97, indx_J97;
    indx_I97.uint_value = 96;
    indx_J97.uint_value = 32;
    float4 indx;
    indx.x = indx_I97.float_value;
    indx.y = indx_J97.float_value;
    indx.z = (362436.0f / 16777216.0f);
    indx.w = 0.0f;

    seedtable[GID + GID_SIZE * 97] = indx;   // (RM_I97,RM_J97,RM_C,0)
}

__kernel void
rmproduce(__global float4 * seedtable,__global float4 * prns, const uint samples)
{
    uint giddst = GID;
    float4 uni;
    Uint_and_Float indx_I97, indx_J97;
    float4 indx = seedtable[GID + GID_SIZE * 97];
    indx_I97.float_value = indx.x;
    indx_J97.float_value = indx.y;
    uint RM_I97 = indx_I97.uint_value;
    uint RM_J97 = indx_J97.uint_value;


    for (uint i=0; i<samples; i++) {
        uni = seedtable[GID + GID_SIZE * RM_I97] - seedtable[GID + GID_SIZE * RM_J97];
	    uni = select(uni, (uni + 1.0f), uni < (float4) 0.0f);
        seedtable[GID + GID_SIZE * RM_I97] = uni;

	    if (RM_I97 == 0) RM_I97 = 97;
	    if (RM_J97 == 0) RM_J97 = 97;
	    RM_I97--;
	    RM_J97--;

        indx.z -= RM_CD;
	    if (indx.z < 0.0f) {indx.z += RM_CM;}
	    uni.x -= indx.z;
	    uni.y -= indx.z;
	    uni.z -= indx.z;
	    uni.w -= indx.z;
	    uni = select(uni, (uni + 1.0f), uni < (float4) 0.0f);

        prns[giddst] = uni;
        giddst += GID_SIZE;
    }

    indx_I97.uint_value = RM_I97;
    indx_J97.uint_value = RM_J97;
    indx.x = indx_I97.float_value;
    indx.y = indx_J97.float_value;
    seedtable[GID + GID_SIZE * 97] = indx;
}

//________________________________________________________________________________________________________ PARK-MILLER PRNG
__kernel void
pm(__global uint4* seed_table, 
                     __global float4* randoms,
                     const uint N)
{
    uint giddst = GID;
    float4 normal = (float4) PM_m_FP;
    uint4 seed = seed_table[GID];
    for (uint i = 0; i < N; i++) {
        uint4 PM_hi = seed / PM_q;
        uint4 PM_lo = seed % PM_q;

        uint4 PM_test_1 = ((uint4) PM_a) * PM_lo;
        uint4 PM_test_2 = ((uint4) PM_r) * PM_hi;

        if (PM_test_1.x > PM_test_2.x) seed.x = PM_test_1.x - PM_test_2.x; else seed.x = PM_test_1.x - PM_test_2.x + (uint) PM_m;
        if (PM_test_1.y > PM_test_2.y) seed.y = PM_test_1.y - PM_test_2.y; else seed.y = PM_test_1.y - PM_test_2.y + (uint) PM_m;
        if (PM_test_1.z > PM_test_2.z) seed.z = PM_test_1.z - PM_test_2.z; else seed.z = PM_test_1.z - PM_test_2.z + (uint) PM_m;
        if (PM_test_1.w > PM_test_2.w) seed.w = PM_test_1.w - PM_test_2.w; else seed.w = PM_test_1.w - PM_test_2.w + (uint) PM_m;

        randoms[giddst] = hgpu_uint4_to_float4(seed) / normal;
        giddst += GID_SIZE;
    }
    seed_table[GID] = seed;
}

//________________________________________________________________________________________________________ XOR7 PRNG
__attribute__((always_inline)) void
xor7_step(uint4* seed1,uint4* seed2)
{
    uint t, y;
        t = (*seed2).w;     t = t ^ (t<<13);    y = t ^ (t<<9);
        t = (*seed2).x;     y^= t ^ (t<<7);
        t = (*seed1).w;     y^= t ^ (t>>3);
        t = (*seed1).y;     y^= t ^ (t>>10);
        t = (*seed1).x;     t = t ^ (t>>7);     y^= t ^ (t<<24);

        (*seed1).xyz = (*seed1).yzw;
        (*seed1).w   = (*seed2).x;
        (*seed2).xyz = (*seed2).yzw;
        (*seed2).w   = y;
}

__kernel void
xor7(__global uint4* seed_table, 
                     __global float4* randoms,
                     const uint N)
{
    uint giddst = GID;
    float4 result;
    uint4 seed1 = seed_table[GID];
    uint4 seed2 = seed_table[GID + GID_SIZE];
    for (uint i = 0; i < N; i++) {
        xor7_step(&seed1,&seed2);
            result.x = (float) seed2.w;
        xor7_step(&seed1,&seed2);
            result.y = (float) seed2.w;
        xor7_step(&seed1,&seed2);
            result.z = (float) seed2.w;
        xor7_step(&seed1,&seed2);
            result.w = (float) seed2.w;

        randoms[giddst] = result / ((float4) XOR7_m_FP);
        giddst += GID_SIZE;
    }
    seed_table[GID] = seed1;
    seed_table[GID + GID_SIZE] = seed2;
}

//________________________________________________________________________________________________________ RANECU
__kernel void
ranecu(__global uint4* seed_table, 
                     __global float4* randoms,
                     const uint N)
{
    uint giddst = GID;

    uint4 seed1 = seed_table[GID];
    uint4 seed2 = seed_table[GID + GID_SIZE];

    uint4 z;
    for (uint i = 0; i < N; i++) {

	    uint4 RANECU_k = seed1 / ((uint4) RANECU_seedP11);

        uint4 RANECU_test_1 = ((uint4) RANECU_seedP13) * (seed1 - RANECU_k * ((uint4) RANECU_seedP11));
        uint4 RANECU_test_2 = RANECU_k * ((uint4) RANECU_seedP12);

        if (RANECU_test_1.x > RANECU_test_2.x) seed1.x = RANECU_test_1.x - RANECU_test_2.x; else seed1.x = RANECU_test_1.x - RANECU_test_2.x + RANECU_icons1;
        if (RANECU_test_1.y > RANECU_test_2.y) seed1.y = RANECU_test_1.y - RANECU_test_2.y; else seed1.y = RANECU_test_1.y - RANECU_test_2.y + RANECU_icons1;
        if (RANECU_test_1.z > RANECU_test_2.z) seed1.z = RANECU_test_1.z - RANECU_test_2.z; else seed1.z = RANECU_test_1.z - RANECU_test_2.z + RANECU_icons1;
        if (RANECU_test_1.w > RANECU_test_2.w) seed1.w = RANECU_test_1.w - RANECU_test_2.w; else seed1.w = RANECU_test_1.w - RANECU_test_2.w + RANECU_icons1;

	    RANECU_k = seed2 / ((uint4) RANECU_seedP21);

        RANECU_test_1 = ((uint4) RANECU_seedP23) * (seed2 - RANECU_k * ((uint4) RANECU_seedP21));
        RANECU_test_2 = RANECU_k * ((uint4) RANECU_seedP22);

        if (RANECU_test_1.x > RANECU_test_2.x) seed2.x = RANECU_test_1.x - RANECU_test_2.x; else seed2.x = RANECU_test_1.x - RANECU_test_2.x + RANECU_icons2;
        if (RANECU_test_1.y > RANECU_test_2.y) seed2.y = RANECU_test_1.y - RANECU_test_2.y; else seed2.y = RANECU_test_1.y - RANECU_test_2.y + RANECU_icons2;
        if (RANECU_test_1.z > RANECU_test_2.z) seed2.z = RANECU_test_1.z - RANECU_test_2.z; else seed2.z = RANECU_test_1.z - RANECU_test_2.z + RANECU_icons2;
        if (RANECU_test_1.w > RANECU_test_2.w) seed2.w = RANECU_test_1.w - RANECU_test_2.w; else seed2.w = RANECU_test_1.w - RANECU_test_2.w + RANECU_icons2;

        if (seed1.x > seed2.x) z.x = seed1.x - seed2.x; else z.x = seed1.x - seed2.x + RANECU_icons3;
        if (seed1.y > seed2.y) z.y = seed1.y - seed2.y; else z.y = seed1.y - seed2.y + RANECU_icons3;
        if (seed1.z > seed2.z) z.z = seed1.z - seed2.z; else z.z = seed1.z - seed2.z + RANECU_icons3;
        if (seed1.w > seed2.w) z.w = seed1.w - seed2.w; else z.w = seed1.w - seed2.w + RANECU_icons3;

        randoms[giddst] = hgpu_uint4_to_float4(z) / ((float4) RANECU_twom31);
        giddst += GID_SIZE;
    }
    seed_table[GID] = seed1;
    seed_table[GID + GID_SIZE] = seed2;
}


#undef RANECU_seedP23
#undef RANECU_seedP22
#undef RANECU_seedP21
#undef RANECU_seedP13
#undef RANECU_seedP12
#undef RANECU_seedP11
#undef RANECU_icons3
#undef RANECU_icons2
#undef RANECU_icons1
#undef RANECU_twom31

#undef XOR7_m_FP

#undef PM_m_FP
#undef PM_m
#undef PM_a
#undef PM_q
#undef PM_r

#undef PM_m_FP
#undef PM_m
#undef PM_a
#undef PM_q
#undef PM_r

#undef RM_CD
#undef RM_CM

#undef RL_RANLUX2
#undef RL_RANLUX3
#undef RL_RANLUX4
#undef RL_zero
#undef RL_icons
#undef RL_itwo24
#undef RL_twom24
#undef RL_twom12
#undef RL_twom24sq
#undef RL_skip
#undef RL_skip_1
#undef RL_skip_4
#undef RL_skip_24

#undef XOR128_m_FP

#endif
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                     
                                                                                                                                                      
                                                                                                                                                      
