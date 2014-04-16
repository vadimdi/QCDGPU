/******************************************************************************
 * @file     random.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Pseudo-random numbers generators library
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

#include "random.h"
#include <math.h>

namespace PRNG_CL{
using PRNG_CL::PRNG;

#define SOURCE_PRNG         "random/random.cl"

#define FNAME_MAX_LENGTH     250  // max length of filename with path


                    PRNG::PRNG(void)	
{
    random = NULL;

    seeds_size               = 0;    // size of input seeds table
    seed_table_size          = 0;    // size of seed table
    randoms_size             = 0;    // size of output buffer for randoms
    randoms_produced         = 0;    // whole number of produced numbers (from first produced number)
    pointer_to_randoms       = NULL;

    PRNG_generator  = PRNG_generator_none;
    PRNG_instances  = 0;     // number of PRNG instances
    PRNG_samples    = 500;   // number of prn samples per one generator
    PRNG_srandtime  = 0;
    PRNG_randseries = 0;     // Type of random series (0: time-dependent series, #_any_#: constant series for #_any_#) 

    PRNG_counter      = 0;   // counter runs of subroutine PRNG_produce 

    PRNG_seeds               = NULL;
    PRNG_seed_id             = 0;    // input seeds ID
    PRNG_seed_table_id       = 0;    // seed table ID
    PRNG_randoms_id          = 0;    // output buffer for randoms ID
    PRNG_seeds4              = NULL;
    PRNG_seed_table_uint4    = NULL;
    PRNG_seed_table_float4   = NULL;
    PRNG_randoms             = NULL;
    PRNG_randoms_kernel_id   = 0;    // kernel ID
    PRNG_randoms_seed_id     = 0;    // kernel for seed ID

    // RANLUX parameters___________________________________
        // Planar values______________________________________
        //static const int	RL_nskip	= 384;	// lux level 4
        //static const int	RL_nskip	= 224;	// lux level 3
        //static const int	RL_nskip	= 96;	// lux level 2
        //static const int	RL_nskip	= 24;	// lux level 1
        //static const int	RL_nskip	= 0;	// lux level 0

        //  0  1   2   3   4  LUXURY level
        // 24 48  97 223 389  original
        // 24 48 120 240 408  planar
        // 24 48 100 224 404
    RL_seed          = 314159265;    // MAX = 2147483647
    RL_nskip         = 24;
    RL_i24           = 23;
    RL_j24           = 9;
    RL_in24          = 0;
    RL_carry         = 0.0f;
    RL_jseed         = RL_seed;

    // RANMAR parameters___________________________________
    RM_I97 = 96;
    RM_J97 = 32;
    RM_seed1 = 1802;
    RM_seed2 = 9373;
    RM_C = (362436.0 / 16777216.0);

    // Park-Miller (GGL) parameters________________________
    PMseedFP = (float) PM_seed;
    PMseed   = PM_seed;

    // XORSeven parameters_________________________________
	for (int i=0; i<8; i++) XOR7_state[i] = 0;
	XOR7_index = 0;

    // XOR128 parameters___________________________________
	XOR128_state.s[0] = 0;
	XOR128_state.s[1] = 0;
	XOR128_state.s[2] = 0;
	XOR128_state.s[3] = 0;

    // RANECU parameters___________________________________
    RANECU_jseed1 = RANECU_seed1;
    RANECU_jseed2 = RANECU_seed2;
}
                    PRNG::~PRNG(void)	
{

}
	
double              PRNG::trunc(double x)
{
	double y, z;
	y = modf(x,&z);
	return z;
}

unsigned int PRNG::convert_generator_to_uint(PRNG::PRNG_generators generator){
    if (generator == PRNG::PRNG_generator_XOR128)  return  1;
    if (generator == PRNG::PRNG_generator_RANLUX0) return  2;
    if (generator == PRNG::PRNG_generator_RANLUX1) return  3;
    if (generator == PRNG::PRNG_generator_RANLUX2) return  4;
    if (generator == PRNG::PRNG_generator_RANLUX3) return  5;
    if (generator == PRNG::PRNG_generator_RANLUX4) return  6;
    if (generator == PRNG::PRNG_generator_RANLUX)  return  7;
    if (generator == PRNG::PRNG_generator_RANMAR)  return  8;
    if (generator == PRNG::PRNG_generator_PM)      return  9;
    if (generator == PRNG::PRNG_generator_XOR7)    return 10;
    if (generator == PRNG::PRNG_generator_RANECU)  return 11;
    return 5; // return RANLUX3 generator otherwise
}
PRNG::PRNG_generators PRNG::convert_uint_to_generator(unsigned int generator){
    if (generator == 1) return PRNG::PRNG_generator_XOR128;
    if (generator == 2) return PRNG::PRNG_generator_RANLUX0;
    if (generator == 3) return PRNG::PRNG_generator_RANLUX1;
    if (generator == 4) return PRNG::PRNG_generator_RANLUX2;
    if (generator == 5) return PRNG::PRNG_generator_RANLUX3;
    if (generator == 6) return PRNG::PRNG_generator_RANLUX4;
    if (generator == 7) return PRNG::PRNG_generator_RANLUX;
    if (generator == 8) return PRNG::PRNG_generator_RANMAR;
    if (generator == 9) return PRNG::PRNG_generator_PM;
    if (generator ==10) return PRNG::PRNG_generator_XOR7;
    if (generator ==11) return PRNG::PRNG_generator_RANECU;
    return PRNG::PRNG_generator_RANLUX3; // return RANLUX3 generator otherwise
}


void                PRNG::parameters_setup(char* parameter,int ivalue,char* text_value){
            if (!strcmp(parameter,"RANDSERIES"))   PRNG_randseries = ivalue;
            if (!strcmp(parameter,"PRNG"))  {
                if (!strcmp(text_value,"RANLUX0")) PRNG_generator   = PRNG_generator_RANLUX0;
                if (!strcmp(text_value,"RANLUX1")) PRNG_generator   = PRNG_generator_RANLUX1;
                if (!strcmp(text_value,"RANLUX2")) PRNG_generator   = PRNG_generator_RANLUX2;
                if (!strcmp(text_value,"RANLUX3")) PRNG_generator   = PRNG_generator_RANLUX3;
                if (!strcmp(text_value,"RANLUX4")) PRNG_generator   = PRNG_generator_RANLUX4;
                if (!strcmp(text_value,"RANLUX"))  PRNG_generator   = PRNG_generator_RANLUX;
                if (!strcmp(text_value,"XOR128"))  PRNG_generator   = PRNG_generator_XOR128;
                if (!strcmp(text_value,"XOR7"))    PRNG_generator   = PRNG_generator_XOR7;
                if (!strcmp(text_value,"RANMAR"))  PRNG_generator   = PRNG_generator_RANMAR;
                if (!strcmp(text_value,"RANECU"))  PRNG_generator   = PRNG_generator_RANECU;
                if (!strcmp(text_value,"PM"))      PRNG_generator   = PRNG_generator_PM;
            }
}
int                 PRNG::print_generator(char* header,int header_size){
    int j = 0;

    j  += sprintf_s(header+j,header_size-j, " random series               : %u\n",PRNG_srandtime); // PRNG_randseries
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANLUX)      j  += sprintf_s(header+j,header_size-j, " PRN generator (%3u)         : RANLUX\n",RL_nskip);
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANLUX0)     j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANLUX0\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANLUX1)     j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANLUX1\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANLUX2)     j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANLUX2\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANLUX3)     j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANLUX3\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANLUX4)     j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANLUX4\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_XOR128)      j  += sprintf_s(header+j,header_size-j, " PRN generator               : XOR128\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANECU)      j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANECU\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_RANMAR)      j  += sprintf_s(header+j,header_size-j, " PRN generator               : RANMAR\n");
    if (PRNG_generator == PRNG_CL::PRNG::PRNG_generator_XOR7)        j  += sprintf_s(header+j,header_size-j, " PRN generator               : XOR7\n");

    return j;
}

void                PRNG::initialize(void)
{
        bool PRNG_instances_autoselect = false;
        if (PRNG_instances==0) {
            PRNG_instances_autoselect = true;
            PRNG_instances = (int) GPU0->GPU_info.max_memory_width;
            PRNG_samples = 1 + (PRNG_samples - 1) / PRNG_instances;
        }
        printf("PRNG instances: %u",PRNG_instances);
        if (PRNG_instances_autoselect) printf(" (auto)");
        printf("\n");
    	// setup initial state: randseries=0(random, system timer-dependent) or other(static)
    	PRNG_srandtime = (unsigned int) time(NULL);
		if (!PRNG_randseries) {
            PRNG_randseries = PRNG_srandtime;  // set the RANDOM seed random
        } else {
            PRNG_srandtime = PRNG_randseries;
        }

        srand(PRNG_srandtime);

        if (PRNG_generator==PRNG_generator_XOR128){
            // XOR128
            XOR128_initialize_CPU();
            XOR128_initialize();
        }
        if ((PRNG_generator==PRNG_generator_RANLUX0)||(PRNG_generator==PRNG_generator_RANLUX1)||
            (PRNG_generator==PRNG_generator_RANLUX2)||(PRNG_generator==PRNG_generator_RANLUX3)||
            (PRNG_generator==PRNG_generator_RANLUX4)||(PRNG_generator==PRNG_generator_RANLUX)){
            // RANLUX
            if (PRNG_generator==PRNG_generator_RANLUX0) RL_nskip = 24;
            if (PRNG_generator==PRNG_generator_RANLUX1) RL_nskip = 48;
            if (PRNG_generator==PRNG_generator_RANLUX2) RL_nskip = 97;
            if (PRNG_generator==PRNG_generator_RANLUX3) RL_nskip = 223;//223;
            if (PRNG_generator==PRNG_generator_RANLUX4) RL_nskip = 389;
            RL_initialize_CPU();
            RL_initialize();
        }
        if (PRNG_generator==PRNG_generator_RANMAR){
            // RANMAR
            RANMAR_initialize_CPU();
            RANMAR_initialize();
        }
        if (PRNG_generator==PRNG_generator_PM){
            // Park-Miller
            PM_initialize_CPU();
            PM_initialize();
        }

        if (PRNG_generator==PRNG_generator_XOR7){
            // XOR7
            XOR7_initialize_CPU();
            XOR7_initialize();
        }

        if (PRNG_generator==PRNG_generator_RANECU){
            // RANECU
            RANECU_initialize_CPU();
            RANECU_initialize();
        }
        produce();
}
void                PRNG::produce(void)
{
        GPU0->kernel_run(PRNG_randoms_kernel_id);
        randoms_produced += PRNG_samples * 4;
        PRNG_counter++;
}
void                PRNG::produce_CPU(float* randoms_cpu)
{
    if (PRNG_generator==PRNG_generator_XOR128) XOR128_produce_CPU(randoms_cpu);
    if ((PRNG_generator==PRNG_generator_RANLUX0)||(PRNG_generator==PRNG_generator_RANLUX1)||
        (PRNG_generator==PRNG_generator_RANLUX2)||(PRNG_generator==PRNG_generator_RANLUX3)||
        (PRNG_generator==PRNG_generator_RANLUX4)||(PRNG_generator==PRNG_generator_RANLUX)){
            RL_produce_CPU(randoms_cpu);
    }
    if (PRNG_generator==PRNG_generator_RANMAR) RANMAR_produce_CPU(randoms_cpu);
    if (PRNG_generator==PRNG_generator_PM)     PM_produce_CPU(randoms_cpu);
    if (PRNG_generator==PRNG_generator_XOR7)   XOR7_produce_CPU(randoms_cpu);
    if (PRNG_generator==PRNG_generator_RANECU) RANECU_produce_CPU(randoms_cpu);
}
void                PRNG::produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    if (PRNG_generator==PRNG_generator_XOR128) XOR128_produce_CPU(randoms_cpu,number_of_prns_CPU);
    if ((PRNG_generator==PRNG_generator_RANLUX0)||(PRNG_generator==PRNG_generator_RANLUX1)||
        (PRNG_generator==PRNG_generator_RANLUX2)||(PRNG_generator==PRNG_generator_RANLUX3)||
        (PRNG_generator==PRNG_generator_RANLUX4)||(PRNG_generator==PRNG_generator_RANLUX)){
            RL_produce_CPU(randoms_cpu,number_of_prns_CPU);
    }
    if (PRNG_generator==PRNG_generator_RANMAR) RANMAR_produce_CPU(randoms_cpu,number_of_prns_CPU);
    if (PRNG_generator==PRNG_generator_PM)     PM_produce_CPU(randoms_cpu,number_of_prns_CPU);
    if (PRNG_generator==PRNG_generator_XOR7)   XOR7_produce_CPU(randoms_cpu,number_of_prns_CPU);
    if (PRNG_generator==PRNG_generator_RANECU) RANECU_produce_CPU(randoms_cpu,number_of_prns_CPU);
}
unsigned int        PRNG::check_seeds(void)
{
    // check seeds
    unsigned int result = 0;
    if ((PRNG_generator==PRNG_generator_RANLUX0)||(PRNG_generator==PRNG_generator_RANLUX1)||
    (PRNG_generator==PRNG_generator_RANLUX2)||(PRNG_generator==PRNG_generator_RANLUX3)||
    (PRNG_generator==PRNG_generator_RANLUX4)||(PRNG_generator==PRNG_generator_RANLUX)){
        result = RL_check_seeds();
        if (result) RL_print_seeds();
    }
    if (PRNG_generator==PRNG_generator_RANMAR){
        result = RANMAR_check_seeds();
        if (result) RANMAR_print_seeds();
    }

    if (result) printf("Seeds check failed: %u errors\n",result);
    else printf("Seeds check passed!\n");

    return result;
}
unsigned int        PRNG::check(void)
{

//    PRNG_randseries and PRNG_samples must be initialized in the parent program
    unsigned int result = 0;
    produce();

    // CPU prn production
    float* randoms_cpu = (float*) calloc(PRNG_samples*4,sizeof(float));
    produce_CPU(randoms_cpu);

    // check results
    cl_float4* pointer_to_randoms = GPU0->buffer_map_float4(PRNG_randoms_id);

    // recheck seeds
    check_seeds();

    int max_output = 512;
    unsigned int offset = PRNG_instances;
    int i_quads = 4;
    if (PRNG_generator == PRNG_generator_RANMAR) {i_quads = 1;}
    if (PRNG_generator == PRNG_generator_PM)     {i_quads = 1;}
    if (PRNG_generator == PRNG_generator_RANECU) {i_quads = 1;}
    for (int i=0; i<PRNG_samples; i++) {
        double prn_cpu0 = 0.0;
        double prn_cpu1 = 0.0;
        double prn_cpu2 = 0.0;
        double prn_cpu3 = 0.0;
        if (i_quads == 4) {
            prn_cpu0 = randoms_cpu[4*i  ];
            prn_cpu1 = randoms_cpu[4*i+1];
            prn_cpu2 = randoms_cpu[4*i+2];
            prn_cpu3 = randoms_cpu[4*i+3];
        } else {
            prn_cpu0 = randoms_cpu[i];
        }

        unsigned int prn_gpu_int0 = GPU0->convert_to_uint((float) pointer_to_randoms[offset*i].s[0]);
        unsigned int prn_gpu_int1 = GPU0->convert_to_uint((float) pointer_to_randoms[offset*i].s[1]);
        unsigned int prn_gpu_int2 = GPU0->convert_to_uint((float) pointer_to_randoms[offset*i].s[2]);
        unsigned int prn_gpu_int3 = GPU0->convert_to_uint((float) pointer_to_randoms[offset*i].s[3]);
        double prn_gpu0 = GPU0->convert_to_double(pointer_to_randoms[offset*i].s[0]);
        double prn_gpu1 = GPU0->convert_to_double(pointer_to_randoms[offset*i].s[1]);
        double prn_gpu2 = GPU0->convert_to_double(pointer_to_randoms[offset*i].s[2]);
        double prn_gpu3 = GPU0->convert_to_double(pointer_to_randoms[offset*i].s[3]);

        if ((i_quads >= 1) && (prn_cpu0-prn_gpu0)) result++;
        if ((i_quads >= 2) && (prn_cpu1-prn_gpu1)) result++;
        if ((i_quads >= 3) && (prn_cpu2-prn_gpu2)) result++;
        if ((i_quads >= 4) && (prn_cpu3-prn_gpu3)) result++;

        if (max_output>0) {
            max_output--;
            if ((i_quads >= 1) && (prn_cpu0-prn_gpu0)) printf("[%2u].x: % e\t (CPU:% e,\t GPU:% e,\t intGPU:%u)\n",i,(prn_cpu0-prn_gpu0),prn_cpu0,prn_gpu0,prn_gpu_int0);
            if ((i_quads >= 2) && (prn_cpu1-prn_gpu1)) printf("[%2u].y: % e\t (CPU:% e,\t GPU:% e,\t intGPU:%u)\n",i,(prn_cpu1-prn_gpu1),prn_cpu1,prn_gpu1,prn_gpu_int1);
            if ((i_quads >= 3) && (prn_cpu2-prn_gpu2)) printf("[%2u].z: % e\t (CPU:% e,\t GPU:% e,\t intGPU:%u)\n",i,(prn_cpu2-prn_gpu2),prn_cpu2,prn_gpu2,prn_gpu_int2);
            if ((i_quads >= 4) && (prn_cpu3-prn_gpu3)) printf("[%2u].w: % e\t (CPU:% e,\t GPU:% e,\t intGPU:%u)\n",i,(prn_cpu3-prn_gpu3),prn_cpu3,prn_gpu3,prn_gpu_int3);
        }
    }

    delete[] randoms_cpu;

    printf("PRNG check result: %u errors (%u samples)\n\n",result,PRNG_samples*4);

    return result;
}
unsigned int        PRNG::check_range(void)
{
    unsigned int result = 0;

    // check results
    cl_float4* pointer_to_randoms = GPU0->buffer_map_float4(PRNG_randoms_id);
    for (unsigned int i=0; i<randoms_size; i++) {
        double prn_gpu0 = GPU0->convert_to_double(pointer_to_randoms[i].s[0]);
        double prn_gpu1 = GPU0->convert_to_double(pointer_to_randoms[i].s[1]);
        double prn_gpu2 = GPU0->convert_to_double(pointer_to_randoms[i].s[2]);
        double prn_gpu3 = GPU0->convert_to_double(pointer_to_randoms[i].s[3]);

        if ((prn_gpu0 <= 0.0) || (prn_gpu0 >= 1.0) || (prn_gpu0 != prn_gpu0)) result++;
        if ((prn_gpu1 <= 0.0) || (prn_gpu1 >= 1.0) || (prn_gpu1 != prn_gpu1)) result++;
        if ((prn_gpu2 <= 0.0) || (prn_gpu2 >= 1.0) || (prn_gpu2 != prn_gpu2)) result++;
        if ((prn_gpu3 <= 0.0) || (prn_gpu3 >= 1.0) || (prn_gpu3 != prn_gpu3)) result++;
    delete(pointer_to_randoms);
    }

    return result;
}

void                PRNG::RL_initialize(void)
{
        char options[1024];
        sprintf_s(options,sizeof(options),"-D RL_skip=%i",RL_nskip);
        char buffer_prng_cl[FNAME_MAX_LENGTH];
        int  j = sprintf_s(buffer_prng_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
             j+= sprintf_s(buffer_prng_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_PRNG);
        random = GPU0->source_read(buffer_prng_cl);
                 GPU0->program_create(random,options);

        seeds_size      = PRNG_instances;
	    seed_table_size = seeds_size * 7;               // 7 = size of Ranlux seed table for each PRNG (in quads)
        randoms_size    = PRNG_instances * PRNG_samples;

	    PRNG_seeds              = (cl_uint*)   calloc(seeds_size,     sizeof(cl_uint));             // Input seeds
        PRNG_randoms            = (cl_float4*) calloc(randoms_size,   sizeof(cl_float4));           // Output randoms
        PRNG_seed_table_float4  = (cl_float4*) calloc(seed_table_size,sizeof(cl_float4));           // Seeds

        if (GPU0->GPU_debug.brief_report){
            printf("seeds_size                  = %u\n",seeds_size);
            printf("randoms_size                = %u\n",randoms_size);
            printf("seed_table_size             = %u\n",seed_table_size);
            printf("PRNG_samples                = %u\n",PRNG_samples);
            printf("PRNG_instances              = %u\n",PRNG_instances);
            printf("PRNG will produce           = %u\n",PRNG_instances * PRNG_samples);
        }

        PRNG_seeds[0] = RL_jseed;                       // setup first thread as CPU

        for (unsigned int i=1; i<seeds_size; i++)
          PRNG_seeds[i] = rand();

        PRNG_seed_id       = GPU0->buffer_init(GPU0->buffer_type_Constant,  seeds_size,      PRNG_seeds,             sizeof(cl_uint));
    GPU0->print_stage("seeds initialized");
        PRNG_seed_table_id = GPU0->buffer_init(GPU0->buffer_type_IO,        seed_table_size, PRNG_seed_table_float4, sizeof(cl_float4));
    GPU0->print_stage("seed_table initialized");
        PRNG_randoms_id    = GPU0->buffer_init(GPU0->buffer_type_IO,        randoms_size,    PRNG_randoms,           sizeof(cl_float4));
    GPU0->print_stage("randoms initialized");

        int argument_id;

        const size_t global_size[]  = {PRNG_instances};

	    PRNG_randoms_seed_id = GPU0->kernel_init("rlseed",1,global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(PRNG_randoms_seed_id,  PRNG_seed_id);
            argument_id = GPU0->kernel_init_buffer(PRNG_randoms_seed_id,  PRNG_seed_table_id);

	    int result = GPU0->kernel_run(PRNG_randoms_seed_id);                // Prepare seeds

        PRNG_randoms_kernel_id = GPU0->kernel_init("rlproduce",1,global_size,NULL);
	        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_seed_table_id);
		    argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_randoms_id);
		    argument_id = GPU0->kernel_init_constant(PRNG_randoms_kernel_id,&PRNG_samples);

        result |= GPU0->buffer_kill(PRNG_seed_id);
}
void                PRNG::RL_initialize_CPU(void)
{
        RL_jseed = rand() % 2147483647;

	    int		RL_k;
	    
	    int RL_jseed_CPU = RL_jseed;
	    
	    for (int i=0; i<24; i++)
	    {
		    RL_k = RL_jseed_CPU / 53668;
		    RL_jseed_CPU = 40014 * (RL_jseed_CPU - RL_k * 53668) - RL_k * 12211;
		    if (RL_jseed_CPU < 0) {RL_jseed_CPU = RL_jseed_CPU + RL_icons;}
		    RL_seeds[i] = ((float) (RL_jseed_CPU % RL_itwo24)) * RL_twom24;
	    }
	    if (RL_seeds[23] == 0) RL_carry = RL_twom24;
}
float               PRNG::RL_produce_one_CPU(void)
{
	float uni;

	if (RL_in24 == 0)
	{
		RL_in24 = 24;
		for (int isk = 0; isk < RL_nskip; isk++)	// TEMPORARY!!!
		{
			uni = RL_seeds[RL_j24] - RL_seeds[RL_i24] - RL_carry;
			if (uni < 0.0f)
			{
				uni = uni + 1.0f;
				RL_carry = RL_twom24;
			}
			else
			{
				RL_carry = 0.0f;
			}
			RL_seeds[RL_i24] = uni;
			if (RL_i24==0) {RL_i24=24;}
			if (RL_j24==0) {RL_j24=24;}
			RL_i24--;
			RL_j24--;
		}
	}

			uni = RL_seeds[RL_j24] - RL_seeds[RL_i24] - RL_carry;
			if (uni < 0.0f)
		{
			uni = uni + 1.0f;
			RL_carry = RL_twom24;
		}
			else
		{
			RL_carry = 0.0f;
		}
			RL_seeds[RL_i24] = uni;
			if (RL_i24==0) {RL_i24=24;}
			if (RL_j24==0) {RL_j24=24;}
			RL_i24--;
			RL_j24--;

	if (uni < RL_twom12)
	{
		uni = uni + RL_twom24 * RL_seeds[RL_j24];
		if (uni == 0.0f) {uni = RL_twom24 * RL_twom24;}
	}
	RL_in24--;

	return uni;
}
unsigned int        PRNG::RL_check_seeds(void)
{
    // check seeds
    unsigned int result = 0;
    unsigned int* seedz = GPU0->buffer_map(PRNG_seed_table_id);
    int index = RL_get_seed_table_index(RL_nskip,randoms_produced);
    for (int i=0; i<6; i++){
        if (GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+0])-RL_seeds[index]) result++;
            index++; if (index >= 24) index = 0;
        if (GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+1])-RL_seeds[index]) result++;
            index++; if (index >= 24) index = 0;
        if (GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+2])-RL_seeds[index]) result++;
            index++; if (index >= 24) index = 0;
        if (GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+3])-RL_seeds[index]) result++;
            index++; if (index >= 24) index = 0;
    }
    return result;
}
void                PRNG::RL_print_seeds(void)
{
    // check seeds
    unsigned int* seedz = GPU0->buffer_map(PRNG_seed_table_id);
    int index = RL_get_seed_table_index(RL_nskip,randoms_produced);
    for (int i=0; i<6; i++){
        double cpu_seed0 = RL_seeds[index];
            index++; if (index >= 24) index = 0;
        double cpu_seed1 = RL_seeds[index];
            index++; if (index >= 24) index = 0;
        double cpu_seed2 = RL_seeds[index];
            index++; if (index >= 24) index = 0;
        double cpu_seed3 = RL_seeds[index];
            index++; if (index >= 24) index = 0;
        double gpu_seed0 = GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+0]);
        double gpu_seed1 = GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+1]);
        double gpu_seed2 = GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+2]);
        double gpu_seed3 = GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+3]);
        printf("\n[%2u]: (CPU)=%e (GPU)=%e",4*i+0,cpu_seed0,gpu_seed0); if (cpu_seed0!=gpu_seed0) printf(" !=");
        printf("\n[%2u]: (CPU)=%e (GPU)=%e",4*i+1,cpu_seed1,gpu_seed1); if (cpu_seed1!=gpu_seed1) printf(" !=");
        printf("\n[%2u]: (CPU)=%e (GPU)=%e",4*i+2,cpu_seed2,gpu_seed2); if (cpu_seed2!=gpu_seed2) printf(" !=");
        printf("\n[%2u]: (CPU)=%e (GPU)=%e",4*i+3,cpu_seed3,gpu_seed3); if (cpu_seed3!=gpu_seed3) printf(" !=");
    }
    printf("\n");
}
int                 PRNG::RL_get_seed_table_index(int skip,int produced){
    int index = 0;
        index = ((24 - (skip * ((int) (20 + produced) / 24)) % 24) % 24);
    return index;
}
void                PRNG::RL_produce_CPU(float* randoms_cpu)
{
    int number_of_prns_CPU = PRNG_samples * 4;
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = RL_produce_one_CPU();
}
void                PRNG::RL_produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = RL_produce_one_CPU();
}

void                PRNG::XOR128_initialize(void)
{
        char buffer_prng_cl[FNAME_MAX_LENGTH];
        int  j = sprintf_s(buffer_prng_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
             j+= sprintf_s(buffer_prng_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_PRNG);
        random = GPU0->source_read(buffer_prng_cl);
                 GPU0->program_create(random);

        seed_table_size         = GPU0->buffer_size_align(PRNG_instances);
        randoms_size            = GPU0->buffer_size_align(PRNG_instances * PRNG_samples);
        PRNG_seed_table_uint4   = (cl_uint4*)  calloc(seed_table_size,sizeof(cl_uint4));
        PRNG_randoms            = (cl_float4*) calloc(randoms_size,sizeof(cl_float4));

        PRNG_seed_table_uint4[0].s[0] = XOR128_state.s[0];    // setup first thread as CPU
        PRNG_seed_table_uint4[0].s[1] = XOR128_state.s[1];
        PRNG_seed_table_uint4[0].s[2] = XOR128_state.s[2];
        PRNG_seed_table_uint4[0].s[3] = XOR128_state.s[3];

        for (unsigned int i=1; i<seed_table_size; i++)
        {
          PRNG_seed_table_uint4[i].s[0] = rand();
          PRNG_seed_table_uint4[i].s[1] = rand();
          PRNG_seed_table_uint4[i].s[2] = rand();
          PRNG_seed_table_uint4[i].s[3] = rand();
        }

        PRNG_seed_table_id = GPU0->buffer_init(GPU0->buffer_type_IO, seed_table_size, PRNG_seed_table_uint4,    sizeof(cl_uint4));
        PRNG_randoms_id    = GPU0->buffer_init(GPU0->buffer_type_IO, randoms_size,    PRNG_randoms,             sizeof(cl_float4));

        int argument_id;
        const size_t global_size[]  = {PRNG_instances};

        PRNG_randoms_kernel_id   = GPU0->kernel_init("xor128",1,global_size,NULL);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_seed_table_id);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_randoms_id);
        argument_id = GPU0->kernel_init_constant(PRNG_randoms_kernel_id,&PRNG_samples);
}
void                PRNG::XOR128_initialize_CPU(void)
{
        XOR128_state.s[0] = rand();
        XOR128_state.s[1] = rand();
        XOR128_state.s[2] = rand();
        XOR128_state.s[3] = rand();
}
unsigned int        PRNG::XOR128_produce_one_uint_CPU(void)
{
        unsigned long t=(XOR128_state.s[0]^(XOR128_state.s[0]<<11));

	    XOR128_state.s[0] = XOR128_state.s[1];
	    XOR128_state.s[1] = XOR128_state.s[2];
	    XOR128_state.s[2] = XOR128_state.s[3];
	    XOR128_state.s[3] = (XOR128_state.s[3]^(XOR128_state.s[3]>>19))^(t^(t>>8));

        return XOR128_state.s[0];
}
float               PRNG::XOR128_produce_one_CPU(void)
{
	    return ((float) XOR128_produce_one_uint_CPU()) / 4294967295.0f;
}
void                PRNG::XOR128_produce_CPU(float* randoms_cpu)
{
    int number_of_prns_CPU = PRNG_samples * 4;
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = XOR128_produce_one_CPU();
}
void                PRNG::XOR128_produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = XOR128_produce_one_CPU();
}

void                PRNG::RANMAR_initialize(void)
{
        char buffer_prng_cl[FNAME_MAX_LENGTH];
        int  j = sprintf_s(buffer_prng_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
             j+= sprintf_s(buffer_prng_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_PRNG);
        random = GPU0->source_read(buffer_prng_cl);
                 GPU0->program_create(random);

        unsigned int half_seed_size = PRNG_instances;
        seeds_size              = half_seed_size * 2;
        seed_table_size         = seeds_size * (97+1);  // 97 - size of seed table, 1 - size of indises
        randoms_size            = PRNG_instances * PRNG_samples;

	    PRNG_seeds4             = (cl_uint4*)  calloc(seeds_size,     sizeof(cl_uint4));    // Input seeds
        PRNG_randoms            = (cl_float4*) calloc(randoms_size,   sizeof(cl_float4));   // Output randoms
        PRNG_seed_table_float4  = (cl_float4*) calloc(seed_table_size,sizeof(cl_float4));   // Seeds

        PRNG_seeds4[0].s[0]              = RM_seed1; // setup first thread as CPU
        PRNG_seeds4[0].s[1]              = RM_seed1; // setup first thread as CPU
        PRNG_seeds4[0].s[2]              = RM_seed1; // setup first thread as CPU
        PRNG_seeds4[0].s[3]              = RM_seed1; // setup first thread as CPU
        PRNG_seeds4[half_seed_size].s[0] = RM_seed2; // setup first thread as CPU
        PRNG_seeds4[half_seed_size].s[1] = RM_seed2; // setup first thread as CPU
        PRNG_seeds4[half_seed_size].s[2] = RM_seed2; // setup first thread as CPU
        PRNG_seeds4[half_seed_size].s[3] = RM_seed2; // setup first thread as CPU

        for (unsigned int i=1; i<seeds_size/2; i++) {
          PRNG_seeds4[i].s[0]                = (rand() % 31328);
          PRNG_seeds4[i].s[1]                = (rand() % 31328);
          PRNG_seeds4[i].s[2]                = (rand() % 31328);
          PRNG_seeds4[i].s[3]                = (rand() % 31328);
          PRNG_seeds4[i+half_seed_size].s[0] = (rand() % 30081);
          PRNG_seeds4[i+half_seed_size].s[1] = (rand() % 30081);
          PRNG_seeds4[i+half_seed_size].s[2] = (rand() % 30081);
          PRNG_seeds4[i+half_seed_size].s[3] = (rand() % 30081);
        }

        PRNG_seed_id       = GPU0->buffer_init(GPU0->buffer_type_Constant,  seeds_size,      PRNG_seeds4,            sizeof(cl_uint4));
    GPU0->print_stage("seeds initialized");
        PRNG_seed_table_id = GPU0->buffer_init(GPU0->buffer_type_IO,        seed_table_size, PRNG_seed_table_float4, sizeof(cl_float4));
    GPU0->print_stage("seed_table initialized");
        PRNG_randoms_id    = GPU0->buffer_init(GPU0->buffer_type_IO,        randoms_size,    PRNG_randoms,           sizeof(cl_float4));
    GPU0->print_stage("randoms initialized");

        int argument_id;

        const size_t global_size[]  = {PRNG_instances};

	    PRNG_randoms_seed_id = GPU0->kernel_init("rmseed",1,global_size,NULL);
            argument_id = GPU0->kernel_init_buffer(PRNG_randoms_seed_id,  PRNG_seed_id);
            argument_id = GPU0->kernel_init_buffer(PRNG_randoms_seed_id,  PRNG_seed_table_id);

	    int result = GPU0->kernel_run(PRNG_randoms_seed_id);                // Prepare seeds

        PRNG_randoms_kernel_id = GPU0->kernel_init("rmproduce",1,global_size,NULL); // local_size
	        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_seed_table_id);
		    argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_randoms_id);
		    argument_id = GPU0->kernel_init_constant(PRNG_randoms_kernel_id,&PRNG_samples);

        result |= GPU0->buffer_kill(PRNG_seed_id);
}
void                PRNG::RANMAR_initialize_CPU(void)
{

    RM_seed1 = rand() % 31328;
    RM_seed2 = rand() % 30081;

	int i = ((RM_seed1 / 177) % 177) + 2;
	int j = (RM_seed1 % 177) + 2;
	int k = ((RM_seed2 / 169) % 178) + 1;
	int l = (RM_seed2 % 169);
	for (int n = 0; n < 97; n++)
	{
		int s = 0;
		int t = 8388608;
		for (int m = 0; m < 24; m++)
		{
			int u = (i * j) % 179;
			u = (u * k) % 179;
			i = j;
			j = k;
			k = u;
			l = (53 * l + 1) % 169;
			if ((l * u) % 64 >= 32) s = s + t; 
			t = t >> 1;
		}
		RM_seeds[n] = ((float) s) / 16777216.0f;
	}
}
float               PRNG::RANMAR_produce_one_CPU(void)
{
	float uni = RM_seeds[RM_I97] - RM_seeds[RM_J97];
	if (uni < 0.0) {uni = uni + 1.0f;}
	RM_seeds[RM_I97] = uni;
	if (RM_I97==0){RM_I97=97;}
	if (RM_J97==0){RM_J97=97;}
	RM_I97--;
	RM_J97--;
	RM_C = (float) RM_C - (float) RM_CD;
	if (RM_C < 0.0) {RM_C = (float) RM_C + (float) RM_CM;}
	uni = uni - RM_C;
	if (uni < 0.0) {uni = uni + 1.0f;}

	return uni;
}
unsigned int        PRNG::RANMAR_check_seeds(void)
{
    // check seeds
    unsigned int result = 0;
    unsigned int* seedz = GPU0->buffer_map(PRNG_seed_table_id);
    int index = 0;
    for (int i=0; i<97; i++){
        if (GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+0])-RM_seeds[index]) result++;
        index++;
    }
    return result;
}
void                PRNG::RANMAR_print_seeds(void)
{
    // check seeds
    unsigned int* seedz = GPU0->buffer_map(PRNG_seed_table_id);
    int index = 0;
    for (int i=0; i<97; i++){
        double cpu_seed0 = RM_seeds[index];
        index++;
        double gpu_seed0 = GPU0->convert_to_float(seedz[(4*i*PRNG_instances)+0]);
        printf("\n[%2u]: (CPU)=%e (GPU)=%e",i,cpu_seed0,gpu_seed0); if (cpu_seed0!=gpu_seed0) printf(" !=");
    }
    printf("\n");
}
void                PRNG::RANMAR_produce_CPU(float* randoms_cpu)
{
    int number_of_prns_CPU = PRNG_samples;
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = RANMAR_produce_one_CPU();
}
void                PRNG::RANMAR_produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = RANMAR_produce_one_CPU();
}

void                PRNG::PM_initialize(void)
{
        char buffer_prng_cl[FNAME_MAX_LENGTH];
        int  j = sprintf_s(buffer_prng_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
             j+= sprintf_s(buffer_prng_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_PRNG);
        random = GPU0->source_read(buffer_prng_cl);
                 GPU0->program_create(random);

        seed_table_size         = PRNG_instances;
        randoms_size            = PRNG_instances * PRNG_samples;
        PRNG_seed_table_uint4   = (cl_uint4*)  calloc(seed_table_size,sizeof(cl_uint4));
        PRNG_randoms            = (cl_float4*) calloc(randoms_size,sizeof(cl_float4));

        PRNG_seed_table_uint4[0].s[0] = PMseed; // setup first thread as CPU
        PRNG_seed_table_uint4[0].s[1] = rand() % 2147483647;
        PRNG_seed_table_uint4[0].s[2] = rand() % 2147483647;
        PRNG_seed_table_uint4[0].s[3] = rand() % 2147483647;

        for (unsigned int i=1; i<seed_table_size; i++) {
          PRNG_seed_table_uint4[i].s[0] = rand() % 2147483647;
          PRNG_seed_table_uint4[i].s[1] = rand() % 2147483647;
          PRNG_seed_table_uint4[i].s[2] = rand() % 2147483647;
          PRNG_seed_table_uint4[i].s[3] = rand() % 2147483647;
        }

        PRNG_seed_table_id = GPU0->buffer_init(GPU0->buffer_type_IO, seed_table_size, PRNG_seed_table_uint4,    sizeof(cl_uint4));
        PRNG_randoms_id    = GPU0->buffer_init(GPU0->buffer_type_IO, randoms_size,    PRNG_randoms,             sizeof(cl_float4));

        int argument_id;
        const size_t global_size[]  = {PRNG_instances};

        PRNG_randoms_kernel_id   = GPU0->kernel_init("pm",1,global_size,NULL);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_seed_table_id);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_randoms_id);
        argument_id = GPU0->kernel_init_constant(PRNG_randoms_kernel_id,&PRNG_samples);
}
void                PRNG::PM_initialize_CPU(void)
{
        PMseed = rand() % 2147483647;
}
float               PRNG::PM_produce_one_CPU(void){
	unsigned int PM_lo, PM_hi;
    int PM_test;

	PM_hi = div(PMseed,PM_q).quot;
	PM_lo = PMseed % PM_q;

	PM_test = PM_a * PM_lo - PM_r * PM_hi;
	if (PM_test > 0) PMseed = PM_test;
	  else PMseed = PM_test + PM_m;

	float rnd = (float) (((double) PMseed) / PM_m_FP);

	return rnd;
}
void                PRNG::PM_produce_CPU(float* randoms_cpu)
{
    int number_of_prns_CPU = PRNG_samples;
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = PM_produce_one_CPU();
}
void                PRNG::PM_produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = PM_produce_one_CPU();
}

void                PRNG::XOR7_initialize(void)
{
        char buffer_prng_cl[FNAME_MAX_LENGTH];
        int  j = sprintf_s(buffer_prng_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
             j+= sprintf_s(buffer_prng_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_PRNG);
        random = GPU0->source_read(buffer_prng_cl);
                 GPU0->program_create(random);

        seed_table_size         = PRNG_instances * 2;
        randoms_size            = PRNG_instances * PRNG_samples;
        PRNG_seed_table_uint4   = (cl_uint4*)  calloc(seed_table_size,sizeof(cl_uint4));
        PRNG_randoms            = (cl_float4*) calloc(randoms_size,sizeof(cl_float4));

        for (unsigned int i=1; i<seed_table_size; i++)
        {
          PRNG_seed_table_uint4[i].s[0] = rand();
          PRNG_seed_table_uint4[i].s[1] = rand();
          PRNG_seed_table_uint4[i].s[2] = rand();
          PRNG_seed_table_uint4[i].s[3] = rand();
        }
printf("SSSSS %i\n", PRNG_seed_table_uint4[1].s[0]);
        PRNG_seed_table_uint4[0].s[0] = XOR7_state[0];    // setup first thread as CPU
        PRNG_seed_table_uint4[0].s[1] = XOR7_state[1];
        PRNG_seed_table_uint4[0].s[2] = XOR7_state[2];
        PRNG_seed_table_uint4[0].s[3] = XOR7_state[3];

        PRNG_seed_table_uint4[PRNG_instances].s[0] = XOR7_state[4];    // setup first thread as CPU
        PRNG_seed_table_uint4[PRNG_instances].s[1] = XOR7_state[5];
        PRNG_seed_table_uint4[PRNG_instances].s[2] = XOR7_state[6];
        PRNG_seed_table_uint4[PRNG_instances].s[3] = XOR7_state[7];

        PRNG_seed_table_id = GPU0->buffer_init(GPU0->buffer_type_IO, seed_table_size, PRNG_seed_table_uint4,    sizeof(cl_uint4));
        PRNG_randoms_id    = GPU0->buffer_init(GPU0->buffer_type_IO, randoms_size,    PRNG_randoms,             sizeof(cl_float4));

        int argument_id;
        const size_t global_size[]  = {PRNG_instances};

        PRNG_randoms_kernel_id   = GPU0->kernel_init("xor7",1,global_size,NULL);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_seed_table_id);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_randoms_id);
        argument_id = GPU0->kernel_init_constant(PRNG_randoms_kernel_id,&PRNG_samples);
}
void                PRNG::XOR7_initialize_CPU(void)
{
        XOR7_state[0] = rand();
        XOR7_state[1] = rand();
        XOR7_state[2] = rand();
        XOR7_state[3] = rand();

        XOR7_state[4]  = rand();
        XOR7_state[5]  = rand();
        XOR7_state[6]  = rand();
        XOR7_state[7]  = rand();
}
float               PRNG::XOR7_produce_one_CPU(void)
{
        unsigned int y, t;
        t = XOR7_state[(XOR7_index+7) & 7];     t = t ^ (t<<13);    y = t ^ (t<<9);
        t = XOR7_state[(XOR7_index+4) & 7];     y^= t ^ (t<<7);
        t = XOR7_state[(XOR7_index+3) & 7];     y^= t ^ (t>>3);
        t = XOR7_state[(XOR7_index+1) & 7];     y^= t ^ (t>>10);
        t = XOR7_state[XOR7_index];             t = t ^ (t>>7);     y^= t ^ (t<<24);
        XOR7_state[XOR7_index] = y;             XOR7_index = (XOR7_index+1) & 7;

        return ((float) ((float) y / XOR7_m_FP));
}
void                PRNG::XOR7_produce_CPU(float* randoms_cpu)
{
    int number_of_prns_CPU = PRNG_samples * 4;
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = XOR7_produce_one_CPU();
}
void                PRNG::XOR7_produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = XOR7_produce_one_CPU();
}

void                PRNG::RANECU_initialize(void)
{
        char buffer_prng_cl[FNAME_MAX_LENGTH];
        int  j = sprintf_s(buffer_prng_cl  ,FNAME_MAX_LENGTH,  "%s",GPU0->cl_root_path);
             j+= sprintf_s(buffer_prng_cl+j,FNAME_MAX_LENGTH-j,"%s",SOURCE_PRNG);
        random = GPU0->source_read(buffer_prng_cl);
                 GPU0->program_create(random);

        seed_table_size         = PRNG_instances * 2;
        randoms_size            = PRNG_instances * PRNG_samples;
        PRNG_seed_table_uint4   = (cl_uint4*)  calloc(seed_table_size,sizeof(cl_uint4));
        PRNG_randoms            = (cl_float4*) calloc(randoms_size,sizeof(cl_float4));

        for (unsigned int i=0; i<seed_table_size; i++)
        {
          PRNG_seed_table_uint4[i].s[0] = rand() % 2147483647;
          PRNG_seed_table_uint4[i].s[1] = rand() % 2147483647;
          PRNG_seed_table_uint4[i].s[2] = rand() % 2147483647;
          PRNG_seed_table_uint4[i].s[3] = rand() % 2147483647;
        }

        PRNG_seed_table_uint4[0].s[0] = RANECU_jseed1;    // setup first thread as CPU

        PRNG_seed_table_uint4[PRNG_instances].s[0] = RANECU_jseed2;    // setup first thread as CPU

        PRNG_seed_table_id = GPU0->buffer_init(GPU0->buffer_type_IO, seed_table_size, PRNG_seed_table_uint4,    sizeof(cl_uint4));
        PRNG_randoms_id    = GPU0->buffer_init(GPU0->buffer_type_IO, randoms_size,    PRNG_randoms,             sizeof(cl_float4));

        int argument_id;
        const size_t global_size[]  = {PRNG_instances};

        PRNG_randoms_kernel_id   = GPU0->kernel_init("ranecu",1,global_size,NULL);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_seed_table_id);
        argument_id = GPU0->kernel_init_buffer(PRNG_randoms_kernel_id,PRNG_randoms_id);
        argument_id = GPU0->kernel_init_constant(PRNG_randoms_kernel_id,&PRNG_samples);
}
void                PRNG::RANECU_initialize_CPU(void)
{
        RANECU_jseed1 = rand() % 2147483647;
        RANECU_jseed2 = rand() % 2147483647;
}
float               PRNG::RANECU_produce_one_CPU(void)
{
	int		RANECU_k;
	//
	RANECU_k = RANECU_jseed1 / RANECU_seedP11;
		RANECU_jseed1 = RANECU_seedP13 * (RANECU_jseed1 - RANECU_k * RANECU_seedP11) - RANECU_k * RANECU_seedP12;
		if (RANECU_jseed1 < 0) {RANECU_jseed1 = RANECU_jseed1 + RANECU_icons1;}

	RANECU_k = RANECU_jseed2 / RANECU_seedP21;
		RANECU_jseed2 = RANECU_seedP23 * (RANECU_jseed2 - RANECU_k * RANECU_seedP21) - RANECU_k * RANECU_seedP22;
		if (RANECU_jseed2 < 0) {RANECU_jseed2 = RANECU_jseed2 + RANECU_icons2;}

	int z = RANECU_jseed1 - RANECU_jseed2;
	if (z < 1) {z = z + RANECU_icons3;}

	float rnd =	(float) ((float) (z)) / ((float) RANECU_twom31);

	return rnd;
}
void                PRNG::RANECU_produce_CPU(float* randoms_cpu)
{
    int number_of_prns_CPU = PRNG_samples * 4;
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = RANECU_produce_one_CPU();
}
void                PRNG::RANECU_produce_CPU(float* randoms_cpu,int number_of_prns_CPU)
{
    for (int i=0; i<number_of_prns_CPU; i++) randoms_cpu[i] = RANECU_produce_one_CPU();
}

}
