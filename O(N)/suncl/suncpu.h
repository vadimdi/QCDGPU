/******************************************************************************
 * @file     suncpu.h
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Procedures for host simulations (header), SU(3) gauge group and O(N) group
 *
 * @section  LICENSE
 *
 * Copyright (c) 2013, Vadim Demchik, Natalia Kolomoyets
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

#ifndef suncpu_h
#define suncpu_h

#include "../clinterface/clinterface.h"
#include "../random/random.h"
#include "../kernel/complex.h"

namespace SUN_CPU{
class SU {

    public:

            typedef struct o_1 {
                double u1;
            } o_1;

            typedef struct su_2 {
                hgpu_complex u1;
                hgpu_complex u2;
                hgpu_complex v1;
                hgpu_complex v2;
            } su_2;

            typedef struct su_3 {
                hgpu_complex u1;
                hgpu_complex u2;
                hgpu_complex u3;
                hgpu_complex v1;
                hgpu_complex v2;
                hgpu_complex v3;
                hgpu_complex w1;
                hgpu_complex w2;
                hgpu_complex w3;
            } su_3;

            typedef struct coords_4{
                int x;
                int y;
                int z;
                int t;
            } coords_4;

            static char**     directions;

           GPU_CL::GPU*     GPU0;                             // pointer to GPU instance
         PRNG_CL::PRNG*     PRNG0;                            // pointer to PRNG instance

        SU(void);
       ~SU(void);
       
           void             lattice_check_cpu(model_CL::model* lat);
           void             lattice_load_cpu(model_CL::model* lat,unsigned int* lattice_pointer);
           double*          lattice_avr_plaquette_cpu(model_CL::model* lat);
           double*          lattice_avr_plaquette_plq_cpu(model_CL::model* lat);
           double*          lattice_avr_Polyakov_loop_cpu(model_CL::model* lat);
           double           lattice_avr_Wilson_loop_cpu(model_CL::model* lat);
           double*          lattice_plaquette_cpu(model_CL::model* lat,coords_4 coords);
           su_3             lattice_get_staple_cpu(model_CL::model* lat,SU::coords_4 coords,int dir1);
           void             lattice_simulate(model_CL::model* lat,unsigned int* lattice_pointer);
           void             lattice_analysis_cpu(model_CL::model* lat);

           void             lattice_coords_print(coords_4 coords);
           unsigned int     lattice_coords_to_gid(model_CL::model* lat,coords_4 coords);
           unsigned int     lattice_coords_to_gid_half(model_CL::model* SUN0,SU::coords_4 coords);
           coords_4         lattice_neighbours_coords(model_CL::model* lat,const SU::coords_4 coord,int dir);
           coords_4         lattice_neighbours_coords_backward(model_CL::model* lat,const coords_4 coord,int dir);
           coords_4         lattice_neighbours_step2_gid(model_CL::model* lat,const coords_4 coord,const coords_4 stepz);
           coords_4         lattice_gid_to_coords(model_CL::model* lat,unsigned int gindex);
           su_3             lattice_table_3(model_CL::model* lat,coords_4 coords,unsigned int gindex,int dir);
           su_3             lattice_get_3(model_CL::model* lat,unsigned int * lattice_table,unsigned int gindex,int dir);
           su_3             lattice_get_raw_3(model_CL::model* lat,unsigned int * lattice_table,unsigned int gindex,int dir);
           void             lattice_store_3(model_CL::model* lat,su_3 matrix,unsigned int gindex,int dir);
           su_3             lattice_staple_3(model_CL::model* lat,unsigned int * lattice_staple,unsigned int gindex);
           su_3             lattice_staple_hermitian3(su_3 m1,su_3 m2,su_3 m3);
           su_3             lattice_staple_hermitian_backward3(su_3 m1,su_3 m2,su_3 m3);
           su_3             lattice_matrix_reconstruct3(su_3 a);
           su_3             lattice_matrix_hermitian(su_3 a);
           su_3             lattice_GramSchmidt_3(su_3 a);
           bool             lattice_matrix_compare3(su_3 a,su_3 b);
           void             lattice_matrix_print(su_3 matrix);
           void             lattice_matrix_print_double(su_3 matrix);
           su_2             lattice_matrix_times2(su_2 a,su_2 b);
           su_3             lattice_matrix_times3(su_3 a,su_3 b);
           su_3             lattice_matrix_add3(su_3 a,su_3 b);
           su_3             lattice_plaquette3(su_3 a,su_3 b,su_3 c,su_3 d);
           double           lattice_retrace(su_3 a);
           double           lattice_imtrace(su_3 a);
           su_2             lattice_unity2(void);
           su_3             lattice_unity3(void);
           su_3             lattice_lambda(int index);
           su_3             lattice_lambda1(void);
           su_3             lattice_lambda2(void);
           su_3             lattice_lambda3(void);
           su_3             lattice_lambda4(void);
           su_3             lattice_lambda5(void);
           su_3             lattice_lambda6(void);
           su_3             lattice_lambda7(void);
           su_3             lattice_lambda8(void);
           su_3             lattice_zero3(void);
           su_2             lattice_heatbath2(model_CL::model* lat,su_2 a,double beta,cl_float4* prns,unsigned int gindex);
           su_3             lattice_heatbath3(model_CL::model* lat,su_3 staple,su_3 U,double beta,cl_float4* prns,unsigned int gindex);
           unsigned int     lattice_odd_gid(model_CL::model* lat,unsigned int gid);

           o_1              lattice_table_o1(model_CL::model* lat,unsigned int gindex);
           void             lattice_store_o1(model_CL::model* lat,o_1 matrix,unsigned int gindex);
           o_1              lattice_get_o1(model_CL::model* lat,unsigned int * lattice_table,unsigned int gindex);
           double           o1_Phi(model_CL::model* lat,o_1 Ux);
           double           o1_to_physical_field(model_CL::model* lat,o_1 Ux);
           double           o1_action(model_CL::model* lat,o_1 Ux,o_1 U0,o_1 U1,o_1 U2,o_1 U3);
           double           o1_correlator(model_CL::model* lat,o_1 Ux,o_1 Uy);

           double*          lattice_action_o1_cpu(model_CL::model* lat);
           double*          lattice_avr_plaquette_plq_o1_cpu(model_CL::model* lat);
           double*          lattice_correlators_o1_cpu(model_CL::model* lat);

    private:
            su_3* lattice_data;                       // lattice data for CPU simulation
            o_1*  lattice_data_o1;                    // lattice data for CPU simulation (O(1) model)

};
}

#endif
