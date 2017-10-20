/******************************************************************************
 * @file     o1cl.cl
 * @author   Vadim Demchik <vadimdi@yahoo.com>
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.0
 *
 * @brief    [QCDGPU]
 *           Defines general procedures for lattice update (O(1) gauge theory)
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

#ifndef O1CL_CL
#define O1CL_CL

                    HGPU_INLINE_PREFIX hgpu_float
o1_log_one_m_Ux(hgpu_float* Ux){
    hgpu_float log_one_m_Ux = log(1.0-(*Ux));
    return log_one_m_Ux;
}

                    HGPU_INLINE_PREFIX hgpu_float
o1_Phi(hgpu_float* log_one_m_Ux,hgpu_float* b,hgpu_float* eta){
    hgpu_float a = (-1.0+0.5*(*eta)*(*log_one_m_Ux))*(*log_one_m_Ux);
    hgpu_float result = 0.5*a*a*(*b);

    return result;
}

#ifdef ON_SUB_SCHEME
                    HGPU_INLINE_PREFIX_VOID void
o1_action(hgpu_float* S,hgpu_float* Ux,hgpu_float4* Uxmu,hgpu_float4* Uxmu_minus, __global const hgpu_float *  lattice_parameters){
#else
                    HGPU_INLINE_PREFIX_VOID void
o1_action(hgpu_float* S,hgpu_float* Ux,hgpu_float4* Uxmu, __global const hgpu_float *  lattice_parameters){
#endif
    hgpu_float z      = lattice_parameters[16];
    hgpu_float zeta   = lattice_parameters[18];
    hgpu_float eta    = lattice_parameters[19];
    hgpu_float b      = lattice_parameters[20];
    hgpu_float sqrt_z_lambda_zeta = lattice_parameters[21];  // 1/4*sqrt(z/(lambda*zeta))

    hgpu_float Ux1,Ux2,Ux3,Ux4;
    hgpu_float log_one_m_Ux_1;
    hgpu_float log_one_m_Ux_2;
    hgpu_float log_one_m_Ux_3;
    hgpu_float log_one_m_Ux_4;

    hgpu_float Phi1,Phi2,Phi3,Phi4;
    hgpu_float logPhi_1,logPhi_2,logPhi_3,logPhi_4;
    hgpu_float term = 0.0;

    hgpu_float log_one_m_Ux = o1_log_one_m_Ux(Ux);
    hgpu_float Phi  = o1_Phi(&log_one_m_Ux,&b,&eta);

    Ux1 = (*Uxmu).x;
    Ux2 = (*Uxmu).y;
    Ux3 = (*Uxmu).z;
    Ux4 = (*Uxmu).w;

    log_one_m_Ux_1 = o1_log_one_m_Ux(&Ux1);
    Phi1     = o1_Phi(&log_one_m_Ux_1,&b,&eta);
    logPhi_1 = log(Phi1/Phi);
    term     = logPhi_1*logPhi_1;

    log_one_m_Ux_2 = o1_log_one_m_Ux(&Ux2);
    Phi2 = o1_Phi(&log_one_m_Ux_2,&b,&eta);
    logPhi_2 = log(Phi2/Phi);
    term    += logPhi_2*logPhi_2;

    log_one_m_Ux_3 = o1_log_one_m_Ux(&Ux3);
    Phi3 = o1_Phi(&log_one_m_Ux_3,&b,&eta);
    logPhi_3 = log(Phi3/Phi);
    term    += logPhi_3*logPhi_3;

    log_one_m_Ux_4 = o1_log_one_m_Ux(&Ux4);
    Phi4 = o1_Phi(&log_one_m_Ux_4,&b,&eta);
    logPhi_4 = zeta*log(Phi4/Phi);
    term    += logPhi_4*logPhi_4;

#ifdef ON_SUB_SCHEME
    Ux1 = (*Uxmu_minus).x;
    Ux2 = (*Uxmu_minus).y;
    Ux3 = (*Uxmu_minus).z;
    Ux4 = (*Uxmu_minus).w;

    log_one_m_Ux_1 = o1_log_one_m_Ux(&Ux1);
    Phi1     = o1_Phi(&log_one_m_Ux_1,&b,&eta);
    logPhi_1 = log(Phi1/Phi);
    term    += logPhi_1*logPhi_1;

    log_one_m_Ux_2 = o1_log_one_m_Ux(&Ux2);
    Phi2 = o1_Phi(&log_one_m_Ux_2,&b,&eta);
    logPhi_2 = log(Phi2/Phi);
    term    += logPhi_2*logPhi_2;

    log_one_m_Ux_3 = o1_log_one_m_Ux(&Ux3);
    Phi3 = o1_Phi(&log_one_m_Ux_3,&b,&eta);
    logPhi_3 = log(Phi3/Phi);
    term    += logPhi_3*logPhi_3;

    log_one_m_Ux_4 = o1_log_one_m_Ux(&Ux4);
    Phi4 = o1_Phi(&log_one_m_Ux_4,&b,&eta);
    logPhi_4 = zeta*log(Phi4/Phi);
    term    += logPhi_4*logPhi_4;

    term = 0.5*term;
#endif

    (*S) =
            // potentail term
            -SKGROUPM1_2 * log(2.0-MIN(1.9999999999,sqrt(2.0*Phi)))
            -log((1.0-eta*log_one_m_Ux)/(1.0-(*Ux)))
            +Phi*(Phi-1.0

            // kinetic term
                +sqrt_z_lambda_zeta*term
            )/z;
}

                    HGPU_INLINE_PREFIX hgpu_float
o1_to_physical_field(gpu_o_1* Ux,hgpu_float* b,hgpu_float* eta){
        hgpu_float oU0 = (*Ux).uv1;

        hgpu_float log_one_m_Ux = o1_log_one_m_Ux(&oU0);
        hgpu_float Phi  = sqrt(2.0*o1_Phi(&log_one_m_Ux,b,eta));

    return Phi;
}

                    HGPU_INLINE_PREFIX hgpu_float
o1_correlator(gpu_o_1* Ux,gpu_o_1* Uy,hgpu_float* b,hgpu_float* eta){
    hgpu_float result = o1_to_physical_field(Ux,b,eta) * o1_to_physical_field(Uy,b,eta);
    return result;
}

#endif
