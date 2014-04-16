/******************************************************************************
 * @file     clinterface.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.4
 *
 * @brief    [QCDGPU]
 *           Interface for OpenCL AMD APP & nVidia SDK environment
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

#include "clinterface.h"

#define MAIN_CPP_FILE       "QCDGPU.cpp"

#ifdef _WIN32
#define PRG_SOURCE_DEBUG    "source.cl"
#define PRG_ERROR_DEBUG     "error.txt"
#define PATH_SEPARATOR      "\\"
#else
#define PRG_SOURCE_DEBUG    "source.cl"
#define PRG_ERROR_DEBUG     "error.txt"
#define PATH_SEPARATOR      "/"
#endif

#define FNAME_MAX_LENGTH    128

#define HASHES_SIZE         30

namespace GPU_CL{
using GPU_CL::GPU;

    char GPU::current_path[FILENAME_MAX]= "\0";     // FILENAME_MAX is definned by <stdio.h>

    GPU::GPU_debug_flags GPU::GPU_debug;

                GPU::GPU(void)
{
    // setup debug_level
    GPU_debug.wait_for_keypress = false;
    GPU_debug.profiling         = false;
    GPU_debug.brief_report      = false;
    GPU_debug.show_stage        = false;
    GPU_debug.local_run         = false;
    GPU_debug.rebuild_binary    = false;
    
    GPU_info.device_name        = NULL;
    GPU_info.local_memory_size  = 0;
    GPU_info.max_constant_size  = 0;
    GPU_info.max_memory_size    = 0;
    GPU_info.max_memory_width   = 0;
    GPU_info.max_workgroup_size = 0;
    GPU_info.memory_align_factor= 0;
    GPU_info.platform_vendor    = GPU::GPU_vendor_None;
    GPU_info.device_vendor      = GPU::GPU_vendor_None;

    GPU_current_kernel          = 0;    // current kernel counter
    GPU_current_buffer          = 0;    // current buffer counter
    GPU_current_program         = 0;    // current program counter
    GPU_active_program          = 0;    // active  program counter

    GPU_platform_id             = 0;    // GPU platform to be used (0 = first, 1 = second, etc)
    GPU_device_id               = 0;    // GPU device to be used   (0 = first, 1 = second, etc)

    GPU_inf_max_n               = 0;    // max number in .inf-filenames
 
    GPU_platforms_number        = 0;    // number of available platforms
    GPU_devices_number          = 0;    // number of available devices
    GPU_total_devices_number    = 0;    // total number of available devices

    GPU_max_work_item_dimensions= 0;    // CL_MAX_WORK_ITEM_DIMENSIONS
    GPU_max_work_item_sizes     = NULL; // CL_MAX_WORK_ITEMS_SIZES

    GPU_platform                = 0;    // utilized platform
    GPU_device                  = 0;    // utilized device
    GPU_context                 = NULL; // utilized context
    GPU_queue                   = NULL; // utilized command queue

    CPU_timers                  = 32;   // total number of reserved timers
    CPU_timer                   = NULL; // setup CPU timers
    CPU_current_timer_id        = 0;    // current timer id

    GPU_limit_max_workgroup_size= 0;    // manually limit max workgroup size

    CPU_timer = (int*) calloc((CPU_timers+1),sizeof(int));

    GPU_kernels = new kernels_hash[HASHES_SIZE];  // Hash for kernels
    GPU_buffers = new buffers_hash[HASHES_SIZE];  // Hash for buffers pointers 
    GPU_programs= new programs_hash[HASHES_SIZE]; // Hash for programs pointers

    cl_root_path = NULL;

    GPU_error = CL_SUCCESS;

    //__ MD5 section ________________________
    MD5_finalized               = false;
}
                GPU::~GPU(void)
{
    delete[] GPU_kernels;
    delete[] GPU_buffers;
    delete[] GPU_programs;

    free(cl_root_path);
    free(CPU_timer);
}

                GPU::kernels_hash::kernels_hash(void)
{
    kernel_name                 = NULL; // ponter to kernel's name declaration
    global_size                 = NULL; // global size
    local_size                  = NULL; // local size
    kernel_preferred_workgroup_size_multiple = 0; // preferred workgroup size
    argument_id                 = 0;
    work_dimensions             = 0;    // work dimensions
    kernel_local_mem_size       = 0;    // local memory size
    program_id                  = 0;    // program_id for GPU_programs array
    kernel_number_of_starts     = 0;    // total number of kernel starts - for deviation calculation
    kernel_elapsed_time         = 0.0;  // total kernel execution time (in nanoseconds)
    kernel_elapsed_time_squared = 0.0;  // total kernel execution time squared (in nanoseconds) - for deviation calculation
}
                GPU::kernels_hash::~kernels_hash(void)
{
      free((void*) (kernel_name));
      free(global_size);
      free(local_size);
}

                GPU::buffers_hash::buffers_hash(void)
{
      buffer_type                       = 0;
      host_ptr                          = NULL;
      size                              = 0;
      size_in_bytes                     = 0;        // buffer size in bytes ( size*sizeof(...) )
      mapped_ptr                        = NULL;     // ptr to corresponding host memory after mapping
                // profiling data __________________________
      buffer_write_elapsed_time         = 0.0;      // total buffer write time (in nanoseconds)
      buffer_write_elapsed_time_squared = 0.0;      // total buffer write time squared (in nanoseconds) - for deviation calculation
      buffer_write_number_of            = 0;        // total number of writes - for deviation calculation
      buffer_read_elapsed_time          = 0.0;      // total buffer read time (in nanoseconds)
      buffer_read_elapsed_time_squared  = 0.0;      // total buffer read time squared (in nanoseconds) - for deviation calculation
      buffer_read_number_of             = 0;        // total number of reads - for deviation calculation
}
                GPU::buffers_hash::~buffers_hash(void)
{
      free(host_ptr);
}

                GPU::programs_hash::programs_hash(void)
{
      source_ptr    = NULL;     // ponter to source's declaration
      options       = NULL;     // pointer to compiling options (-D... in clBuildProgram)
      md5           = NULL;     // md5 for source code
      datetime      = NULL;     // datetime of kernel compiling
      device        = NULL;     // device
      build_log     = NULL;     // program build log
      source_length = 0;
}
                GPU::programs_hash::~programs_hash(void)
{
      free((void*) (source_ptr));
      free((void*) (options));
      free((void*) (md5));
      free((void*) (datetime));
      free((void*) (device));
      free((void*) (build_log));
}

char*           GPU::trim(char* str)
{
    int start_str=0;
    int end_str=(int) strlen(str) - 1;
    if (end_str<=0) return str;
    bool start_flag = true;
    bool end_flag   = true;
    while ((start_flag) && (start_str<(int) strlen(str))) {
        if (str[start_str]==32) start_str++;
        else start_flag = false;
    }
    while ((end_flag) && (end_str>0)) {
        if (str[end_str]==32) end_str--;
        else end_flag = false;
    }
    int str_len = end_str - start_str + 1;
    if (str_len<0) str_len=0;
    char* result = (char*) calloc(str_len+1,sizeof(char));
    int j = 0;
    for (int i=start_str;i<end_str+1;i++){
        result[j] =str[i];
        j++;
    }
    result[j]=0;
    return result;
}

void            GPU::OpenCL_Check_Error(cl_int CL_Error_code,const char * CL_Error_description)
{
    using GPU_CL::GPU;
    if (CL_Error_code != CL_SUCCESS){
        printf("ERROR %i: (%s)\n", CL_Error_code, CL_Error_description);
        if (GPU_current_kernel>0)
        {
            char* clBuildLog;
            size_t clBuildLog_size, clBuildLog_actual_size;
            FILE *stream;
            char buffer[250];
            int j;

                clGetProgramInfo(GPU_programs[GPU_active_program].program,CL_PROGRAM_SOURCE,0,NULL,&clBuildLog_size);
                clBuildLog = (char*) calloc(clBuildLog_size+1, sizeof(char));
                clGetProgramInfo(GPU_programs[GPU_active_program].program,CL_PROGRAM_SOURCE,clBuildLog_size+1,clBuildLog,&clBuildLog_actual_size);

                j = sprintf_s(buffer  ,sizeof(buffer),  "%s",PRG_SOURCE_DEBUG);
                    fopen_s(&stream,buffer,"w+");
                    if(stream) {
                       fprintf(stream,clBuildLog);
                       if ( fclose(stream) ) {printf( "The file was not closed!\n" ); }
                    }

                clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_STATUS, 
                                        0, NULL, &clBuildLog_size );
                        clBuildLog = (char*) calloc(clBuildLog_size+1, sizeof(char));
                        clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_STATUS, 
                                        clBuildLog_size+1, clBuildLog, &clBuildLog_actual_size );
                        printf("\nBuild status:\n%s\n", clBuildLog);
                free(clBuildLog);

                    clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_OPTIONS, 
                                        0, NULL, &clBuildLog_size );
                        clBuildLog = (char*) calloc(clBuildLog_size+1, sizeof(char));
                        clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_OPTIONS, 
                                        clBuildLog_size+1, clBuildLog, &clBuildLog_actual_size );
                        printf("\nBuild options:\n%s\n", clBuildLog);
                free(clBuildLog);


                    clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_LOG, 
                                    0, NULL, &clBuildLog_size );
                        clBuildLog = (char*) calloc(clBuildLog_size+1, sizeof(char));
                        clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_LOG, 
                                    clBuildLog_size+1, clBuildLog, &clBuildLog_actual_size );
                        printf("\nLog:\n%s\n", clBuildLog);

                    j  = sprintf_s(buffer  ,sizeof(buffer),  "%s",PRG_ERROR_DEBUG);
                    fopen_s(&stream,buffer,"w+");
                    if(stream) {
                       fprintf(stream,clBuildLog);
                       if ( fclose(stream) ) {printf( "The file was not closed!\n" ); }
                    }
                free(clBuildLog);
        }
        device_finalize(-1);
        exit(-1);
        }
}
cl_uint         GPU::clGetDeviceInfoUint(cl_device_id device,cl_device_info inf)
{
        size_t result;
        OpenCL_Check_Error(clGetDeviceInfo(device,inf,sizeof(result), &result, NULL),(char*) "clGetDeviceInfo failed");
        return (cl_uint) result;
}
cl_ulong        GPU::clGetDeviceInfoUlong(cl_device_id device,cl_device_info device_info)
{
        cl_ulong result;
        OpenCL_Check_Error(clGetDeviceInfo(device,device_info,sizeof(result), &result, NULL),(char*) "clGetDeviceInfo failed");
        return result;
}
// ___ converting _________________________________________________________________________________
float           GPU::convert_to_float(unsigned int value)
{
        _Uint_and_Float u2fconv;
        u2fconv.uint_value[0]=value;
        return u2fconv.float_value;
}
float           GPU::convert_to_float(double value)
{
        return (float) value;
}
double          GPU::convert_to_double(float value)
{
        return (double) value;
}
double          GPU::convert_to_double(unsigned int value)
{
        _Uint_and_Float u2fconv;
        u2fconv.uint_value[0]=value;
        return (double) u2fconv.float_value;
}
double          GPU::convert_to_double(unsigned int value_LOW,unsigned int value_HIGH)
{
    Int_to_Double test;
        test.int_value[0] = value_LOW;
        test.int_value[1] = value_HIGH;
    return test.double_value;
}
unsigned int    GPU::convert_to_uint_HIGH(double x)
{
    Int_to_Double test;
        test.double_value= x;
    return test.int_value[1];
}
unsigned int    GPU::convert_to_uint_LOW(double x)
{
    Int_to_Double test;
        test.double_value = x;
    return test.int_value[0];
}
unsigned int    GPU::convert_to_uint(float value)
{
        Uint_and_Float u2fconv;
        u2fconv.float_value = value;
        return u2fconv.uint_value[0];
}
// ___ device _____________________________________________________________________________________
int             GPU::device_initialize(void)
{
    // initialize CPU timers by time of device initialization time
    int current_time = clock();
    for (int i = 0; i < CPU_timers; i++) CPU_timer[i] = current_time;

    // setup current path for loading .cl files
    if (!GetCurrentDir(current_path, FILENAME_MAX)) {
        GPU_error = GPU_error_device_initialization_failed;
        return 0;
    }

    // get .cl root path
    if (!cl_root_path){
        int j = 0;
        cl_root_path = (char*) calloc(FNAME_MAX_LENGTH,sizeof(char));
        char* j_prev = cl_root_path;
        char* j_cur  = cl_root_path;
        bool flag_found = false;

        while((!flag_found)&&(j_cur)){
            j = sprintf_s(cl_root_path,FNAME_MAX_LENGTH,"%s%s",current_path,PATH_SEPARATOR);
            j_cur = strstr(j_prev+1,PATH_SEPARATOR);
            if((j_cur) && (j_cur>j_prev)) {
                j_prev=j_cur;
                j = sprintf_s(j_cur+1,FNAME_MAX_LENGTH-(j_cur+1-cl_root_path),"%s",MAIN_CPP_FILE);

                if (is_file_exist(cl_root_path)) {
                    sprintf_s(j_cur,FNAME_MAX_LENGTH-(j_cur+1-cl_root_path),"%s",PATH_SEPARATOR);
                    flag_found = true;
                }
            }
        }
#ifdef _WIN32
        str_char_replace(cl_root_path,92,47);   // replace "\\" symbols with "/" in root cl path
#endif
    }


    // select desired device
    device_select(GPU_platform_id,GPU_device_id);

    // setup GPU_workgroupsize
    cl_uint GPU_max_work_item_dimensions = clGetDeviceInfoUint(GPU_device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);
    GPU_max_work_item_sizes = (size_t*) calloc((GPU_max_work_item_dimensions + 1),sizeof(size_t));
    OpenCL_Check_Error(clGetDeviceInfo(GPU_device,CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*GPU_max_work_item_dimensions, (void*) GPU_max_work_item_sizes, NULL),"clGetDeviceInfo failed");

    GPU_info.max_workgroup_size  = clGetDeviceInfoUint(GPU_device,CL_DEVICE_MAX_WORK_GROUP_SIZE);
    GPU_info.memory_align_factor = clGetDeviceInfoUint(GPU_device,CL_DEVICE_MAX_WORK_GROUP_SIZE);

    // create context
    GPU_context = clCreateContext(NULL,1,&GPU_device,NULL, NULL, &GPU_error);
    OpenCL_Check_Error(GPU_error,"clCreateContext failed");

    cl_command_queue_properties profiling_properties = 0;
    if (GPU_debug.profiling) profiling_properties = (CL_QUEUE_PROFILING_ENABLE);    // enable profiling for debuging
    GPU_queue = clCreateCommandQueue(GPU_context,GPU_device,profiling_properties,&GPU_error);
    OpenCL_Check_Error(GPU_error,"clCreateCommandQueue failed");


    GPU_info.local_memory_size  =          clGetDeviceInfoUlong(GPU_device,CL_DEVICE_LOCAL_MEM_SIZE);
    GPU_info.max_constant_size  =          clGetDeviceInfoUlong(GPU_device,CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE);
    GPU_info.max_memory_size    =          clGetDeviceInfoUlong(GPU_device,CL_DEVICE_MAX_MEM_ALLOC_SIZE);
    GPU_info.max_memory_width   = (size_t) clGetDeviceInfoUlong(GPU_device,CL_DEVICE_IMAGE3D_MAX_WIDTH);

    if (!GPU_info.max_memory_width)  GPU_info.max_memory_width  = 32;

    GPU_info.device_name = device_get_name(GPU_device);

    // GPU_debug.rebuild_binary
    if (GPU_debug.rebuild_binary) {
        // kill all .inf and .bin files
        int err = 0;
        int file_idx = 1;
        while (!err) {
            err = inf_file_delete(file_idx);
            file_idx++;
        }
    } else
    GPU_inf_max_n = inf_file_get_max_n();

    if (GPU_info.device_vendor == GPU_vendor_Intel) GPU_info.max_workgroup_size = 64;
    if ((GPU_limit_max_workgroup_size) && (GPU_limit_max_workgroup_size<GPU_info.max_workgroup_size)) GPU_info.max_workgroup_size = GPU_limit_max_workgroup_size;

    return GPU_devices_number;
}
int             GPU::device_finalize(int error_code)
{
    free(GPU_max_work_item_sizes);

    // clean GPU_kernels and programms

    for (int i=1; i<=GPU_current_program; i++){
            if (GPU_programs[i].program)
                clReleaseProgram(GPU_programs[i].program);
    }

    // clean GPU_buffers
    for (int i=1; i<GPU_current_buffer; i++) buffer_kill(i);

    // clean command queue and context
    if (GPU_queue) clReleaseCommandQueue(GPU_queue);
    if (GPU_context) clReleaseContext(GPU_context);

    if (GPU_debug.wait_for_keypress) {printf("\nPress any key to exit...\n"); _getch();}

    return error_code;
}
bool            GPU::device_auto_select(int platform_vendor,int vendor)
{
    bool supported_platform = false;    // is there supported platform?
    if (platform_vendor==GPU_vendor_any) supported_platform = true; // if vendor_any then do not filter platfroms
    bool supported_device   = false;    // is there supported device?
    cl_uint platforms_number;           // number of available platforms
    cl_uint devices_number;             // number of available devices
    cl_platform_id platform = NULL;     // supported platform
    cl_device_id device = NULL;         // supported device

    // platform initialization
    OpenCL_Check_Error(clGetPlatformIDs(0,NULL,&platforms_number),"clGetPlatformIDs failed");
    if(platforms_number==0){
       printf("There are no any available OpenCL platforms\n");
       exit(0);
    }

    // query available platforms
    char infobuf[4096];
    cl_platform_id* platforms = (cl_platform_id*) calloc(platforms_number,sizeof(cl_platform_id));
    OpenCL_Check_Error(clGetPlatformIDs(platforms_number,platforms,NULL),"clGetPlatformIDs failed");
    for(unsigned int i=0; i<platforms_number; i++)
    {
       platform = platforms[i];
       OpenCL_Check_Error(clGetPlatformInfo(platforms[i],CL_PLATFORM_VENDOR,sizeof(infobuf),infobuf,NULL),"clGetPlatformInfo failed");
       if (platform_vendor==GPU_vendor_Apple)
          supported_platform = (strstr(infobuf,"APPLE")==NULL) ? false : true;
       else if (platform_vendor==GPU_vendor_nVidia)
          supported_platform = (strstr(infobuf,"NVIDIA")==NULL) ? false : true;
       else if (platform_vendor==GPU_vendor_AMD)
          supported_platform = (strstr(infobuf,"Advanced Micro Devices")==NULL) ? false : true;
       else if (platform_vendor==GPU_vendor_Intel)
          supported_platform = (strstr(infobuf,"Intel")==NULL) ? false : true;

       if (supported_platform)
       {
          devices_number = 0;
          cl_int get_number_of_devices_result = clGetDeviceIDs(platform,CL_DEVICE_TYPE_ALL,0, 0, &devices_number);

          if (get_number_of_devices_result!=CL_DEVICE_NOT_FOUND) OpenCL_Check_Error(get_number_of_devices_result, "clGetDeviceIDs failed");
          if (devices_number>0)
          {
             // query available devices
             cl_device_id* devices = (cl_device_id*) calloc(devices_number,sizeof(cl_device_id));
             OpenCL_Check_Error(clGetDeviceIDs(platform,CL_DEVICE_TYPE_ALL,devices_number, devices, &devices_number),"clGetDeviceIDs failed");
             for(unsigned int t=0; t<devices_number; t++)
             {
                device = devices[t];
                OpenCL_Check_Error(clGetDeviceInfo(device,CL_DEVICE_VENDOR,sizeof(infobuf), &infobuf, NULL),"clGetDeviceInfo failed");
                if (vendor==GPU_vendor_Apple)
                   supported_device = (strstr(infobuf,"APPLE")==NULL) ? false : true;
                if (vendor==GPU_vendor_nVidia)
                   supported_device = (strstr(infobuf,"NVIDIA")==NULL) ? false : true;
                if (vendor==GPU_vendor_AMD)
                   supported_device = (strstr(infobuf,"Advanced Micro Devices")==NULL) ? false : true;
                if (vendor==GPU_vendor_Intel)
                   supported_device = (strstr(infobuf,"Intel")==NULL) ? false : true;
                if (vendor==GPU_vendor_any)
                   supported_device = true;

                if (supported_device)
                {
                   // setup supported platform and device
                   GPU_platform_id = i;
                   GPU_device_id   = t;
                   i = platforms_number;
                   t = devices_number;
                }
             }
             if (devices) delete[] devices;
          }
       }
    }
    return supported_device;
}

bool            GPU::device_select(unsigned int platform_id,unsigned int device_id)
{
        // platform initialization
        OpenCL_Check_Error(clGetPlatformIDs(0,NULL,&GPU_platforms_number),"clGetPlatformIDs failed");
        if(GPU_platforms_number==0){
           printf("There are no any available OpenCL platforms\n");
           exit(GPU_error_no_platform);
        }
        if (GPU_platforms_number<(platform_id+1)){
           printf("There is no OpenCL platform %u\n",platform_id);
           exit(GPU_error_no_platform);
        }
        char infobuf[4096];
        cl_platform_id* GPU_platforms = (cl_platform_id*) calloc(GPU_platforms_number,sizeof(cl_platform_id));
        OpenCL_Check_Error(clGetPlatformIDs(GPU_platforms_number,GPU_platforms,NULL),"clGetPlatformIDs failed");
        GPU_platform=GPU_platforms[platform_id];

    OpenCL_Check_Error(clGetPlatformInfo(GPU_platform,CL_PLATFORM_VENDOR,sizeof(infobuf),infobuf,NULL),"clGetPlatformInfo failed");
    if      (strstr(infobuf,"APPLE")!=NULL)                  { GPU_info.platform_vendor = GPU_vendor_Apple; }
    else if (strstr(infobuf,"NVIDIA")!=NULL)                 { GPU_info.platform_vendor = GPU_vendor_nVidia; }
    else if (strstr(infobuf,"Advanced Micro Devices")!=NULL) { GPU_info.platform_vendor = GPU_vendor_AMD; }
    else if (strstr(infobuf,"Intel")!=NULL)                  { GPU_info.platform_vendor = GPU_vendor_Intel; }

    //check total number of devices
    GPU_total_devices_number = 0;
    cl_uint devices_on_platform = 0; 
    for (unsigned int i=0; i<GPU_platforms_number; i++){
        devices_on_platform = 0;
        cl_int get_number_of_devices_result = clGetDeviceIDs(GPU_platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &devices_on_platform);
        if (get_number_of_devices_result!=CL_DEVICE_NOT_FOUND) OpenCL_Check_Error(get_number_of_devices_result, "clGetDeviceIDs failed");
        GPU_total_devices_number += devices_on_platform;
    }

    // device initialization
    cl_int get_number_of_devices_result = clGetDeviceIDs(GPU_platform,CL_DEVICE_TYPE_ALL, 0, NULL, &GPU_devices_number);
    if (get_number_of_devices_result!=CL_DEVICE_NOT_FOUND) OpenCL_Check_Error(get_number_of_devices_result, "clGetDeviceIDs failed");
    if(GPU_devices_number==0){
       printf("There are no any available OpenCL GPU devices\n");
       exit(GPU_error_no_device);
    }
    if (GPU_devices_number<(device_id+1)){
       printf("There is no OpenCL device %u\n",device_id);
       exit(GPU_error_no_device);
    }

        cl_device_id* GPU_devices = (cl_device_id*) calloc(GPU_devices_number,sizeof(cl_device_id));
        OpenCL_Check_Error(clGetDeviceIDs(GPU_platform,CL_DEVICE_TYPE_ALL,GPU_devices_number, GPU_devices, &GPU_devices_number),"clGetDeviceIDs failed");
        GPU_device = GPU_devices[device_id];

    OpenCL_Check_Error(clGetDeviceInfo(GPU_device,CL_DEVICE_VENDOR,sizeof(infobuf), &infobuf, NULL),"clGetDeviceInfo failed");
    if      (strstr(infobuf,"APPLE")!=NULL)                  { GPU_info.device_vendor = GPU_vendor_Apple; }
    else if (strstr(infobuf,"NVIDIA")!=NULL)                 { GPU_info.device_vendor = GPU_vendor_nVidia; }
    else if (strstr(infobuf,"Advanced Micro Devices")!=NULL) { GPU_info.device_vendor = GPU_vendor_AMD; }
    else if (strstr(infobuf,"Intel")!=NULL)                  { GPU_info.device_vendor = GPU_vendor_Intel; }

    if(GPU_platforms) delete[] GPU_platforms;
    if(GPU_devices)   delete[] GPU_devices;

    GPU_platform_id = platform_id;
    GPU_device_id   = device_id;

    return true;
}

char*           GPU::device_get_name(cl_device_id device){
    size_t result_length = 4096;
    size_t result_actual_size;

        char* result = (char*) calloc(result_length,sizeof(char));
        OpenCL_Check_Error(clGetDeviceInfo(device,CL_DEVICE_NAME,result_length, (void*) result, &result_actual_size),"clGetDeviceInfo failed");
        result = trim(result);
        result = (char*) realloc(result,result_actual_size * sizeof(char));

    return result;
}

char*           GPU::platform_get_name(cl_platform_id platform){
    size_t result_length = 4096;
    size_t result_actual_size;

        char* result = (char*) calloc(result_length,sizeof(char));
        OpenCL_Check_Error(clGetPlatformInfo(platform,CL_PLATFORM_NAME,result_length, (void*) result, &result_actual_size),"clGetPlatformInfo failed");
        result = trim(result);
        result = (char*) realloc(result,result_actual_size * sizeof(char));

    return result;
}

// ___ source _____________________________________________________________________________________
char*           GPU::source_read(const char* file_name)
{
        // create a CL program using the kernel source
        FILE * cl_kernels_file;
        char buffer[1024];

        sprintf_s(buffer,sizeof(buffer),"%s",file_name);
        fopen_s(&cl_kernels_file,buffer,"r");
           if(!cl_kernels_file)
           {
              printf("File: %s\n",buffer);
              OpenCL_Check_Error(1,".cl-file not found");
              return NULL;
           }
           fseek (cl_kernels_file , 0 , SEEK_END);
           unsigned int flen = ftell(cl_kernels_file);
           fseek (cl_kernels_file , 0 , SEEK_SET);
           char* cl_kernels_source = (char*) calloc(flen + 2, sizeof(char));
           unsigned int freadlen = (unsigned int) fread( cl_kernels_source, sizeof(char), flen, cl_kernels_file );
           if (freadlen==0) return NULL;
        if ( fclose(cl_kernels_file) ) printf( "The file was not closed!\n" );
    return cl_kernels_source;
}

char*           GPU::source_add(char* source, const char* file_name)
{
    int error_code = 0;
    char* temporary_source = source_read(file_name);
    if (temporary_source==NULL) return source;

    const char* source_separator = "\n\n\n";  // separator between two source codes
    int temporary_source_length = (int) strlen(temporary_source) + 1;
    int source_length = (int) strlen(source)+1;
    int result_source_length = source_length + temporary_source_length + 10;
    char* result_source  = (char*) calloc(result_source_length, sizeof(char));
    char* result_source2 = (char*) calloc(result_source_length, sizeof(char));
    error_code  = strncpy_s(result_source, result_source_length * sizeof(char), source, source_length);
    error_code += (int) (strlen(result_source) + strlen(source_separator) - _snprintf_s(result_source2, result_source_length * sizeof(char), _TRUNCATE, "%s%s", result_source, source_separator));
    error_code += (int) (strlen(result_source2) + strlen(temporary_source) - _snprintf_s(result_source, result_source_length * sizeof(char), _TRUNCATE, "%s%s", result_source2, temporary_source));
    if (error_code) return source;

    free(temporary_source);
    free(result_source2);

    return result_source;
}

// ___ program ____________________________________________________________________________________
int             GPU::program_create(const char* source){
    return program_create(source,NULL);
}

int             GPU::program_create(const char* source,const char* options){
    start_timer_CPU(10);

    FILE * cl_program_file = NULL;
    char buffer[FNAME_MAX_LENGTH];
    char buffer_inf[FNAME_MAX_LENGTH];

    GPU_current_program++;
    GPU_active_program = GPU_current_program;

    int GPU_active_file = 0;

    // setup reserve kernel's source
    int source_length = (int) strlen(source) + 1;
    char* temporary_source = (char*) calloc(source_length + 1, sizeof(char));
    
    strncpy_s(temporary_source, source_length, source, source_length);
    GPU_programs[GPU_active_program].source_ptr     = temporary_source;
    GPU_programs[GPU_active_program].source_length  = source_length;
    if (options!=NULL){
       int options_length = (int) strlen(options) + 1;
       char* temporary_options = (char*) calloc(options_length + 1, sizeof(char));
       strncpy_s(temporary_options, options_length, options, options_length);
        GPU_programs[GPU_active_program].options        = temporary_options;
    } else {
       GPU_programs[GPU_active_program].options        = NULL;
    }

    GPU_programs[GPU_active_program].md5      = MD5(source);
    GPU_programs[GPU_active_program].device   = device_get_name(GPU_device);
    GPU_programs[GPU_active_program].platform = platform_get_name(GPU_platform);
    GPU_programs[GPU_active_program].datetime = get_current_datetime();

    for (int i = 1; i <= GPU_inf_max_n; i++){
       sprintf_s(buffer_inf,FNAME_MAX_LENGTH,"program%u.inf",i);
       // get .inf-file
       int parameters_items = 0;
       GPU_init_parameters* parameters = get_init_file(buffer_inf);
       if (parameters==NULL) {delete[] &buffer_inf; return 1;}
       bool inf_number  = false;
       bool inf_md5     = false;
       bool inf_device  = false;
       bool inf_platform= false;
       bool inf_options = false;
       int  inf_idx_md5     = 0;
       int  inf_idx_device  = 0;
       int  inf_idx_platform= 0;
       int  inf_idx_options = 0;
       bool parameters_flag = false;
       while(!parameters_flag){
          if (!strcmp(parameters[parameters_items].Variable,"NUMBER"))  {
             if(parameters[parameters_items].iVarVal == GPU_active_program) inf_number = true;
           }
          if (!strcmp(parameters[parameters_items].Variable,"MD5"))  {
             inf_md5 = true;
             inf_idx_md5 = parameters_items;
           }
          if (!strcmp(parameters[parameters_items].Variable,"DEVICE"))  {
             inf_device = true;
             inf_idx_device = parameters_items;
           }
          if (!strcmp(parameters[parameters_items].Variable,"PLATFORM"))  {
             inf_platform = true;
             inf_idx_platform = parameters_items;
           }
          if (!strcmp(parameters[parameters_items].Variable,"OPTIONS"))  {
             inf_options = true;
             inf_idx_options = parameters_items;
           }
          parameters_flag = parameters[parameters_items].final;
          parameters_items++;
        }
       // check parameters        
       if (inf_number){
            // check MD5
            if (!strcmp(parameters[inf_idx_md5].txtVarVal,GPU_programs[GPU_active_program].md5)){
                // check device and options
                if ((!strcmp(parameters[inf_idx_device].txtVarVal,  GPU_programs[GPU_active_program].device))&&
                    (!strcmp(parameters[inf_idx_platform].txtVarVal,GPU_programs[GPU_active_program].platform))&&
                    (((GPU_programs[GPU_active_program].options==NULL)&&(!inf_options))||
                     ((GPU_programs[GPU_active_program].options!=NULL)&&(!strcmp(parameters[inf_idx_options].txtVarVal,GPU_programs[GPU_active_program].options)))))
                    {
                        GPU_active_file = i;
                        i = GPU_inf_max_n;
                    }
            // skip current .inf-file
            } else {
                // MD5 does not match
                // kill current .inf-file and rename last .inf-file into current
                inf_file_delete(i);
                if (i < GPU_inf_max_n) inf_file_rename(GPU_inf_max_n,i);    // check if file is last
                GPU_inf_max_n--;
                i--;
            }            
        }
     }

    // if GPU_active_file == 0 then create files
    if (GPU_active_file == 0) GPU_active_file = GPU_inf_max_n + 1;

    int j = sprintf_s(buffer,    FNAME_MAX_LENGTH,"program%u.bin",GPU_active_file);
        j = sprintf_s(buffer_inf,FNAME_MAX_LENGTH,"program%u.inf",GPU_active_file);
    if (GPU_active_file <= GPU_inf_max_n)  {
        // try to open .bin-file
        fopen_s(&cl_program_file,buffer,"rb");
     }

    if((GPU_active_file > GPU_inf_max_n)||(!cl_program_file)){
        printf("\nprogram%u.bin is being compiled... \n",GPU_active_program);
        GPU_programs[GPU_active_program].program = clCreateProgramWithSource(GPU_context, 1,&GPU_programs[GPU_active_program].source_ptr, NULL, &GPU_error);
        OpenCL_Check_Error(GPU_error,"clCreateProgramWithSource failed");
        OpenCL_Check_Error(clBuildProgram(GPU_programs[GPU_active_program].program, 1, &GPU_device, GPU_programs[GPU_active_program].options, NULL, NULL),"clBuildProgram failed");

        cl_uint num_devices;
        OpenCL_Check_Error(clGetProgramInfo(GPU_programs[GPU_active_program].program, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint), &num_devices, NULL),"clGetProgramInfo1 failed");

        cl_device_id* devices = (cl_device_id*) calloc(num_devices, sizeof(cl_device_id));
        OpenCL_Check_Error(clGetProgramInfo(GPU_programs[GPU_active_program].program, CL_PROGRAM_DEVICES, num_devices * sizeof(cl_device_id), devices, 0),"clGetProgramInfo2 failed");

        size_t* binary_sizes = (size_t*) calloc(num_devices, sizeof(size_t));    
        OpenCL_Check_Error(clGetProgramInfo(GPU_programs[GPU_active_program].program, CL_PROGRAM_BINARY_SIZES, num_devices * sizeof(size_t), binary_sizes, NULL),"clGetProgramInfo3 failed");

        char** binary = (char**) calloc(num_devices, sizeof(char*));
        for( unsigned int i=0; i<num_devices; ++i) binary[i]= (char*) malloc(binary_sizes[i]);

        OpenCL_Check_Error(clGetProgramInfo(GPU_programs[GPU_active_program].program, CL_PROGRAM_BINARIES, num_devices * sizeof(size_t), binary, NULL),"clGetProgramInfo4 failed");

                size_t size;
                OpenCL_Check_Error(clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size),"clGetProgramBuildInfo failed");
                if (size>4) {
                   GPU_programs[GPU_active_program].build_log = (char*) calloc(size,sizeof(char));
                   OpenCL_Check_Error(clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_LOG, size, GPU_programs[GPU_active_program].build_log, NULL),"clGetProgramBuildInfo failed");
                   if (GPU_debug.brief_report) printf("Program buid log: [%s]\n", GPU_programs[GPU_active_program].build_log);
                 }

        unsigned int idx = 0;
        while( idx<num_devices && devices[idx] != GPU_device ) ++idx;

        // save inf file
        fopen_s(&cl_program_file,buffer_inf,"wb");
        if ( (cl_program_file) && (idx < num_devices) && (binary_sizes[idx]>0) ){
            fprintf(cl_program_file,"NUMBER=%u\n",GPU_active_program);
            fprintf(cl_program_file,"MD5=%s\n",GPU_programs[GPU_active_program].md5);
            if (GPU_programs[GPU_active_program].options!=NULL)
                fprintf(cl_program_file,"OPTIONS=%s\n",GPU_programs[GPU_active_program].options);
            fprintf(cl_program_file,"DEVICE=%s\n",GPU_programs[GPU_active_program].device);
            fprintf(cl_program_file,"PLATFORM=%s\n",GPU_programs[GPU_active_program].platform);
            fprintf(cl_program_file,"DATE=%s",GPU_programs[GPU_active_program].datetime);
            if ( fclose(cl_program_file) ) printf( "The file was not closed!\n" );
        }

        printf("program%u.bin compilation done (%f seconds)!\n",GPU_active_program,get_timer_CPU(10));

        // save binary file
        fopen_s(&cl_program_file,buffer,"wb");
        if ( (cl_program_file) && (idx < num_devices) && (binary_sizes[idx]>0) ){
           fwrite(binary[idx],1,binary_sizes[idx],cl_program_file);
           if ( fclose(cl_program_file) ) printf( "The file was not closed!\n" );
        }
        free(devices);
        free(binary_sizes);
        for( unsigned int i=0; i<num_devices; ++i) free(binary[i]);
        free(binary);
        GPU_inf_max_n++;
    } else {

        // load binary file
        fseek (cl_program_file, 0, SEEK_END);
        const size_t binary_size = ftell(cl_program_file);
        rewind(cl_program_file);
        unsigned char* binary;
        binary = (unsigned char*) malloc (binary_size);
        fread(binary, 1, binary_size, cl_program_file);
        fclose(cl_program_file);

        cl_int status;
        GPU_programs[GPU_active_program].program = clCreateProgramWithBinary(GPU_context, 1, &GPU_device, &binary_size, (const unsigned char**)&binary, &status, &GPU_error);
        OpenCL_Check_Error(status,"clCreateProgramWithBinary failed");
        OpenCL_Check_Error(GPU_error,"clCreateProgramWithBinary failed");

        OpenCL_Check_Error(clBuildProgram(GPU_programs[GPU_active_program].program, 1, &GPU_device, GPU_programs[GPU_active_program].options, NULL, NULL),"clBuildProgram failed");

                size_t size;
                OpenCL_Check_Error(clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size),"clGetProgramBuildInfo failed");
                if (size>4) {
                   GPU_programs[GPU_active_program].build_log = (char*) calloc(size,sizeof(char));
                   OpenCL_Check_Error(clGetProgramBuildInfo(GPU_programs[GPU_active_program].program, GPU_device, CL_PROGRAM_BUILD_LOG, size, GPU_programs[GPU_active_program].build_log, NULL),"clGetProgramBuildInfo failed");
                   if (GPU_debug.brief_report) printf("Program buid log: [%s]\n", GPU_programs[GPU_active_program].build_log);
                 }
        free(binary);
    }

    // setup program

    return GPU_active_program;
}

int             GPU::program_set_active(int program_id){
    int result = -1;    // default: error
    if ((program_id>0)&&(program_id<=GPU_current_program)) {
        result = 0;
        GPU_active_program = program_id;
    }
    return result;
}

int             GPU::program_get_active(void){
    return GPU_active_program;
}

// ___ kernel _____________________________________________________________________________________
int             GPU::kernel_init(const char* kernel_name, unsigned int work_dimensions, const size_t* global_size, const size_t* local_size)
{
    GPU_current_kernel++;

    // setup program
    GPU_kernels[GPU_current_kernel].program_id  = GPU_active_program;
    GPU_kernels[GPU_current_kernel].program     = GPU_programs[GPU_active_program].program;

    // setup kernel
    GPU_kernels[GPU_current_kernel].kernel = clCreateKernel( GPU_kernels[GPU_current_kernel].program, kernel_name, &GPU_error );
    OpenCL_Check_Error(GPU_error,"clCreateKernel failed");

    // setup reserve kernel's name
    unsigned int temporary_kernel_name_length = (unsigned int) strlen(kernel_name)+1;
    char* temporary_kernel_name = (char*) calloc(temporary_kernel_name_length, sizeof(char));
    strncpy_s(temporary_kernel_name, temporary_kernel_name_length * sizeof(char), kernel_name, temporary_kernel_name_length);
    GPU_kernels[GPU_current_kernel].kernel_name = temporary_kernel_name;

    // setup reserve kernel's argument counter
    GPU_kernels[GPU_current_kernel].argument_id = 0;

    OpenCL_Check_Error(clGetKernelWorkGroupInfo(GPU_kernels[GPU_current_kernel].kernel,GPU_device,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,sizeof(GPU_kernels[GPU_current_kernel].kernel_preferred_workgroup_size_multiple),&GPU_kernels[GPU_current_kernel].kernel_preferred_workgroup_size_multiple,NULL),"clGetKernelWorkGroupInfo failed");
    OpenCL_Check_Error(clGetKernelWorkGroupInfo(GPU_kernels[GPU_current_kernel].kernel,GPU_device,CL_KERNEL_LOCAL_MEM_SIZE,sizeof(GPU_kernels[GPU_current_kernel].kernel_local_mem_size),&GPU_kernels[GPU_current_kernel].kernel_local_mem_size,NULL),"clGetKernelWorkGroupInfo failed");
    size_t kernel_work_group_size = 0;
    OpenCL_Check_Error(clGetKernelWorkGroupInfo(GPU_kernels[GPU_current_kernel].kernel,GPU_device,CL_KERNEL_WORK_GROUP_SIZE,sizeof(size_t),&kernel_work_group_size,NULL),"clGetKernelWorkGroupInfo failed");
    kernel_work_group_size = (unsigned int) (1<<((int) floor(log((double) kernel_work_group_size)/log(2.0))));
    
    if (GPU_info.device_vendor == GPU::GPU_vendor_Intel) kernel_work_group_size =  64;
    
    if((GPU_limit_max_workgroup_size)&&(GPU_limit_max_workgroup_size < kernel_work_group_size)) kernel_work_group_size = GPU_limit_max_workgroup_size;

    // setup kernel's size
    size_t* temporary_local_size  = (size_t*) calloc(work_dimensions+1,sizeof(size_t));

    if ((local_size) && (local_size[0]!=0)) temporary_local_size[0] = local_size[0];
    else temporary_local_size[0] = kernel_work_group_size;

    size_t* temporary_global_size = (size_t*) calloc(work_dimensions+1,sizeof(size_t));
    for (unsigned int i=0; i<work_dimensions; i++) temporary_global_size[i] = global_size[i];

    GPU_kernels[GPU_current_kernel].global_size                 = temporary_global_size;
    GPU_kernels[GPU_current_kernel].local_size                  = temporary_local_size;
    GPU_kernels[GPU_current_kernel].work_dimensions             = work_dimensions;

    // setup profiling data __________________________________________
    GPU_kernels[GPU_current_kernel].kernel_start                = 0;
    GPU_kernels[GPU_current_kernel].kernel_finish               = 0;
    GPU_kernels[GPU_current_kernel].kernel_elapsed_time         = 0.0;
    GPU_kernels[GPU_current_kernel].kernel_elapsed_time_squared = 0.0;
    GPU_kernels[GPU_current_kernel].kernel_number_of_starts     = 0;

    return GPU_current_kernel;
}

int             GPU::kernel_init_buffer(int kernel_id,int buffer_id)
{
    if (GPU_buffers[buffer_id].buffer_type==buffer_type_LDS)
     {
       OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, GPU_kernels[kernel_id].argument_id, GPU_buffers[buffer_id].size_in_bytes, NULL),"clSetKernelArg failed");
     } else {
       OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, GPU_kernels[kernel_id].argument_id, sizeof(GPU_buffers[buffer_id].buffer), (void*) &GPU_buffers[buffer_id].buffer),"clSetKernelArg failed");
     }
    GPU_kernels[kernel_id].argument_id++;

    return GPU_kernels[kernel_id].argument_id;
}

int             GPU::kernel_init_constant(int kernel_id,int*    host_ptr)
{
        OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, GPU_kernels[kernel_id].argument_id, sizeof(int), (void*) host_ptr),"clSetKernelArg failed");
        GPU_kernels[kernel_id].argument_id++;
        return GPU_kernels[kernel_id].argument_id;
}

int             GPU::kernel_init_constant_reset(int kernel_id,int* host_ptr, int argument_id)
{
        OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, argument_id, sizeof(int), (void*) host_ptr),"clSetKernelArg failed");
        return argument_id;
}

int             GPU::kernel_init_constant(int kernel_id,cl_uint4* host_ptr)
{
        OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, GPU_kernels[kernel_id].argument_id, sizeof(cl_uint4), (void*) host_ptr),"clSetKernelArg failed");
        GPU_kernels[kernel_id].argument_id++;
        return GPU_kernels[kernel_id].argument_id;
}

int             GPU::kernel_init_constant(int kernel_id,float*  host_ptr)
{
        OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, GPU_kernels[kernel_id].argument_id, sizeof(float), (void*) host_ptr),"clSetKernelArg failed");
        GPU_kernels[kernel_id].argument_id++;
        return GPU_kernels[kernel_id].argument_id;
}

int             GPU::kernel_init_constant(int kernel_id,double* host_ptr)
{
        OpenCL_Check_Error(clSetKernelArg(GPU_kernels[kernel_id].kernel, GPU_kernels[kernel_id].argument_id, sizeof(double), (void*) host_ptr),"clSetKernelArg failed");
        GPU_kernels[kernel_id].argument_id++;
        return GPU_kernels[kernel_id].argument_id;
}
int             GPU::kernel_run(int kernel_id)
{
    cl_event kernel_event;
    if (GPU_debug.profiling){
        // run with profiling
        cl_ulong kernel_start, kernel_finish;

        OpenCL_Check_Error(clEnqueueNDRangeKernel(GPU_queue,GPU_kernels[kernel_id].kernel,GPU_kernels[kernel_id].work_dimensions,NULL,GPU_kernels[kernel_id].global_size,GPU_kernels[kernel_id].local_size, 0, NULL, &kernel_event),"clEnqueueNDRangeKernel failed");
        OpenCL_Check_Error(clWaitForEvents(1, &kernel_event),"clWaitForEvents failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &kernel_finish, 0),"clGetEventProfilingInfo failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(kernel_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &kernel_start,  0),"clGetEventProfilingInfo failed");
        double elapsed_time = (double) (kernel_finish-kernel_start);
        GPU_kernels[kernel_id].kernel_elapsed_time          += elapsed_time;
        GPU_kernels[kernel_id].kernel_elapsed_time_squared  += elapsed_time * elapsed_time;
        GPU_kernels[kernel_id].kernel_start                  = kernel_start;
        GPU_kernels[kernel_id].kernel_finish                 = kernel_finish;
        GPU_kernels[kernel_id].kernel_number_of_starts++;
    } else {
        // run without profiling
        size_t local_workgroup_size;
        OpenCL_Check_Error(clGetKernelWorkGroupInfo(GPU_kernels[kernel_id].kernel,GPU_device,CL_KERNEL_WORK_GROUP_SIZE,sizeof(size_t),&local_workgroup_size,NULL),"clGetKernelWorkGroupInfo failed");
        
        OpenCL_Check_Error(clEnqueueNDRangeKernel(GPU_queue,GPU_kernels[kernel_id].kernel,GPU_kernels[kernel_id].work_dimensions,NULL,GPU_kernels[kernel_id].global_size,GPU_kernels[kernel_id].local_size, 0, NULL, &kernel_event),"clEnqueueNDRangeKernel failed");
        OpenCL_Check_Error(clWaitForEvents(1, &kernel_event),"clWaitForEvents failed");
    }
    OpenCL_Check_Error(clFinish(GPU_queue),"clFinish failed");
    return kernel_id;
}

int             GPU::kernel_get_worksize(int kernel_id){
        size_t result;

    if (GPU_kernels[kernel_id].local_size==NULL)
        OpenCL_Check_Error(clGetKernelWorkGroupInfo(GPU_kernels[kernel_id].kernel,GPU_device,CL_KERNEL_WORK_GROUP_SIZE,sizeof(result),&result,NULL),"clGetKernelWorkGroupInfo failed");
    else
        result = GPU_kernels[kernel_id].local_size[0];
    return (int) result;
}

GPU::GPU_time_deviation  GPU::kernel_get_execution_time(int kernel_id){
    GPU_time_deviation execution_time = time_get_deviation(GPU_kernels[kernel_id].kernel_elapsed_time,GPU_kernels[kernel_id].kernel_elapsed_time_squared,GPU_kernels[kernel_id].kernel_number_of_starts);
    return execution_time;
}

// ___ buffer _____________________________________________________________________________________
int             GPU::buffer_init(int buffer_type, int size, void* host_ptr, int size_of)
{
    GPU_current_buffer++;

    // setup profiling data __________________________________________________
    GPU_buffers[GPU_current_buffer].buffer_write_start                  = 0;
    GPU_buffers[GPU_current_buffer].buffer_write_finish                 = 0;
    GPU_buffers[GPU_current_buffer].buffer_write_elapsed_time           = 0.0;
    GPU_buffers[GPU_current_buffer].buffer_write_elapsed_time_squared   = 0.0;
    GPU_buffers[GPU_current_buffer].buffer_write_number_of              = 0;

    GPU_buffers[GPU_current_buffer].buffer_read_start                   = 0;
    GPU_buffers[GPU_current_buffer].buffer_read_finish                  = 0;
    GPU_buffers[GPU_current_buffer].buffer_read_elapsed_time            = 0.0;
    GPU_buffers[GPU_current_buffer].buffer_read_elapsed_time_squared    = 0.0;
    GPU_buffers[GPU_current_buffer].buffer_read_number_of               = 0;

        cl_mem_flags flags = 0;
        cl_mem_flags flags2 = 0;
        switch (buffer_type) {
                case buffer_type_Input:      { flags = CL_MEM_READ_ONLY;  flags2 = CL_MEM_COPY_HOST_PTR; break;}
                case buffer_type_Constant:   { flags = CL_MEM_READ_ONLY;  flags2 = CL_MEM_COPY_HOST_PTR; break;}
                case buffer_type_IO:         { flags = CL_MEM_READ_WRITE; flags2 = CL_MEM_COPY_HOST_PTR; break;}
                case buffer_type_Output:     { flags = CL_MEM_WRITE_ONLY;                                break;}
                case buffer_type_Global:     { flags = CL_MEM_READ_WRITE;                                break;}
        }
        GPU_buffers[GPU_current_buffer].buffer_type = buffer_type;
        GPU_buffers[GPU_current_buffer].size = size;
        GPU_buffers[GPU_current_buffer].size_in_bytes = size * size_of;
        GPU_buffers[GPU_current_buffer].host_ptr = host_ptr;
        GPU_buffers[GPU_current_buffer].mapped_ptr = NULL;

    if (buffer_type!=buffer_type_LDS){
       flags = flags | flags2;

       GPU_buffers[GPU_current_buffer].buffer = clCreateBuffer( GPU_context,flags,size * size_of,host_ptr, &GPU_error );
       OpenCL_Check_Error(GPU_error,"clCreateKernel failed");
       GPU_error = buffer_write(GPU_current_buffer);
    } else {
       GPU_buffers[GPU_current_buffer].buffer = NULL;
    }

    return GPU_current_buffer;
}

int             GPU::buffer_write(int buffer_id)
{
    cl_event buffer_event;
    cl_ulong buffer_write_start, buffer_write_finish;
    OpenCL_Check_Error(clEnqueueWriteBuffer( GPU_queue, GPU_buffers[buffer_id].buffer, CL_TRUE, 0, GPU_buffers[buffer_id].size_in_bytes, GPU_buffers[buffer_id].host_ptr, 0, NULL, &buffer_event),"clEnqueueWriteBuffer failed");
    if (GPU_debug.profiling){
        OpenCL_Check_Error(clWaitForEvents(1, &buffer_event),"clWaitForEvents failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &buffer_write_finish, 0),"clGetEventProfilingInfo failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &buffer_write_start,  0),"clGetEventProfilingInfo failed");
        double elapsed_time = (double) (buffer_write_finish-buffer_write_start);
        GPU_buffers[buffer_id].buffer_write_elapsed_time           += elapsed_time;
        GPU_buffers[buffer_id].buffer_write_start                   = buffer_write_start;
        GPU_buffers[buffer_id].buffer_write_finish                  = buffer_write_finish;
        GPU_buffers[buffer_id].buffer_write_elapsed_time_squared   += elapsed_time*elapsed_time;
        GPU_buffers[buffer_id].buffer_write_number_of++;
    }
    return (int) GPU_error;
}

unsigned int*   GPU::buffer_map(int buffer_id)
{
    cl_uint *ptr;
    cl_event buffer_event;
    cl_ulong buffer_read_start, buffer_read_finish;
    ptr = (cl_uint *) clEnqueueMapBuffer( GPU_queue,GPU_buffers[buffer_id].buffer,CL_TRUE,CL_MAP_READ,0,GPU_buffers[buffer_id].size_in_bytes,0, NULL, &buffer_event, &GPU_error );
    OpenCL_Check_Error(GPU_error,"clEnqueueMapBuffer failed");
    if (GPU_debug.profiling){
        OpenCL_Check_Error(clWaitForEvents(1, &buffer_event),"clWaitForEvents failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &buffer_read_finish, 0),"clGetEventProfilingInfo failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &buffer_read_start,  0),"clGetEventProfilingInfo failed");
        GPU_buffers[buffer_id].buffer_read_elapsed_time   += (double) (buffer_read_finish-buffer_read_start);
        GPU_buffers[buffer_id].buffer_read_start           = buffer_read_start;
        GPU_buffers[buffer_id].buffer_read_finish          = buffer_read_finish;
        GPU_buffers[buffer_id].buffer_read_number_of++;
    }
    GPU_buffers[buffer_id].mapped_ptr = ptr;
    return ptr;
}

cl_float4*      GPU::buffer_map_float4(int buffer_id)
{
    cl_float4* ptr;
    cl_event buffer_event;
    cl_ulong buffer_read_start, buffer_read_finish;
    ptr = (cl_float4*) clEnqueueMapBuffer( GPU_queue,GPU_buffers[buffer_id].buffer,CL_TRUE,CL_MAP_READ,0,GPU_buffers[buffer_id].size_in_bytes,0, NULL, &buffer_event, &GPU_error );
    OpenCL_Check_Error(GPU_error,"clEnqueueMapBuffer failed");
    if (GPU_debug.profiling){
        OpenCL_Check_Error(clWaitForEvents(1, &buffer_event),"clWaitForEvents failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &buffer_read_finish, 0),"clGetEventProfilingInfo failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &buffer_read_start,  0),"clGetEventProfilingInfo failed");
        GPU_buffers[buffer_id].buffer_read_elapsed_time   += (double) (buffer_read_finish-buffer_read_start);
        GPU_buffers[buffer_id].buffer_read_start           = buffer_read_start;
        GPU_buffers[buffer_id].buffer_read_finish          = buffer_read_finish;
        GPU_buffers[buffer_id].buffer_read_number_of++;
    }
    GPU_buffers[buffer_id].mapped_ptr = (unsigned int*) ptr;
    return ptr;
}

cl_float4*      GPU::buffer_read_float4(int buffer_id)
{
    cl_event buffer_event;
    cl_ulong buffer_read_start, buffer_read_finish;
    cl_float4* ptr = (cl_float4*) GPU_buffers[buffer_id].mapped_ptr;
    GPU_error = clEnqueueReadBuffer(GPU_queue,GPU_buffers[buffer_id].buffer,CL_TRUE,0,GPU_buffers[buffer_id].size_in_bytes,ptr,0,NULL,&buffer_event);
    OpenCL_Check_Error(GPU_error,"clEnqueueReadBuffer failed");
    if (GPU_debug.profiling){
        OpenCL_Check_Error(clWaitForEvents(1, &buffer_event),"clWaitForEvents failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &buffer_read_finish, 0),"clGetEventProfilingInfo failed");
        OpenCL_Check_Error(clGetEventProfilingInfo(buffer_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &buffer_read_start,  0),"clGetEventProfilingInfo failed");
        GPU_buffers[buffer_id].buffer_read_elapsed_time   += (double) (buffer_read_finish-buffer_read_start);
        GPU_buffers[buffer_id].buffer_read_start           = buffer_read_start;
        GPU_buffers[buffer_id].buffer_read_finish          = buffer_read_finish;
        GPU_buffers[buffer_id].buffer_read_number_of++;
    }
    return ptr;
}

int             GPU::buffer_kill(int buffer_id)
{
    cl_int result = CL_SUCCESS;
    if (GPU_buffers[buffer_id].buffer) {
        result = clReleaseMemObject(GPU_buffers[buffer_id].buffer);
        GPU_buffers[buffer_id].buffer = NULL;
    }
    else result = CL_INVALID_MEM_OBJECT;
    return result;
}

GPU::GPU_time_deviation  GPU::buffer_write_get_time(int buffer_id){
    GPU_time_deviation execution_time = time_get_deviation(GPU_buffers[buffer_id].buffer_write_elapsed_time,GPU_buffers[buffer_id].buffer_write_elapsed_time_squared,GPU_buffers[buffer_id].buffer_write_number_of);
    return execution_time;
}

GPU::GPU_time_deviation  GPU::buffer_read_get_time(int buffer_id){
    GPU_time_deviation execution_time = time_get_deviation(GPU_buffers[buffer_id].buffer_read_elapsed_time,GPU_buffers[buffer_id].buffer_read_elapsed_time_squared,GPU_buffers[buffer_id].buffer_read_number_of);
    return execution_time;
}

unsigned int    GPU::buffer_size_align(unsigned int size)
{
    // align buffer size to memory_align_factor
    unsigned int result = buffer_size_align(size,(unsigned int) GPU_info.max_workgroup_size);
    return result;
}

unsigned int    GPU::buffer_size_align(unsigned int size,unsigned int step)
{
    return (1 + (size-1) / step) * step;
}

// ___ print ______________________________________________________________________________________
void            GPU::print_available_hardware(void)
{
    printf("OpenCL platforms available: %u\n",  GPU_platforms_number);
    printf("OpenCL devices available:   %u\n",  GPU_total_devices_number);
    printf("Active OpenCL platform: %s\n",      platform_get_name(GPU_platform));
    printf("Active OpenCL device:   %s\n\n",    GPU_info.device_name);
}

void            GPU::print_stage(const char * stage){
                if (GPU_debug.show_stage) printf(">>> Runtime stage >>> %s\n",stage);
}

void            GPU::print_memory_utilized(void){
    double result = 0;  // in bytes
    for (int i=1; i<GPU_current_buffer; i++) {
            result += (double) GPU_buffers[i].size_in_bytes;
    }
    printf("Used GPU memory, MB: %f\n",result/1024./1024.);
}

int             GPU::print_time_detailed(void){
    GPU_time_deviation elapsed_time, elapsed_time_read, elapsed_time_write;
    printf("--------------------------------------------------------\n");
    for (int i=1; i<=GPU_current_kernel; i++){
        elapsed_time = kernel_get_execution_time(i);
        if(GPU_debug.brief_report) {
            printf("[%2u] kernel \"%s\" (N=%u): %f ms\n",i,GPU_kernels[i].kernel_name,GPU_kernels[i].kernel_number_of_starts,elapsed_time.mean*1.E-6f);
        } else {
            printf("[%2u] kernel \"%s\" executed %u times:\n\t execution time: %f (+/-%f) ms\n",i,GPU_kernels[i].kernel_name,GPU_kernels[i].kernel_number_of_starts,elapsed_time.mean*1.E-6f,elapsed_time.deviation*1.E-6f);
        }
    }
    printf("--------------------------------------------------------\n");

    for (int i=1; i<=GPU_current_buffer; i++) {
        elapsed_time_write = buffer_write_get_time(i);
        elapsed_time_read  = buffer_read_get_time(i);
        double buffer_size = ((double) GPU_buffers[i].size_in_bytes);
        if(GPU_debug.brief_report) {
            printf("[%2u] %f Mbytes:\n",i,buffer_size/1024.0/1024.);
        } else {
            printf("[%2u] buffer size=%f Mbytes:\n",i,buffer_size/1024.0/1024.);
        }
        if ((GPU_buffers[i].buffer_write_number_of)>0){
            if(GPU_debug.brief_report) {
                printf("\t\tH>D (N=%u): %f Mbytes/sec\n",GPU_buffers[i].buffer_write_number_of,time_megabytes_per_second(elapsed_time_write.mean,buffer_size));
            } else {
                printf("\t host to device %u times\n\t\t elapsed time: %f (+/-%f) ms (%f Mbytes/sec)\n",GPU_buffers[i].buffer_write_number_of,elapsed_time_write.mean*1.E-6f,elapsed_time_write.deviation*1.E-6f,time_megabytes_per_second(elapsed_time_write.mean,buffer_size));
            }
        }
        if (GPU_buffers[i].buffer_read_number_of>0){
            if(GPU_debug.brief_report) {
                printf("\t\tD>H (N=%u): %f Mbytes/sec\n",GPU_buffers[i].buffer_read_number_of,time_megabytes_per_second(elapsed_time_read.mean,buffer_size));
            } else {
                printf("\t device to host %u times:\n\t\t elapsed time: %f (+/-%f) ms (%f Mbytes/sec)\n",GPU_buffers[i].buffer_read_number_of,elapsed_time_read.mean*1.E-6f,elapsed_time_read.deviation*1.E-6f,time_megabytes_per_second(elapsed_time_read.mean,buffer_size));
            }
}
    }
    printf("--------------------------------------------------------\n");


    return 0;
}

int             GPU::print_mapped_buffer_uint(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %u\n", i, ptr[i]);
    }
    return 0;
}

int             GPU::print_mapped_buffer_uint2(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %u \t %u\n", i, ptr[2*i], ptr[2*i+1]);
    }
    return 0;
}

int             GPU::print_mapped_buffer_uint3(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %u \t %u \t %u\n", i, ptr[3*i], ptr[3*i+1], ptr[3*i+2]);
    }
    return 0;
}

int             GPU::print_mapped_buffer_uint4(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %u \t %u \t %u \t %u\n", i, ptr[4*i], ptr[4*i+1], ptr[4*i+2], ptr[4*i+3]);
    }
    return 0;
}

int             GPU::print_mapped_buffer_float(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f\n", i, convert_to_float(ptr[i]));
    }
    return 0;
}

int             GPU::print_mapped_buffer_float2(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f \t %f\n", i, convert_to_float(ptr[2*i]), convert_to_float(ptr[2*i+1]));
    }
    return 0;
}

int             GPU::print_mapped_buffer_float3(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f \t %f \t %f\n", i, convert_to_float(ptr[3*i]), convert_to_float(ptr[3*i+1]), convert_to_float(ptr[3*i+2]));
    }
    return 0;
}

int             GPU::print_mapped_buffer_float4(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f \t %f \t %f \t %f\n", i, convert_to_float(ptr[4*i]), convert_to_float(ptr[4*i+1]), convert_to_float(ptr[4*i+2]), convert_to_float(ptr[4*i+3]));
    }
    return 0;
}

int             GPU::print_mapped_buffer_float4(int buffer_id,unsigned int number_of_elements,unsigned int offset){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f \t %f \t %f \t %f\n", i, convert_to_float(ptr[4*(i + offset)]), convert_to_float(ptr[4*(i + offset)+1]), convert_to_float(ptr[4*(i + offset)+2]), convert_to_float(ptr[4*(i + offset)+3]));
    }
    return 0;
}

int             GPU::print_mapped_buffer_double(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f\n", i, convert_to_double(ptr[i*2],ptr[i*2+1]));
    }
    return 0;
}

int             GPU::print_mapped_buffer_double2(int buffer_id,unsigned int number_of_elements){
    unsigned int* ptr = GPU_buffers[buffer_id].mapped_ptr;
    if (ptr==NULL) return GPU_error_no_buffer;
    for (unsigned int i = 0; i<number_of_elements; i++){
        printf("[%4u] %f  \t  %f\n", i, convert_to_double(ptr[4*i],ptr[4*i+1]), convert_to_double(ptr[4*i+2],ptr[4*i+3]));
    }
    return 0;
}

// ___ timers _____________________________________________________________________________________
int             GPU::start_timer_CPU(void)
{
    CPU_current_timer_id++;
    int id = CPU_current_timer_id;
    if (id>=CPU_timers) id = 0;
    CPU_timer[id] = clock();
    return id;
}

int             GPU::start_timer_CPU(int timer)
{
    int id = timer;
    if (id>=CPU_timers) id = 0;
    CPU_timer[id] = clock();
    return id;
}

int             GPU::stop_timer_CPU(int timer)
{
    int id = timer;
    if (id>=CPU_timers) id = 0;
    CPU_timer[id] = clock() - CPU_timer[id];
    return CPU_timer[id];
}

double          GPU::timer_in_seconds_CPU(int timer)
{
    int id = timer;
    if (id >= CPU_timers) id = 0;
    return CPU_timer[id] * 1.0 / CLOCKS_PER_SEC;
}

double          GPU::get_timer_CPU(int timer)
{
    int id = timer;
    if (id >= CPU_timers) id = 0;
    return (clock() - CPU_timer[id]) * 1.0 / CLOCKS_PER_SEC;
}

char*           GPU::get_current_datetime(void){
    time_t tim;
    time(&tim);
    char* result = (char*) calloc(256,sizeof(char));
#ifdef _WIN32
    ctime_s(result, 26, &tim);
#else
    sprintf(result,"%s",ctime((const time_t*) &tim));
#endif
    return result;
}

// ___ MD5 ________________________________________________________________________________________
inline unsigned int GPU::MD5_rol(unsigned int x, int s) {
  return ((x << s) | (x >> (32-s)));
}

void            GPU::MD5_init(void){
    MD5_finalized = false;
    MD5_bytes[0] = 0;
    MD5_bytes[1] = 0;

    MD5_state[0] = 0x67452301;
    MD5_state[1] = 0xefcdab89;
    MD5_state[2] = 0x98badcfe;
    MD5_state[3] = 0x10325476;

    MD5_s[ 0] = 7;  MD5_s[ 1] = 12; MD5_s[ 2] = 17; MD5_s[ 3] = 22;
    MD5_s[16] = 5;  MD5_s[17] = 9;  MD5_s[18] = 14; MD5_s[19] = 20;
    MD5_s[32] = 4;  MD5_s[33] = 11; MD5_s[34] = 16; MD5_s[35] = 23;
    MD5_s[48] = 6;  MD5_s[49] = 10; MD5_s[50] = 15; MD5_s[51] = 21;

    for (int i=0; i<4; i++)
    for (int j=1; j<4; j++) {
        MD5_s[16 * i + 4 * j + 0] = MD5_s[16 * i    ];
        MD5_s[16 * i + 4 * j + 1] = MD5_s[16 * i + 1];
        MD5_s[16 * i + 4 * j + 2] = MD5_s[16 * i + 2];
        MD5_s[16 * i + 4 * j + 3] = MD5_s[16 * i + 3];

    }

    for (int i=0; i<64; i++)
        MD5_T[i] = (unsigned int) floor(abs(sin((double) i + 1)) * pow(2.0,32));

    memset(MD5_buffer, 0, sizeof(MD5_buffer));
}

void            GPU::MD5_finalize(void) {
    unsigned char MD5_padding[64];
    MD5_padding[0] = 128;
    for (int i=1; i<64; i++) MD5_padding[i]=0;

    if (!MD5_finalized) {
        unsigned char bits[8];
        MD5_setword(MD5_bytes, bits, 8);

    unsigned int index = MD5_bytes[0] / 8 % 64;
    unsigned int padLen = (index < 56) ? (56 - index) : (120 - index);
    MD5_update(MD5_padding, padLen);
    MD5_update(bits, 8);
    MD5_setword(MD5_state, MD5_result, 16);

    // Zeroize sensitive information.
    memset(MD5_buffer, 0, sizeof(MD5_buffer));
    memset(MD5_bytes, 0, sizeof(MD5_bytes));

    MD5_finalized=true;
  }
}

void            GPU::MD5_getword(const unsigned char* buffer, unsigned int* x, unsigned int len) {
    unsigned int i = 0;
    for (unsigned int j = 0; j < len; j += 4) {
        x[i] = ((unsigned int) buffer[j    ]) |
              (((unsigned int) buffer[j + 1]) <<  8) |
              (((unsigned int) buffer[j + 2]) << 16) |
              (((unsigned int) buffer[j + 3]) << 24);
        i++;
    }
}

void            GPU::MD5_setword(const unsigned int* x, unsigned char* buffer, unsigned int len) {
    unsigned int i = 0;
    for (unsigned int j = 0; j < len; j += 4) {
        buffer[j    ] = (x[i]      ) & 0xff;
        buffer[j + 1] = (x[i] >>  8) & 0xff;
        buffer[j + 2] = (x[i] >> 16) & 0xff;
        buffer[j + 3] = (x[i] >> 24) & 0xff;
        i++;
    }
}

void            GPU::MD5_step(const unsigned char* MD5_block) {
    unsigned int a = MD5_state[0];
    unsigned int b = MD5_state[1];
    unsigned int c = MD5_state[2];
    unsigned int d = MD5_state[3];

    unsigned int f,g,temp;
    unsigned int w[MD5_blocksize];

    MD5_getword(MD5_block, w, MD5_blocksize);

    for(int i=0; i<MD5_blocksize; i++){
        if (i < 16) {
            f = (b & c) | ((~b) & d);
            g = i;
        } else if ((16 <= i) && (i < 32)) {
            f = (d & b) | ((~d) & c);
            g = (5 * i + 1) % 16;
        } else if ((32 <= i) && (i < 48)) {
            f = b ^ c ^ d;
            g = (3 * i + 5) % 16;
        } else if ((48 <= i) && (i < 64)) {
            f = c ^ ((~d) | b);
            g = (7 * i) % 16;
        } else {
            f = 0;
            g = 0;
        }
        temp = d;
        d = c;
        c = b;
        b = b + MD5_rol((a + f + MD5_T[i] + w[g]) , MD5_s[i]);
        a = temp;
    }

    MD5_state[0] += a;
    MD5_state[1] += b;
    MD5_state[2] += c;
    MD5_state[3] += d;

    memset(w, 0, sizeof w);
}

void            GPU::MD5_update(const unsigned char* input, unsigned int len) {
    unsigned int index = MD5_bytes[0] / 8 % MD5_blocksize;

    if ((MD5_bytes[0] += (len << 3)) < (len << 3)) MD5_bytes[1]++;
    MD5_bytes[1] += (len >> 29);

    unsigned int firstpart = 64 - index;
    unsigned int i;

    if (len >= firstpart) {
        memcpy(&MD5_buffer[index], input, firstpart);
        MD5_step(MD5_buffer);

        for (i = firstpart; i + MD5_blocksize <= len; i += MD5_blocksize) MD5_step(&input[i]);
        index = 0;
    } else
        i = 0;

    memcpy(&MD5_buffer[index], &input[i], len-i);
}

void            GPU::MD5_update(const char* input, unsigned int len)
{
    MD5_update((const unsigned char*) input, len);
}

char*           GPU::MD5_getresult(void){
  char* buf = (char*) calloc(34,sizeof(int));
  if (MD5_finalized) {
    for (int i=0; i<16; i++)
        sprintf_s(buf+i*2, sizeof(buf), "%02x", MD5_result[i]);
  }
        buf[32]=0;

  return buf;
}

char*           GPU::MD5(const char* str) {
    unsigned int len = (unsigned int) strlen(str);

    MD5_init();
    MD5_update(str, len);
    MD5_finalize();

    return MD5_getresult();
}

// ___ other ______________________________________________________________________________________
GPU::GPU_time_deviation  GPU::time_get_deviation(double elapsed_time, double elapsed_time_squared, int number)
{
    GPU_time_deviation execution_time;
    execution_time.mean     = elapsed_time;
    execution_time.number   = number;
    execution_time.deviation = number < 2 ? 0.0 : sqrt(abs(elapsed_time_squared - pow(elapsed_time,2)) / ((double) number));
    return execution_time;
}

double          GPU::time_gigabytes_per_second(double elapsed_time,double size){
    if (elapsed_time==0.0)
        return 0.0;
    else
        return (size/elapsed_time/(1.024*1.024*1.024)); // converting bytes->gigabytes
}

double          GPU::time_megabytes_per_second(double elapsed_time,double size){
    if (elapsed_time==0.0)
        return 0.0;
    else
        return (size/elapsed_time/(1.024*1.024*1.024)*1024.0); // converting bytes->megabytes
}

int             GPU::str_char_replace(char* str, char search, char replace){
      int result = 0;
        for (int i = 0; i<strlen(str); i++)
            if (str[i]==search){
                str[i] = replace;
                result++;
            }
      return result;
}

GPU::GPU_init_parameters*    GPU::get_init_file(char finitf[])
{
    FILE *stream;
    int struc_quant = 16;
    int struc_length = 0;

    GPU_init_parameters* result = (GPU_init_parameters*) calloc(struc_quant, sizeof(GPU_init_parameters));

    char line[250];
    char buffer[250];
    int j,j2;
    char Variable[250];
    char txtVal[250];
    int iVarVal;
    double fVarVal;

    j = sprintf_s(buffer,sizeof(buffer),"%s",finitf);

    fopen_s(&stream,buffer,"r");
    if(stream)
    {
    while(fgets( line, 250, stream ) != NULL)
      {
        int istart1 = 0;
        char ch=line[istart1];
        while((ch==' ')||(ch=='\t')) ch=line[++istart1];

        unsigned int istart2 = (unsigned int) strcspn(line,"=");
        unsigned int istart3 = istart2;
        if (istart2<strlen(line)){
            ch=line[--istart2];
            while((ch==' ')||(ch=='\t')) ch=line[--istart2];

            int j2=0;
            for (unsigned int j=istart1; j<=istart2; j++)
            {
                Variable[j2]=line[j];
                j2++;
            }
            Variable[j2]=0;
        }
        unsigned int ifinish = (unsigned int) strcspn(line,"#");
        ch=line[--ifinish];
        while((ch==' ')||(ch=='\t')||(ch=='\n')) ch=line[--ifinish];

        ch=line[++istart3];
        while((ch==' ')||(ch=='\t')) ch=line[++istart3];

            int j3=0;
            for (unsigned int j=istart3; j<=ifinish; j++)
            {
                txtVal[j3]=line[j];
                j3++;
            }
            txtVal[j3]=0;

            sscanf_s(txtVal,"%d", &iVarVal);
            sscanf_s(txtVal,"%lf", &fVarVal);

            if (struc_length>0) result[struc_length-1].final = false;
            j2 = sprintf_s(result[struc_length].Variable,sizeof(result[struc_length].Variable),"%s",Variable);
            j2 = sprintf_s(result[struc_length].txtVarVal,sizeof(result[struc_length].Variable),"%s",txtVal);
            result[struc_length].iVarVal = iVarVal;
            result[struc_length].fVarVal = fVarVal;
            result[struc_length].final   = true;

            struc_length++;
            if (struc_length % struc_quant == 0)
                result = (GPU_init_parameters*) realloc(result, (struc_quant + struc_length) * sizeof(GPU_init_parameters));
      }
      fclose( stream );
    }

    if (struc_length == 0) free(result);

    return result;
}

void            GPU::make_start_file(char* path)
{
        FILE *stream;

        char buffer[FNAME_MAX_LENGTH*2];
        int j = sprintf_s(buffer  ,FNAME_MAX_LENGTH*2,"%s",path);
        j += sprintf_s(buffer+j,FNAME_MAX_LENGTH*2-j,"%s","finish.txt");
        bool flag = true;

        while (flag)
        {
           fopen_s(&stream,buffer,"r");
           if (stream)
            {
               fclose(stream);
            }
           else
            {
               flag=false;
            }
           Sleep(500);
        }

        j = sprintf_s(buffer  ,FNAME_MAX_LENGTH*2,"%s",path);
        j += sprintf_s(buffer+j,FNAME_MAX_LENGTH*2-j,"%s","start.txt");

        fopen_s(&stream,buffer,"w+");
        if(stream)
         {
            fprintf(stream, "%s\n","task started");
         }

   if( stream)
    {
       if ( fclose( stream ) )
        {
           printf( "The file was not closed!\n" );
        }
    }
}

void            GPU::make_finish_file(char* path)
{
    if(path) {
       FILE *stream;

       char buffer[FNAME_MAX_LENGTH*2];
       int j = sprintf_s(buffer  ,FNAME_MAX_LENGTH*2,"%s",path);
       j += sprintf_s(buffer+j,FNAME_MAX_LENGTH*2-j,"%s","finish.txt");

       fopen_s(&stream,buffer,"w+");
          if(stream)
           {
              fprintf(stream, "%s\n","task done");
           }

          if(stream)
           {
              if ( fclose( stream ) )
               {
                  printf( "The file was not closed!\n" );
               }
           }
    }
}

int             GPU::inf_file_delete(int index){
        int err;
        char* buffer_inf = (char*) calloc(FNAME_MAX_LENGTH,sizeof(char));

        // kill .inf-file
        int j  = sprintf_s(buffer_inf,FNAME_MAX_LENGTH,"program%u.inf",index);
        err = remove(buffer_inf);
        if (err) {delete[] buffer_inf; return err;}
        j = sprintf_s(buffer_inf,FNAME_MAX_LENGTH,"program%u.bin",index);
        // kill .bin-file
        err = remove(buffer_inf);
        delete[] buffer_inf;
        return err;
}

int             GPU::inf_file_rename(int index_old,int index_new){
        int j;
        int err;
        char* buffer_old = (char*) calloc(FNAME_MAX_LENGTH,sizeof(char));
        char* buffer_new = (char*) calloc(FNAME_MAX_LENGTH,sizeof(char));

            j = sprintf_s(buffer_old,FNAME_MAX_LENGTH,"program%u.inf",index_old);
            j = sprintf_s(buffer_new,FNAME_MAX_LENGTH,"program%u.inf",index_new);
        // rename .inf-file
        err = rename(buffer_old,buffer_new);
        if (err) {delete[] buffer_old; delete[] buffer_new; return err;}

            j = sprintf_s(buffer_old,FNAME_MAX_LENGTH,"program%u.bin",index_old);
            j = sprintf_s(buffer_new,FNAME_MAX_LENGTH,"program%u.bin",index_new);
        // rename .bin-file
        err = rename(buffer_old,buffer_new);
        delete[] buffer_old;
        delete[] buffer_new;

        return err;
}

int             GPU::inf_file_get_max_n(void){
        char* buffer_inf = (char*) calloc(FNAME_MAX_LENGTH,sizeof(char));
        int j = 0;

        int result = 1;
        bool flag = true;
        while(flag){
            j = sprintf_s(buffer_inf,FNAME_MAX_LENGTH,"program%u.inf",result);
            flag = is_file_exist(buffer_inf);
            result++;
        }
        result -= 2;
        free(buffer_inf);
        return result;
}

bool            GPU::is_file_exist(char* path){
        FILE * is_file;
        bool flag = true;
         fopen_s(&is_file,path,"r");
            if (!is_file) flag = false;
            else fclose(is_file);
        return flag;
}

}
