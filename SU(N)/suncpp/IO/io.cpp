/******************************************************************************
 * @file     io.cpp
 * @author   Vadim Demchik <vadimdi@yahoo.com>,
 * @author   Natalia Kolomoyets <rknv7@mail.ru>
 * @version  1.6
 *
 * @brief    [QCDGPU]
 *           Defines procedures for work with files
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

#include "io.h"

#define FNAME_MAX_LENGTH    128

init_parameters*    get_init_file(char finitf[]){
    FILE *stream;
    int struc_quant = 16;
    int struc_length = 0;

    init_parameters* result = (init_parameters*) calloc(struc_quant, sizeof(init_parameters));
    char line[LENGTH];
    char buffer[LENGTH];
    int j,j2;
    char Variable[LENGTH];
    char txtVal[LENGTH];
    int iVarVal;
    double fVarVal;

    j = sprintf_s(buffer,sizeof(buffer),"%s",finitf);

    fopen_s(&stream,buffer,"r");
    if(stream)
    {
    while(fgets( line, LENGTH, stream ) != NULL)
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
                result = (init_parameters*) realloc(result, (struc_quant + struc_length) * sizeof(init_parameters));
      }
      fclose( stream );
    }

    if (struc_length == 0) free(result);
    return result;
}

void make_start_file(char* path){
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

void make_finish_file(char* path){
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
