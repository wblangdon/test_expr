/*WBL 14 Feb 2016 "$Revision: 1.3 $"
Based on OpenCV version: 3.1.0-dev Apr 20  2016
opencv_ssbse2016/david/davidrobertwhite-gicow-a9fc22557c67/opencv/modules/core/src/system.cpp
  Compiles ok with gcc version 4.8.5 20150623, gcc -c time.c
*/
#include <time.h>
#include "tick.h"

/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009, Willow Garage Inc., all rights reserved.
// Copyright (C) 2015, Itseez Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

int64 getTickCount(void)
{
#if defined WIN32 || defined _WIN32 || defined WINCE
    LARGE_INTEGER counter;
    QueryPerformanceCounter( &counter );
    return (int64)counter.QuadPart;
#elif defined __linux || defined __linux__
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return (int64)tp.tv_sec*1000000000 + tp.tv_nsec;
#elif defined __MACH__ && defined __APPLE__
    return (int64)mach_absolute_time();
#else
    struct timeval tv;
    struct timezone tz;
    gettimeofday( &tv, &tz );
    return (int64)tv.tv_sec*1000000 + tv.tv_usec;
#endif
}

double getTickFrequency(void)
{
#if defined WIN32 || defined _WIN32 || defined WINCE
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return (double)freq.QuadPart;
#elif defined __linux || defined __linux__
    return 1e9;
#elif defined __MACH__ && defined __APPLE__
    static double freq = 0;
    if( freq == 0 )
    {
        mach_timebase_info_data_t sTimebaseInfo;
        mach_timebase_info(&sTimebaseInfo);
        freq = sTimebaseInfo.denom*1e9/sTimebaseInfo.numer;
    }
    return freq;
#else
    return 1e6;
#endif
}

#if defined __GNUC__ && (defined __i386__ || defined __x86_64__ || defined __ppc__)
#if defined(__i386__)

int64 getCPUTickCount(void)
{
    int64 x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#elif defined(__x86_64__)

int64 getCPUTickCount(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return (int64)lo | ((int64)hi << 32);
}

#elif defined(__ppc__)

int64 getCPUTickCount(void)
{
    int64 result = 0;
    unsigned upper, lower, tmp;
    __asm__ volatile(
                     "0:                  \n"
                     "\tmftbu   %0           \n"
                     "\tmftb    %1           \n"
                     "\tmftbu   %2           \n"
                     "\tcmpw    %2,%0        \n"
                     "\tbne     0b         \n"
                     : "=r"(upper),"=r"(lower),"=r"(tmp)
                     );
    return lower | ((int64)upper << 32);
}

#else

#error "RDTSC not defined"

#endif

#elif defined _MSC_VER && defined WIN32 && defined _M_IX86

int64 getCPUTickCount(void)
{
    __asm _emit 0x0f;
    __asm _emit 0x31;
}

#else

//#ifdef HAVE_IPP
//int64 getCPUTickCount(void)
//{
//    return ippGetCpuClocks();
//}
//#else
int64 getCPUTickCount(void)
{
    return getTickCount();
}
//#endif

#endif

