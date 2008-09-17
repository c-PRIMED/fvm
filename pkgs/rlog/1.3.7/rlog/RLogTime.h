/*****************************************************************************
 * Author:   Vadim Zeitlin <vadim@wxwidgets.org>
 *
 *****************************************************************************
 * Copyright (c) 2004 Vadim Zeitlin
 *
 * This library is free software; you can distribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (LGPL), as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the LGPL in the file COPYING for more
 * details.
 *
 */

#ifndef _rlog_time_incl
#define _rlog_time_incl

#include <rlog/common.h>

/*! @file time.h
    @brief Defines functions for getting current time and benchmarking.
*/

#ifdef _WIN32

#include <windows.h>

#define HAVE_QUERYPERFORMANCECOUNTER 1

typedef __int64 rlog_time_interval;

#if HAVE_QUERYPERFORMANCECOUNTER

typedef LARGE_INTEGER rlog_time;

#define RLOG_TIME_UNIT "clock cycles"

inline
void rlog_get_time(rlog_time *pt)
{
    QueryPerformanceCounter(pt);
}

inline
rlog_time_interval rlog_time_diff(const rlog_time& end, const rlog_time& start)
{
    long long llEnd, llStart;
    memcpy(&llEnd, &end, sizeof(long long));
    memcpy(&llStart, &start, sizeof(long long));
    return llEnd - llStart;
}

#else // !HAVE_QUERYPERFORMANCECOUNTER

typedef FILETIME rlog_time;

#define RLOG_TIME_UNIT "usec"

inline
void rlog_get_time(rlog_time *pt)
{
    GetSystemTimeAsFileTime(pt);
}

inline
rlog_time_interval rlog_time_diff(const rlog_time& end, const rlog_time& start)
{
    ULONGLONG ullEnd, ullStart;
    memcpy(&ullEnd, &end, sizeof(ULONGLONG));
    memcpy(&ullStart, &start, sizeof(ULONGLONG));
    return 10*(ullEnd - ullStart);
}

#endif // HAVE_QUERYPERFORMANCECOUNTER

inline
void sleep(int seconds)
{
    ::Sleep(seconds * 1000);
}

#else // Unix

#include <sys/time.h>
#include <unistd.h> // for sleep()

#if RLOG_TIME_TSC

#include <stdint.h>

typedef uint64_t rlog_time;
typedef int64_t rlog_time_interval;

#define RLOG_TIME_UNIT "clock cycles"

inline void rlog_get_time(uint64_t *pt)
{
    asm volatile("RDTSC" : "=A" (*pt));
}

inline
rlog_time_interval rlog_time_diff( const rlog_time &end, const rlog_time &start )
{
    return end - start;
}

#else // !HAVE_TSC

#include <unistd.h>

typedef timeval rlog_time;
typedef long rlog_time_interval;

#define RLOG_TIME_UNIT "usec"

inline
void rlog_get_time(rlog_time *pt)
{
    gettimeofday( pt, 0 );
}

inline
rlog_time_interval rlog_time_diff( const rlog_time &end, const rlog_time &start )
{
    return (end.tv_sec - start.tv_sec) * 1000 * 1000 + 
	(end.tv_usec - start.tv_usec);
}

#endif // HAVE_TSC/!HAVE_TSC

#endif // Win32/Unix

#endif // _rlog_time_incl
