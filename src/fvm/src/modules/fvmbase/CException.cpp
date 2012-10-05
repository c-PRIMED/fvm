// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "CException.h"

#ifdef __GNUC__
#include <execinfo.h> 
#include <cxxabi.h>
#include <stdlib.h>
#endif

using namespace std;



CException::CException(const string& what) :
  runtime_error(what)
{
#if 1
#ifdef __GNUC__
    void *traces[100];
    cout << "Backtrace... " << endl;
    int traceCount = backtrace(traces,20);
    char **traceStrings = backtrace_symbols(traces, traceCount);

    if(traceStrings != NULL)
    {
        for(int x = 0; x < traceCount; x++)
        {
            if (traceStrings[x] == NULL) { break; }
            int status;
            string s(traceStrings[x]);
            size_t funcNameBegin = s.find("(");
            size_t funcNameEnd = s.find("+");
            if (funcNameBegin != string::npos &&
                funcNameEnd != string::npos)
            {
                string funcNameMangled(s.begin()+funcNameBegin+1,
                                       s.begin()+funcNameEnd);
                
                string fileName(s.begin(),s.begin()+funcNameBegin);
                
                string offset(s.begin() + funcNameEnd, s.end());
                
                char *realName = abi::__cxa_demangle(funcNameMangled.c_str(), 0, 0, &status);
                
                std::cout << x << " :  " << realName   << endl;
                free(realName);
            }
            else
                std::cout << x << " :  " << s   << endl;
        }
    }
    cout << "Backtrace...end " << endl;
    cout << flush;
    free(traceStrings);
#endif
#endif
    
}

