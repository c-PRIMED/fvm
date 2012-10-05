// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _RLOGINTERFACE_H_
#define _RLOGINTERFACE_H_

#include <rlog/rlog.h>
#include "CException.h"

using namespace std;


class RLogInterface
{
public:

  static rlog::RLogChannel *ctorChannel;
  static rlog::RLogChannel *dtorChannel;
  static rlog::RLogChannel *infoChannel;
  static rlog::RLogChannel *warningChannel;
  static rlog::RLogChannel *errorChannel;
  
};

#define logCtor() _rMessage(LOGID, RLogInterface::ctorChannel, "constructing %s (%p)", \
                               getTypeName().c_str(), this)

#define logCtorVerbose(str,...) _rMessage(LOGID, RLogInterface::ctorChannel, \
                                          "constructing %s (%p) "str,    \
                                          getTypeName().c_str(), this, ##__VA_ARGS__)

#define logDtor() _rMessage(LOGID, RLogInterface::dtorChannel, "destroying %s (%p)", \
                               getTypeName().c_str(), this)

#define logDtorVerbose(str,...) _rMessage(LOGID, RLogInterface::dtorChannel, \
                                          "destroying %s (%p) "str,    \
                                          getTypeName().c_str(), this, ##__VA_ARGS__)

#define logInfo(...) _rMessage(LOGID, RLogInterface::infoChannel, ##__VA_ARGS__)
#define logWarning(...) _rMessage(LOGID, RLogInterface::warningChannel, ##__VA_ARGS__)
#define logError(...) _rMessage(LOGID, RLogInterface::errorChannel, ##__VA_ARGS__)



#define DEFINE_TYPENAME(T)                            \
  static string getTypeName() {return T;}             \

#endif
