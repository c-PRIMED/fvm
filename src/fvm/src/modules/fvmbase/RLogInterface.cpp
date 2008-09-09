#include "RLogInterface.h"

rlog::RLogChannel* RLogInterface::ctorChannel = RLOG_CHANNEL("cdtor/ctor");
rlog::RLogChannel* RLogInterface::dtorChannel = RLOG_CHANNEL("cdtor/dtor");
rlog::RLogChannel* RLogInterface::infoChannel = RLOG_CHANNEL("info");
rlog::RLogChannel* RLogInterface::warningChannel = RLOG_CHANNEL("info/warning");
rlog::RLogChannel* RLogInterface::errorChannel = RLOG_CHANNEL("info/error");
