#ifndef _MACROPARAMETERS_H_
#define _MACROPARAMETERS_H_


#include "Field.h"

struct MacroParameters
{
  MacroParameters(const string baseName);
  Field velocity;
  Field pressure;
  Field temperature;
  Field viscosity;
  Field density;
 
};

#endif
