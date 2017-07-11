/** @file common.h
 *  @brief commonly used stuff
 *
 *  This file includes commonly used classes and functions.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

#define EXPORT_SAVEVTK true /* export solution to VTK */ //TODO: old and unnecessary?
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */ //TODO: old and unnecessary?

#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iostream>
#include <iomanip>
#include <typeinfo> 
#include <cmath>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <stack>
#include <limits>

#include "general/common/initialize.h"
#include "general/common/timer.h"
#include "general/common/memorycheck.h"
#include "general/common/powercheck.h"
#include "general/common/globalmanager.h"
#include "general/common/consoleoutput.h"
#include "general/common/consoleinput.h"
#include "general/common/logging.h"
#include "general/common/mvnrnd.h"
#include "general/common/shortinfo.h"
#include "general/common/decomposition.h"


#endif
