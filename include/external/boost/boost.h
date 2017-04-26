#ifndef PASCINFERENCE_BOOST_H
#define	PASCINFERENCE_BOOST_H

#include "pascinference.h"

#ifdef USE_CUDA
	#undef _GLIBCXX_ATOMIC_BUILTINS
#endif

/* load console parameters with boost */
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#endif
