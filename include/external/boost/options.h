/** @file options.h
 *  @brief define command line options using boost/program_options
 *
 *  Includes the definition of all command line options with description.
 *  Call program with --help to see these options.
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_OPTIONS_BOOST_H
#define	PASC_OPTIONS_BOOST_H

#include <boost/program_options.hpp>

namespace pascinference {

/** @brief add all options of common library files using boost/program_options
*
* @param description boost program options instance
* @param console_nmb_cols number of columns in console
*/
extern void add_options(boost::program_options::options_description *description, int console_nmb_cols);


}



#endif
