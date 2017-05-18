/** @file spg_fs.h
 *  @brief For manipulation with Armijo condition in GLL
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SPG_FS_H
#define	PASC_SPG_FS_H

#include "general/common/timer.h"
#include "general/common/logging.h"

namespace pascinference {
namespace solver {

/** \class SPG_fs
 *  \brief generalized Armijo condition
 *
 *  For manipulation with fs - function values for generalized Armijo condition used in SPGQP.
*/
class SPG_fs {
	private:
		int m; /**< the length of list */
		double *fs_list; /**< the list with function values */
		int last_idx;

	public: 
		/** @brief constructor
		*
		* @param new_m length of lists
		*/
		SPG_fs(int new_m);

		/** @brief destructor
		*
		*/
		~SPG_fs();

		/** @brief set all values to given one
		*
		* At the begining of computation, all values are the same, set them using this function.
		* 
		* @param fx function value
		*/
		void init(double fx);

		/** @brief return maximum value from the list
		*
		*/
		double get_max();		

		/** @brief return length of lists
		*
		*/
		int get_size();
		
		/** @brief update list - add new value and remove oldest one
		*
		*/
		void update(double new_fx);
		
		/** @brief print content of the lists
		*
		* @param output where to print
		*/
		void print(ConsoleOutput &output);
};

}
} /* end of namespace */

#endif
