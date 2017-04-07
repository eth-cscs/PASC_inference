/** @file consoleoutput.cpp
 *  @brief manipulation with print to console
 *
 *  Several utils for manipulation with console print. Defines also the offset of the print and the begin of the print with proc id.
 * 
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_CONSOLEOUTPUT_H
#define	PASC_COMMON_CONSOLEOUTPUT_H

#include <iostream>
#include <string>

#ifdef USE_PETSC
	#include "petsc.h"
	#include "petscsys.h"

	#include "external/petsc/algebra/vector/petscvector.h"
#endif

#include "common/globalmanager.h"


namespace pascinference {
namespace common {

/** @class OffsetClass
 *  @brief to control the space at the beginning of print
 * 
 *  For manipulating with the space at the begining of each output line.
 * 
 */ 
class OffsetClass {
	private:
		int size; /**< number of spaces in offset */

	public:
		/** @brief constructor
		* 
		*  Set initial size of offset to zero equal to zero.
		*
		*/
		OffsetClass();
		
		/** @brief increase the size of offset
		*
		*/
		void push();

		/** @brief decrease the size of offset
		*
		*/
		void pop();

		friend std::ostream &operator<<(std::ostream &output, const OffsetClass &my_offset);

};

extern OffsetClass offset; /**< global instance of output offset */


#ifdef USE_PETSC

/** @class ConsoleOutput
 *  @brief print only on master or on all processors
 * 
 *  Overloaded std::ostream cout function. Print based on rank.
 * 
 */ 
class ConsoleOutput : public std::ostream {
	private:

		/* Write a stream buffer that prefixes each line with Plop */
		class ConsoleOutputBuf: public std::stringbuf{
			private:
				std::ostream&   output;
			public:
				int print_rank; /**< which rank to print, -1 if all **/
				int rank; /**< rank of this process */

				ConsoleOutputBuf(std::ostream& str);

				~ConsoleOutputBuf();

				virtual int sync ();
				
		};

		ConsoleOutputBuf buffer; /**< instance of output buffer */

	public:

		/** @brief constructor from given output stream and rank
		*
		* @param std output stream (for example std::cout)
		*/
		ConsoleOutput(std::ostream& str, int rank = -1);

		/** @brief increase the size of offset
		*
		*/
		void push();

		/** @brief decrease the size of offset
		*
		*/
		void pop();

		/** @brief synchronize the output on all processes
		*
		*  In the case of PETSC call PetscSynchronizedFlush.
		*/
		void synchronize();
		
};

#else 

/** @class ConsoleOutput
 *  @brief print only on master or on all processors
 * 
 *  Overloaded std::ostream cout function. Print based on rank.
 * 
 */ 
class ConsoleOutput : public std::ostream {
	private:

		/* Write a stream buffer that prefixes each line with Rank */
		class ConsoleOutputBuf: public std::stringbuf{
			private:
				std::ostream&   output;
			public:
				int rank; /**< rank of this process */
				bool rankset; /**< the rank was already obtained */
				int print_rank; /**< which rank to print, -1 if all **/

				ConsoleOutputBuf(std::ostream& str);

				~ConsoleOutputBuf();

				/** @brief set rank of this processor
				*
				*/
				void set_rank();
				
				virtual int sync ( );
				
		};

		ConsoleOutputBuf buffer; /**< instance of output buffer */

	public:

		/** @brief constructor from given output stream and rank
		*
		* @param str output stream (for example std::cout)
		* @param rank rank of the processor which prints, if -1 then all print
		*/
		ConsoleOutput(std::ostream& str, int rank = -1);

		/** @brief increase the size of offset
		*
		*/
		void push();

		/** @brief decrease the size of offset
		*
		*/
		void pop();

		/** @brief synchronize the output on all processes
		*
		*/
		void synchronize();		
};

#endif

static ConsoleOutput coutMaster(std::cout,0); /**< instance of output console stream on master */
static ConsoleOutput coutAll(std::cout); /**< instance of output console stream */


}
} /* end of namespace */

#endif
