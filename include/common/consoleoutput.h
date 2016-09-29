#ifndef PASC_COMMON_CONSOLEOUTPUT_H
#define	PASC_COMMON_CONSOLEOUTPUT_H

namespace pascinference {
namespace common {

/** @class OffsetClass
 *  @brief output space
 * 
 *  For manipulating with the space in the begining of output line.
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
		OffsetClass() {
			size = 0;
		}
		
		/** @brief increase the size of offset
		*
		*/
		void push(){
			size += 3;
		}

		/** @brief decrease the size of offset
		*
		*/
		void pop(){
			size -= 3;
			if(size < 0) size = 0;
		}

		friend std::ostream &operator<<(std::ostream &output, const OffsetClass &my_offset);

};

std::ostream &operator<<(std::ostream &output, const OffsetClass &my_offset){
	int i;
	for(i=0;i<my_offset.size;i++){
		output << " ";
	}
	return output;
}

OffsetClass offset; /**< global instance of output offset */


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

				ConsoleOutputBuf(std::ostream& str):output(str){
				}

				~ConsoleOutputBuf(){
				}

				virtual int sync ( ){
					if(this->rank == print_rank || this->print_rank == -1){
						std::stringstream output_string;

						/* write here also a rank of processor */
						if(GlobalManager.get_size() > 1){
							output_string << "[" << GlobalManager.get_rank() << "] " << offset << str();
						} else {
							output_string << offset << str();
						}
						str("");

//						output_string << output.rdbuf();

						if(this->print_rank == 0 || GlobalManager.get_size() <= 1){
							/* master prints */
							TRY( PetscPrintf(PETSC_COMM_WORLD, output_string.str().c_str()) );
						} else {
							/* all prints */
							TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, output_string.str().c_str()) );
						}

						output_string.str("");
						output_string.clear();

					} else {
					    str("");
					}

//					output.flush();
					return 0;
				}
				
		};

		ConsoleOutputBuf buffer; /**< instance of output buffer */

	public:

		/** @brief constructor from given output stream and rank
		*
		* @param std output stream (for example std::cout)
		*/
		ConsoleOutput(std::ostream& str, int rank = -1) : std::ostream(&buffer), buffer(str) {
			buffer.print_rank = rank;
			buffer.rank = GlobalManager.get_rank();

		}

		/** @brief increase the size of offset
		*
		*/
		void push(){
			offset.push();
		}

		/** @brief decrease the size of offset
		*
		*/
		void pop(){
			offset.pop();
		}

		void synchronize(){
			if(buffer.print_rank == -1){
				TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
			}
		}

		
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

				ConsoleOutputBuf(std::ostream& str):output(str){
				}

				~ConsoleOutputBuf(){
				}

				/** @brief set rank of this processor
				*
				*/
				void set_rank(){
					if(!rankset){
						#ifdef USE_PETSCVECTOR
							/* can be set after initialize of petsc */
							if(petscvector::PETSC_INITIALIZED){
								TRY(PetscBarrier(NULL));
						
								MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
								rankset = true;
						
							}	
						#else
							rankset = true; /* if it is not with petsc, then this is always master */
						#endif
					}
				}
				
				virtual int sync ( ){
					set_rank();
					if(this->rank == print_rank || this->print_rank == -1){
						#ifdef USE_PETSCVECTOR
							/* write here also a rank of processor */
							output << "[" << this->rank << "] " << offset << str();

						#else
							output << offset << str();
						#endif
						str("");
						output.flush();
					} else {
					    str("");
					}

					#ifdef USE_PETSCVECTOR
						TRY(PetscBarrier(NULL)); /* barrier after each cout */
					#endif

					return 0;
				}
				
		};

		ConsoleOutputBuf buffer; /**< instance of output buffer */

	public:

		/** @brief constructor from given output stream and rank
		*
		* @param std output stream (for example std::cout)
		*/
		ConsoleOutput(std::ostream& str, int rank = -1) : std::ostream(&buffer), buffer(str) {
			buffer.rankset = false;
			buffer.set_rank();
			buffer.print_rank = rank;
		}

		/** @brief increase the size of offset
		*
		*/
		void push(){
			offset.push();
		}

		/** @brief decrease the size of offset
		*
		*/
		void pop(){
			offset.pop();
		}

		void synchronize(){
		}
		
};

#endif

static ConsoleOutput coutMaster(std::cout,0); /**< instance of output console stream on master */
static ConsoleOutput coutAll(std::cout); /**< instance of output console stream */


}
} /* end of namespace */

#endif
