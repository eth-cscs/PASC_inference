/** @file globalmanager.h
 *  @brief Manipulation with information about processes - MPI_COMM size and rank.
 *
 *  @author Lukas Pospisil
 */


#ifndef PASC_COMMON_GLOBALMANAGER_H
#define	PASC_COMMON_GLOBALMANAGER_H

namespace pascinference {
namespace common {

/** \class GlobalManagerClass
 *  \brief Manipulation with processes informations - MPI environment.
 *
 *  Works also without MPI. In this case rank=0 and size=1.
*/
class GlobalManagerClass {
	private:
		int rank;			/**< the rank of this process */
		int size;			/**< number of processes */
		bool initialized;	/**< the rank and size were already obtained */

		/** @brief set rank and number of processes
		 * 
		 * If PetscVector is used, then MPI_Comm_rank and MPI_Comm_size called.
		 * Works also without PetscVector, in this case the problem is one-process, rank=0 and size=1.
		 * 
		 */
		void init(){
			if(!initialized){
			#ifdef USE_PETSC
				/* can be set after initialize of petsc */
				if(petscvector::PETSC_INITIALIZED){
					TRYCXX(PetscBarrier(NULL));
						
					MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
					MPI_Comm_size(MPI_COMM_WORLD, &this->size);

					initialized = true;
						
				}	
			#else
				initialized = true; /* if it is not with petsc, then this is always master */
				this->rank = 0;
				this->size = 1; /* one processor */
			#endif
			}
		}

	public:
		/** @brief default constructor
		 * 
		 * Calls init(), set rank and size.
		 * 
		 */
		GlobalManagerClass(){
			this->init();
		}

		/** @brief return rank of this process
		 * 
		 */
		int get_rank(){
			this->init();
			return this->rank;
		}
		
		/** @brief return number of processes
		 * 
		 */
		int get_size(){
			this->init();
			return this->size;
		}
};

static GlobalManagerClass GlobalManager; /**< for manipulation with rank and size of MPI */

}
} /* end of namespace */

#endif
