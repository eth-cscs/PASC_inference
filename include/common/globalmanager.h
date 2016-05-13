#ifndef PASC_COMMON_GLOBALMANAGER_H
#define	PASC_COMMON_GLOBALMANAGER_H

namespace pascinference {
	
class GlobalManagerClass {
	private:
		int rank;
		int size;
		bool initialized; /**< the rank and size were already obtained */

		void init(){
			if(!initialized){
			#ifdef USE_PETSCVECTOR
				/* can be set after initialize of petsc */
				if(petscvector::PETSC_INITIALIZED){
					TRY(PetscBarrier(NULL));
						
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
		GlobalManagerClass(){
			this->init();
		}

		int get_rank(){
			this->init();
			return this->rank;
		}
		
		int get_size(){
			this->init();
			return this->size;
		}
};

static GlobalManagerClass GlobalManager; /**< for manipulation with rank and size of MPI */

} /* end of namespace */

#endif
