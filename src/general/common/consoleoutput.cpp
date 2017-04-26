#include "general/common/consoleoutput.h"

namespace pascinference {
namespace common {

OffsetClass::OffsetClass() {
	size = 0;
}
		
void OffsetClass::push() {
	size += 3;
}

void OffsetClass::pop(){
	size -= 3;
	if(size < 0) size = 0;
}

std::ostream &operator<<(std::ostream &output, const OffsetClass &my_offset){
	int i;
	for(i=0;i<my_offset.size;i++){
		output << " ";
	}
	return output;
}

OffsetClass offset;


#ifdef USE_PETSC

ConsoleOutput::ConsoleOutputBuf::ConsoleOutputBuf(std::ostream& str):output(str){
}

ConsoleOutput::ConsoleOutputBuf::~ConsoleOutputBuf(){
}

int ConsoleOutput::ConsoleOutputBuf::sync(){
	if(this->rank == print_rank || this->print_rank == -1){
		std::stringstream output_string;

		/* write here also a rank of processor */
		if(GlobalManager.get_size() > 1){
			output_string << "[" << GlobalManager.get_rank() << "] " << offset << str();
		} else {
			output_string << offset << str();
		}
		str("");

//		output_string << output.rdbuf();
		if(this->print_rank == 0 || GlobalManager.get_size() <= 1){
			/* master prints */
			TRYCXX( PetscPrintf(PETSC_COMM_WORLD, output_string.str().c_str()) );
		} else {
			/* all prints */
			TRYCXX( PetscSynchronizedPrintf(PETSC_COMM_WORLD, output_string.str().c_str()) );
		}

		output_string.str("");
		output_string.clear();

	} else {
	    str("");
	}

//	output.flush();
	return 0;
}

ConsoleOutput::ConsoleOutput(std::ostream& str, int rank) : std::ostream(&buffer), buffer(str){
	buffer.print_rank = rank;
	buffer.rank = GlobalManager.get_rank();
}

void ConsoleOutput::push(){
	offset.push();
}

void ConsoleOutput::pop(){
	offset.pop();
}

void ConsoleOutput::synchronize(){
	if(buffer.print_rank == -1){
		TRYCXX( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
	}
}

#else 

ConsoleOutput::ConsoleOutputBuf::ConsoleOutputBuf(std::ostream& str):output(str){

}

ConsoleOutput::ConsoleOutputBuf::~ConsoleOutputBuf(){
}

void ConsoleOutput::ConsoleOutputBuf::set_rank(){
	if(!rankset){
		#ifdef USE_PETSC
			/* can be set after initialize of petsc */
			if(petscvector::PETSC_INITIALIZED){
				TRYCXX(PetscBarrier(NULL));
						
				MPI_Comm_rank(MPI_COMM_WORLD, &this->rank);
				rankset = true;
				
			}	
		#else
			rankset = true; /* if it is not with petsc, then this is always master */
		#endif
	}
}
				
int ConsoleOutput::ConsoleOutputBuf::sync(){
	set_rank();
	if(this->rank == print_rank || this->print_rank == -1){
		#ifdef USE_PETSC
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

	#ifdef USE_PETSC
		TRYCXX(PetscBarrier(NULL)); /* barrier after each cout */
	#endif

	return 0;
}
				
ConsoleOutput::ConsoleOutput(std::ostream& str, int rank) : std::ostream(&buffer), buffer(str) {
	buffer.rankset = false;
	buffer.set_rank();
	buffer.print_rank = rank;
}

void ConsoleOutput::push(){
	offset.push();
}

void ConsoleOutput::pop(){
	offset.pop();
}

void ConsoleOutput::synchronize(){
}

#endif

}
} /* end of namespace */

