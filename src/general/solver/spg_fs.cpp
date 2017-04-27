#include "general/solver/spg_fs.h"

namespace pascinference {
namespace solver {
	
/* constructor */
SPG_fs::SPG_fs(int new_m){
	this->m = new_m;
	this->fs_list = (double*)malloc(this->m*sizeof(double));
}

SPG_fs::~SPG_fs(){
	free(this->fs_list);
}

/* init the list with function values using one initial fx */
void SPG_fs::init(double fx){
	LOG_FUNC_BEGIN

	for(int i=0; i<this->m;i++){
		this->fs_list[i] = fx;
	}
	this->last_idx = 0;

	LOG_FUNC_END
}

/* get the size of the list */
int SPG_fs::get_size(){
	return this->m;
}

/* get the value of max value in the list */
double SPG_fs::get_max(){
	LOG_FUNC_BEGIN
	
	int max_idx = 0;
	double max_value = this->fs_list[max_idx];
	for(int i=0;i<this->m;i++){
		if(this->fs_list[i] > max_value){
			max_idx = i;
			max_value = this->fs_list[max_idx];
		}
	}

	LOG_FUNC_END

	return max_value;
}

/* update the list by new value - pop the first and push the new value (FIFO) */
void SPG_fs::update(double new_fx){
	LOG_FUNC_BEGIN

	this->last_idx++;
	if(this->last_idx >= this->m){
		this->last_idx = 0;
	}
	
	this->fs_list[this->last_idx] = new_fx;

	LOG_FUNC_END
}

/* print the content of the list */
void SPG_fs::print(ConsoleOutput &output)
{
	output << "[ ";
	/* for each component go throught the list */
	for(int i=this->last_idx;i<this->last_idx+this->m;i++){
		if(i < this->m){
			output << this->fs_list[i];
		} else {
			output << this->fs_list[i - this->m];
		}
		
		if(i < this->last_idx+this->m-1){ /* this is not the last node */
				output << ", ";
		}
	}
	output << " ]";
}

}
} /* end namespace */

