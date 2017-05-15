#include "external/seqarrayvector/algebra/vector/seqarrayvector.h"

namespace seqarrayvector {

seqarrayvector_all_type all;

SeqArrayVector::SeqArrayVector(){
	this->inner_size = -1;
}


SeqArrayVector::SeqArrayVector(int n){
	this->inner_size = n;
	this->inner_array = new double(this->inner_size);
}


SeqArrayVector::SeqArrayVector(double *values, int n){
	this->inner_size = n;
	this->inner_array = values; /* owner of this array is somewhere else */

}


SeqArrayVector::SeqArrayVector(const SeqArrayVector &vec){
	/* there is duplicate... this function has to be called as less as possible */
	this->inner_size = vec.size();
	this->inner_array = new double(this->inner_size);

	double *inner_array2 = vec.get_array();

	for(int i=0;i<this->inner_size;i++){
		this->inner_array[i] = inner_array2[i];
	}
	
}

SeqArrayVector::~SeqArrayVector(){
	/* if there is any inner vector, then destroy it */
	if(this->inner_size > 0){
		free(this->inner_array);
	}
}


void SeqArrayVector::set(double new_value){
	for(int i=0;i<this->inner_size;i++){
		this->inner_array[i] = new_value;
	}
}


void SeqArrayVector::set(int index, double new_value){
	if(index < inner_size){
		this->inner_array[index] = new_value;
	} else {
		//TODO: give error
	}
}

void SeqArrayVector::load_csv(std::string filename){
	//TODO
}

void SeqArrayVector::save_csv(std::string filename){
	//TODO
}

double* SeqArrayVector::get_array() const { 
	return this->inner_array;
}

int SeqArrayVector::size() const{
	return this->inner_size;
}


double SeqArrayVector::get(int i)
{
	return this->inner_array[i];
}


void SeqArrayVector::scale(double alpha){
	for(int i=0;i<this->inner_size;i++){
		this->inner_array[i] *= alpha;
	}
}


std::ostream &operator<<(std::ostream &output, const SeqArrayVector &vector)		
{
	double *inner_array = vector.get_array();
	
	output << "[";
	for(int i=0; i<vector.size(); i++){
		output << inner_array[i];
		if(i < vector.size()-1) output << ", ";
	}
	output << "]";
			
	return output;
}

/* vec1 = vec2, assignment operator (set vector) */
SeqArrayVector &SeqArrayVector::operator=(const SeqArrayVector &vec2){
	/* check for self-assignment by comparing the address of the implicit object and the parameter */
	/* vec1 = vec1 */
    if (this == &vec2){
        return *this;
	}

	/* vec1 is not initialized yet */
	if (this->inner_size <= 0){
		this->inner_size = vec2.size();
		this->inner_array = new double(this->inner_size);
	}

	//TODO: this->inner_size == vec2.size() ?

	/* copy values */
	double *inner_array2 = vec2.get_array();
	for(int i=0;i<this->inner_size;i++){
		this->inner_array[i] = inner_array2[i];
	}
	
	return *this;	
}

/* vec1 = scalar_value <=> vec1(all) = scalar_value, assignment operator */
SeqArrayVector &SeqArrayVector::operator=(double scalar_value){
	this->set(scalar_value);
	return *this;	
}

/* return subvector to be able to overload vector(index) = new_value */ 
SeqArrayVector SeqArrayVector::operator()(int index) const
{   
	double *new_array = new double(1);
	new_array[0] = this->inner_array[index];
	
	return SeqArrayVector(new_array,1);
}

/* return subvector vector(index_begin:index_end), i.e. components with indexes: [index_begin, index_begin+1, ..., index_end] */ 
SeqArrayVector SeqArrayVector::operator()(int index_begin, int index_end) const
{   
	double *new_array = new double(index_end-index_begin);

	for(int i=0;i<index_end-index_begin;i++){
		new_array[i] = this->inner_array[index_begin+i];
	}
	
	return SeqArrayVector(new_array, index_end-index_begin);
}

/* define SeqArrayVector(all) */
SeqArrayVector SeqArrayVector::operator()(seqarrayvector_all_type all_type) const{
	return *this;
} 


/* vec1 *= alpha */
void operator*=(SeqArrayVector &vec1, double alpha)
{
	vec1.scale(alpha);
}

/* vec1 += comb */
void operator+=(const SeqArrayVector &vec1, SeqArrayVector vec2)
{
	double *inner_array1 = vec1.get_array();
	double *inner_array2 = vec2.get_array();
	for(int i=0;i<vec1.inner_size;i++){
		inner_array1[i] += inner_array2[i]; 		
	}
}

/* vec1 -= comb */
void operator-=(const SeqArrayVector &vec1, SeqArrayVector vec2)
{
	double *inner_array1 = vec1.get_array();
	double *inner_array2 = vec2.get_array();
	for(int i=0;i<vec1.size();i++){
		inner_array1[i] -= inner_array2[i]; 		
	}
}

/* dot = dot(vec1,vec2) */
double dot(const SeqArrayVector &vec1, const SeqArrayVector &vec2)
{
	double dot_value=0;

	double *inner_array1 = vec1.get_array();
	double *inner_array2 = vec2.get_array();
	for(int i=0;i<vec1.size();i++){
		dot_value += inner_array1[i]*inner_array2[i];
	}

	return dot_value;
}

/* norm = norm_2(vec1) */
double norm(const SeqArrayVector &vec1)
{
	double dot_value = dot(vec1,vec1);
	return sqrt(dot_value);
}

/* max = max(vec1) */
double max(const SeqArrayVector &vec1)
{
	double max_value = -std::numeric_limits<double>::max(); /* - Inf */

	double *inner_array1 = vec1.get_array();
	for(int i=0;i<vec1.size();i++){
		if(inner_array1[i] > max_value){
			max_value = inner_array1[i];
		}
	}

	return max_value;
}

/* min = min(vec1) */
double min(const SeqArrayVector &vec1)
{
	double min_value = std::numeric_limits<double>::max(); /* + Inf */

	double *inner_array1 = vec1.get_array();
	for(int i=0;i<vec1.size();i++){
		if(inner_array1[i] < min_value){
			min_value = inner_array1[i];
		}
	}

	return min_value;
}

/* sum = sum(vec1) */
double sum(const SeqArrayVector &vec1)
{
	double sum_value = 0.0;
	
	double *inner_array1 = vec1.get_array();
	for(int i=0;i<vec1.size();i++){
		sum_value += inner_array1[i];
	}
	
	return sum_value;
}

/* vec3 = vec1./vec2 */
const SeqArrayVector operator/(const SeqArrayVector &vec1, const SeqArrayVector &vec2)
{
	//TODO
	return vec1;
}

SeqArrayVector mul(const SeqArrayVector &vec1, const SeqArrayVector &vec2)
{
	//TODO
	return vec1;

}



} /* end of namespace */

