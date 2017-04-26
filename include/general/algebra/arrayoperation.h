/** @file arrayoperation.h
 *  @brief class for manipulation with arrays
 *
 *  Defines some basic functions for manipulaton with arrays.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_ARRAYOPERATION_H
#define	PASC_COMMON_ARRAYOPERATION_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

namespace pascinference {
namespace algebra {

/** @brief print content of array to string
*
*  @param my_array pointer to array
*  @param my_size the size of array (number of elements)
*/
template<class MyType>
std::string print_array(MyType *my_array, int my_size){
	std::ostringstream out;
	out << "[";
	for(int i=0;i<my_size;i++){
		out << my_array[i];
		if(i<my_size-1) out << ",";
	}	
	out << "]";
	return out.str();
}

/** @brief print content of array with subarrays to string
*
*  @param my_array pointer to array
*  @param my_size1 the number of subarrays
*  @param my_size2 the size of subarray
*/
template<class MyType>
std::string print_array(MyType *my_array, int my_size1, int my_size2){
	std::ostringstream out;
	out << "{";
	for(int i=0;i<my_size1;i++){
		out << "[";
		for(int j=0;j<my_size2;j++){
			out << my_array[i*my_size2+j];
			if(j<my_size2-1) out << ",";
		}
		out << "]";
		if(i<my_size1-1) out << ",";
	}	
	out << "}";
	return out.str();
}


extern bool parse_strings_to_doubles(int K, int Km, std::vector<std::string> Theta_list, double *Theta_solution);


/** @brief print content of vector to string
*
*  @param my_vector vector to print
*/
template<class MyType>
std::string print_vector(std::vector<MyType> &my_vector){
	std::ostringstream out;
	out << "[";
	for(int i=0;i<my_vector.size();i++){
		out << my_vector[i];
		if(i<my_vector.size()-1) out << ",";
	}	
	out << "]";
	return out.str();
}


/** @brief print content of array
*
*  @param output where to print
*  @param my_size the size of array (number of elements)
*  @param my_array pointer to array
*/
template<class MyType>
void print_array(std::ostream &output, int my_size, MyType *my_array){
	output << "[";
	for(int i=0;i<my_size;i++){
		output << my_array[i];
		if(i<my_size-1) output << ",";
	}	
	output << "]";
}


/** @brief compute sum of elements in array
*
*  @param my_size the size of array (number of elements)
*  @param my_array pointer to array
*/
template<class MyType>
MyType sum_array(int my_size, const MyType *my_array){
	MyType sum = 0;
	for(int i=0;i<my_size;i++){
		sum += my_array[i];
	}	
	return sum;
}


/** @brief find largest elements in array
*
*  @param my_size the size of array (number of elements)
*  @param my_array pointer to array
*/
template<class MyType>
MyType max_array(int my_size, const MyType *my_array){
	MyType return_value = my_array[0];
	for(int i=1;i<my_size;i++){
		if(my_array[i] > return_value){
			return_value = my_array[i];
		}
	}	
	return return_value;
}


/** @brief find the index of largest elements in array
*
*  @param my_size the size of array (number of elements)
*  @param my_array pointer to array
*/
template<class MyType>
int max_arg_array(int my_size, const MyType *my_array){
	MyType max_value = my_array[0];
	int return_value;
	for(int i=1;i<my_size;i++){
		if(my_array[i] > return_value){
			max_value = my_array[i];
			return_value = i;
		}
	}	
	return return_value;
}


/** @brief set all elements in array to given value
*
*  @param my_size the size of array (number of elements)
*  @param my_array pointer to array
*  @param my_value new value of all elements
*/
template<class MyType>
void set_value_array(int my_size, MyType *my_array, MyType my_value){
	for(int i=0;i<my_size;i++){
		my_array[i] = my_value;
	}	
}


/** @brief find largest difference between cosequent elements in array
*
*  @param my_size the size of array (number of elements) > 1
*  @param my_array pointer to array
*/
template<class MyType>
MyType max_diff_array(int my_size, const MyType *my_array){
	MyType return_value = my_array[1] - my_array[0];
	for(int i=2;i<my_size;i++){
		if((my_array[i]-my_array[i-1]) > return_value){
			return_value = my_array[i]-my_array[i-1];
		}
	}	
	return return_value;
}


/** @brief compute sum of some elements in array
*
*  @param start index of starting element of subarray
*  @param end index of ending element of subarray
*  @param my_array pointer to array
*/
template<class MyType>
MyType sum_subarray(int start, int end, const MyType *my_array){
	MyType sum = 0;
	for(int i=start;i<=end;i++){
		sum += my_array[i];
	}	
	return sum;
}


/** @brief compute dot product of two arrays
*
*  @param size the size of arrays (number of elements)
*  @param my_array1 pointer to array
*  @param my_array2 pointer to array
*/
template<class MyType>
MyType dot_arrays(int size, const MyType *my_array1, const MyType *my_array2){
	MyType sum = 0;
	for(int i=0;i<=size;i++){
		sum += my_array1[i]*my_array2[i];
	}	
	return sum;
}


/** @brief compare two object and return larger one
*
*  @param a1 object
*  @param a2 object
*/
template<class MyType>
MyType max(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 > a1){
		return_value = a2;
	}
	return return_value;
}


/** @brief compare two object and return smaller one
*
*  @param a1 object
*  @param a2 object
*/
template<class MyType>
MyType min(const MyType a1, const MyType a2){
	MyType return_value = a1;
	if(a2 < a1){
		return_value = a2;
	}
	return return_value;
}


/** @brief compute point-wise multiplication of two arrays
*
*  my_array[i] = my_array1[i]*my_array2[i]
* 
*  @param my_size the size of arrays (number of elements)
*  @param my_array pointer to output array
*  @param my_array1 pointer to input array
*  @param my_array2 pointer to input array
*/
template<class MyType>
void mult_pw_array(int my_size, MyType *my_array, const MyType *my_array1, const MyType *my_array2){
	for(int i=0;i<my_size;i++){
		my_array[i] = my_array1[i]*my_array2[i];
	}	
}

void myround(double in, double *out);
std::string printbool(bool input);
void arg_parse(const char *args, int *argc, char ***argv);

std::vector<std::string> split(const std::string &s, char delim);

template<typename Out> void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

}
} /* end of namespace */

#endif
