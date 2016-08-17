/** @file arrayoperation.h
 *  @brief class for manipulation with arrays
 *
 *  Defines some basic functions for manipulaton with arrays.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_ARRAYOPERATION_H
#define	PASC_COMMON_ARRAYOPERATION_H

namespace pascinference {

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


} /* end of namespace */

#endif
