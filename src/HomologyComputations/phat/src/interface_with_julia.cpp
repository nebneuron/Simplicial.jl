	/*
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2018 Swansea University UK
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 
#include "../include/phat/compute_persistence_pairs.h"

// main data structure (choice affects performance)
#include "../include/phat/representations/vector_vector.h"

// algorithm (choice affects performance)
#include "../include/phat/algorithms/standard_reduction.h"
#include "../include/phat/algorithms/chunk_reduction.h"
#include "../include/phat/algorithms/row_reduction.h"
#include "../include/phat/algorithms/twist_reduction.h"

#include <cassert>

/**
 * This is an auxiliary procedure. For datailed descripion of input parameters
 * please consult the following procedures:
 * compute_Z2_Betti_numbers_from_sparse_matrix and
 * compute_persistence_from_sparse_matrix 
**/ 
template <typename PHAT_reduction_algorithm = phat::twist_reduction , typename PHAT_collumn_representation = phat::bit_tree_pivot_column >
phat::persistence_pairs reduce_boundary_matrix_with_phat
 ( const std::vector< std::vector< unsigned > >& sparse_boundary_matrix, 
   const std::vector< unsigned >& dimensions_of_cells )
 {
	 assert( sparse_boundary_matrix.size() == dimensions_of_cells.size() );
	 
	//first define a boundary matrix with the chosen internal representation
    phat::boundary_matrix<PHAT_collumn_representation > boundary_matrix;

    //set the number of columns (equal to the size of compute_persistence_from_sparse_matrix)
    boundary_matrix.set_num_cols( sparse_boundary_matrix.size() );
    
    //set the dimension of the cell that a column represents:
    for ( size_t i = 0 ; i != dimensions_of_cells.size() ; ++i )
    {
		boundary_matrix.set_dim( i, dimensions_of_cells[i] );
	}    

    //input the boundary matrix into phat format.
    std::vector< phat::index > temp_col;    
    for ( size_t i = 0 ; i != sparse_boundary_matrix.size() ; ++i )
    {
		for ( size_t j = 0 ; j != sparse_boundary_matrix[i].size() ; ++j )
		{
			temp_col.push_back( sparse_boundary_matrix[i][j] );
		}		
		std::sort( temp_col.begin() , temp_col.end() );
		boundary_matrix.set_col( i, temp_col );
		temp_col.clear();
	}

    // define the object to hold the resulting persistence pairs
    phat::persistence_pairs pairs;

    //reduce the matrix.
    phat::compute_persistence_pairs< PHAT_reduction_algorithm >( pairs, boundary_matrix );
    
    //sort the persistence pairs by birth index.
    pairs.sort();
    
    return pairs;    
}

 
 /**
  * This is the compute_Z2_Betti_numbers_from_sparse_matrix procedure. The input parameter
  * is a std::vector< std::vector< unsigned > >  sparse boundary matrix. 
  * It is assumed that every cell in the considered complex is enumerated and the enumeration 
  * is compatible with filtration. If there are no filtration we assume that the bounday elements
  * of a cell having index i all have indices < i. 
  * 
  * The second parameter is a vector of nonnegative integeres encoding the dimensions of the following		
  * cells.
  * 
  * The value returned by the procedure is a vector of Betti numbers.
 **/ 
 std::vector< unsigned >  compute_Z2_Betti_numbers_from_sparse_matrix
												( const std::vector< std::vector< unsigned > >& sparse_boundary_matrix, 
												  const std::vector< unsigned >& dimensions_of_cells )
{
	bool dbg = false;
	if ( dbg )std::cerr << "Using compute_Z2_Betti_numbers_from_sparse_matrix procedure \n";
	
	phat::persistence_pairs pairs = reduce_boundary_matrix_with_phat(sparse_boundary_matrix, dimensions_of_cells);

	if ( dbg )std::cerr << "Done reducing boundary matrix. \n";
	
	
	//check the dimension of the complex:
	unsigned dimension_of_the_complex = 0;
	for ( size_t i = 0 ; i != dimensions_of_cells.size() ; ++i )
	{
		if ( dimension_of_the_complex < dimensions_of_cells[i] )dimension_of_the_complex = dimensions_of_cells[i];
	}	
	
	//prepare the output data.
	std::vector< unsigned > result( dimension_of_the_complex+1,0 );
	
	//We need to check which cells are not paired. They will constitute the Betti numbers.
	std::vector< bool > which_cells_were_not_reduced( sparse_boundary_matrix.size() , false );	
	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
	{
        which_cells_were_not_reduced[ pairs.get_pair( idx ).first ] = true;
        which_cells_were_not_reduced[ pairs.get_pair( idx ).second ] = true;
	}
	
	
	//and now we compute the Betti numbers.
	for ( size_t i = 0 ; i != which_cells_were_not_reduced.size() ; ++i )
	{		
		if ( which_cells_were_not_reduced[i] )continue;
		++result[ dimensions_of_cells[i] ];
	}
		
	return result;    
}//compute_Z2_Betti_numbers_from_sparse_matrix


 /**
  * This is the compute_Z2_Betti_numbers_from_sparse_matrix procedure. The input parameter
  * is a std::vector< std::vector< unsigned > >  sparse boundary matrix. 
  * It is assumed that every cell in the considered complex is enumerated and the enumeration 
  * is compatible with filtration. If there are no filtration we assume that the bounday elements
  * of a cell having index i all have indices < i. 
  * 
  * The second parameter is a vector of nonnegative integeres encoding the dimensions of the following
  * cells.
  * 
  * The value returned by the procedure is a vector of all pair of indices. To retrive persistence
  * pairs out of it, one need to change a pair <index_1,index_2> into <filtration of index_1, filtration of index_2>
 **/  
std::vector< std::pair< unsigned,unsigned > >  compute_persistence_from_sparse_matrix
												( const std::vector< std::vector< unsigned > >& sparse_boundary_matrix, 
												  const std::vector< unsigned >& dimensions_of_cells )
{
	phat::persistence_pairs pairs = reduce_boundary_matrix_with_phat(sparse_boundary_matrix, dimensions_of_cells);

	//prepare the result.
	std::vector< std::pair< unsigned, unsigned > > result;
	result.reserve( sparse_boundary_matrix.size()/2 );
	
	//check the maximal dimension involved:
	std::vector< bool > which_cells_were_not_reduced( sparse_boundary_matrix.size() , false );	
	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
	{
		result.push_back( std::pair<unsigned,unsigned>(pairs.get_pair( idx ).first , pairs.get_pair( idx ).second) );
		which_cells_were_not_reduced[ pairs.get_pair( idx ).first ] = true;
		which_cells_were_not_reduced[ pairs.get_pair( idx ).second ] = true;
	}
	
	for ( size_t i = 0 ; i != which_cells_were_not_reduced.size() ; ++i )
	{
		if ( which_cells_were_not_reduced[i] )continue;//in this case, this cell is an element of a pair.
		result.push_back( std::pair< unsigned,unsigned >( i , std::numeric_limits<unsigned>::max() ) );
	}
	
	return result;    
}//compute_persistence_from_sparse_matrix


/**
 * For some reason Julia do not like vector< pair >, so I will use this:
**/ 
std::vector< unsigned >  compute_persistence_from_sparse_matrix_return_vector
												( const std::vector< std::vector< unsigned > >& sparse_boundary_matrix, 
												  const std::vector< unsigned >& dimensions_of_cells )
{
	std::vector< std::pair< unsigned,unsigned > > intervals =  compute_persistence_from_sparse_matrix
												(sparse_boundary_matrix, dimensions_of_cells );
	std::vector< unsigned > result;
	result.reserve( 2*intervals.size() );
	for ( size_t i = 0 ; i != intervals.size() ; ++i )										
	{
		result.push_back( intervals[i].first );
		result.push_back( intervals[i].second );
	}
	return result;    
}//compute_persistence_from_sparse_matrix_return_vector




template <class RandomIt, class OutputIt>
OutputIt compute_betti_numbers_not_optimal(const std::uint64_t number_of_cells,
					 const std::uint64_t dimension, 
					 RandomIt xbegin, RandomIt xend,
                     OutputIt rbegin) {   
  
  
  std::vector< std::vector< unsigned > > sparse_boundary_matrix;
  sparse_boundary_matrix.reserve(number_of_cells);
  std::vector< unsigned > dimensions_of_cells;
  dimensions_of_cells.reserve(number_of_cells);
  
  std::size_t position = 0;				  
  for (std::size_t i = 0; i < number_of_cells; i++) 
  {
	  //each cell is encoded in a vector by the following string:
	  //cell id (to check if it is equal i+1)
	  //dimension
	  //and then the sequence of boundary elements
	  //-1 indicate the end of a cell.
	  
	  //first we check if the cell id  is the right one:
	  size_t cell_id = (size_t)xbegin[position];
	  assert ( i+1 == cell_id );
	  ++position;
	  
	  //we read the dimension of a cell and store it:
	  unsigned dimenson = (unsigned)xbegin[position];
	  dimensions_of_cells.push_back( dimenson );
	  ++position;
	  
	  //and now we read the vector of boundary
	  std::vector< unsigned > boundary;
	  boundary.reserve(20);
	  while ( xbegin[position] != -1 )
	  {
		  //-1, since we convert from Julia to C++ style.
		  boundary.push_back( (unsigned)xbegin[position]-1 );
		  ++position;
	  }
	  ++position;
	  sparse_boundary_matrix.push_back( boundary );
		
  }
    
  std::vector< unsigned > result = compute_Z2_Betti_numbers_from_sparse_matrix
								   ( sparse_boundary_matrix, dimensions_of_cells );
  
  for ( size_t i = 0 ; i != result.size() ; ++i )
  {
	   *rbegin++ = result[i];
  }
  return rbegin;
}//compute_betti_numbers_not_optimal
















/**
 * This is an auxiliary procedure to create Phat boundary matrix and to 
 * reduce it. It is assumed to take as an input the encoded vector representing 
 * boundary matrix. Detailed description of the content of the vector can 
 * be found in the procedure body below. 
**/ 
template 
<class RandomIt, 
class OutputIt, 
typename PHAT_reduction_algorithm = phat::twist_reduction, 
typename PHAT_collumn_representation = phat::bit_tree_pivot_column >
std::pair<phat::persistence_pairs , std::vector<unsigned> >
reduce_boundary_matrix_with_phat_C_style_input
(const std::uint64_t number_of_cells,
const std::uint64_t dimension, 
RandomIt xbegin, RandomIt xend,
OutputIt rbegin )
 {
	 bool dbg = false;
	 
	//first define a boundary matrix with the chosen internal representation
	phat::boundary_matrix<PHAT_collumn_representation > boundary_matrix;

	//set the number of columns (equal to the size of compute_persistence_from_sparse_matrix)
	boundary_matrix.set_num_cols( number_of_cells );	 
	std::vector<unsigned> dimensions;
	dimensions.reserve( number_of_cells );
	 
	std::size_t position = 0;				  
	std::vector< phat::index > temp_col;    
	for (std::size_t i = 0; i < number_of_cells; i++) 
	{
		//each cell is encoded in a vector by the following string:
		//cell id (to check if it is equal i+1)
		//dimension
		//and then the sequence of boundary elements
		//-1 indicate the end of a cell.

		//first we check if the cell id  is the right one:
		size_t cell_id = (size_t)xbegin[position];
		assert ( i+1 == cell_id );
		++position;

		//we read the dimension of a cell and store it:
		unsigned dimension_ = (unsigned)xbegin[position];	 
		++position;

		dimensions.push_back( dimension_ );
		boundary_matrix.set_dim( i, dimension_ );
		
		if ( dbg )
		{
			std::cout << "This is a cell havng id : " << cell_id << " and dimension : " << dimension_ << " and the following boundary elements : " << std::endl;
		}

		//and now we read the vector of boundary
		std::vector< unsigned > boundary;
		boundary.reserve(20);
		while ( xbegin[position] != -1 )
		{
		  //-1, since we convert from Julia to C++ style.		  
		  temp_col.push_back( (unsigned)xbegin[position]-1 );
		  if ( dbg )
		  {
			  std::cout << (unsigned)xbegin[position] << std::endl;
		  }
		  ++position;
		}	  
		boundary_matrix.set_col( i, temp_col );
		temp_col.clear();

		++position;	
	}
 

    // define the object to hold the resulting persistence pairs
    phat::persistence_pairs pairs;

    //reduce the matrix.
    phat::compute_persistence_pairs< PHAT_reduction_algorithm >( pairs, boundary_matrix );
    
    //sort the persistence pairs by birth index.
    pairs.sort();
    
    return std::make_pair(pairs,dimensions);    
}
























template <class RandomIt, class OutputIt>
OutputIt compute_betti_numbers(const std::uint64_t number_of_cells,
					 const std::uint64_t dimension, 
					 RandomIt xbegin, RandomIt xend,
                     OutputIt rbegin) {   
     bool dbg = false;
     
     if ( dbg )
     {
		std::cerr << "We are using compute_betti_numbers function. \n";
	 }
     
	std::pair<phat::persistence_pairs,std::vector<unsigned> > output_ = 
	reduce_boundary_matrix_with_phat_C_style_input
	(number_of_cells, dimension, xbegin, xend, rbegin);


	phat::persistence_pairs pairs = output_.first;
	std::vector<unsigned> dimensions = output_.second;


	//prepare the output data.
	std::vector< unsigned > result( dimension+1,0 );
	
	if ( dbg )std::cout << "Here is the number of cells : " << number_of_cells << std::endl;
	
	//We need to check which cells are not paired. They will constitute the Betti numbers.
	std::vector< bool > which_cells_were_not_reduced( number_of_cells , false );	
	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
	{
		if ( dbg )
		{
			std::cout << "We have the following pair : " << pairs.get_pair( idx ).first << " , " << pairs.get_pair( idx ).second << std::endl;
		}		
        which_cells_were_not_reduced[ pairs.get_pair( idx ).first ] = true;
        which_cells_were_not_reduced[ pairs.get_pair( idx ).second ] = true;
	}
	
	//and now we compute the Betti numbers.
	for ( size_t i = 0 ; i != which_cells_were_not_reduced.size() ; ++i )
	{
		if ( which_cells_were_not_reduced[i] )continue;
		if ( dbg )std::cout << "The cell number : " << i << " was not reduced. Its dimension is : " << dimensions[i] << std::endl;
		++result[ dimensions[i] ];
	}
   
	for ( size_t i = 0 ; i != result.size() ; ++i )
	{
	   *rbegin++ = result[i];
	   //std::cerr << "result[i] : " << result[i] << std::endl;
	}
	return rbegin;
}//compute_betti_numbers






template <class RandomIt, class OutputIt>
OutputIt compute_persistence_intervals(const std::uint64_t number_of_cells,
const std::uint64_t dimension,
RandomIt xbegin, RandomIt xend,
OutputIt rbegin) 
{  
    bool dbg = false;
     
    if ( dbg )
    {
		std::cerr << "We are using compute_persistence_intervals function. \n";
	}
     
	std::pair<phat::persistence_pairs,std::vector<unsigned> > output_ =
	reduce_boundary_matrix_with_phat_C_style_input
	(number_of_cells, dimension, xbegin, xend, rbegin);

	phat::persistence_pairs pairs = output_.first;

	//check the maximal dimension involved:
	std::vector< bool > which_cells_were_not_reduced( number_of_cells , false );
	for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
	{
		which_cells_were_not_reduced[ pairs.get_pair( idx ).first ] = true;
		which_cells_were_not_reduced[ pairs.get_pair( idx ).second ] = true;
		*rbegin++ =  pairs.get_pair( idx ).first;
		*rbegin++ = pairs.get_pair( idx ).second;
		if ( dbg )
		{
			std::cout << "( " << pairs.get_pair( idx ).first << " , " << pairs.get_pair( idx ).second << " ) ";
		}
	}
	*rbegin++ = -1;
	for ( size_t i = 0 ; i != which_cells_were_not_reduced.size() ; ++i )
	{
		if ( which_cells_were_not_reduced[i] )continue;//in this case, this cell is an element of a pair.
		*rbegin++ = i;		
		if ( dbg )
		{
			std::cout << "Infinite interval starting at  " << i << " ";
		}
	}
	return rbegin;
}//compute_persistence_intervals











extern "C" 
{
	void compute_betti_numbers_not_optimal(const std::uint64_t number_of_cells /* use fixed-size integers */,
					 const std::uint64_t dimension /* use fixed-size integers */,
					 const std::uint64_t *xbegin /* use fixed-size integers */,
					 std::int64_t *rbegin /* use fixed-size integers */) 
	{
		//Here are various functions that allow to carry on the computations:
		compute_betti_numbers_not_optimal
		(number_of_cells, 
		dimension,
		xbegin, 
		xbegin + number_of_cells, 
		rbegin);
	}
	
	
	void compute_betti_numbers(const std::uint64_t number_of_cells /* use fixed-size integers */,
					 const std::uint64_t dimension /* use fixed-size integers */,
					 const std::uint64_t *xbegin /* use fixed-size integers */,
					 std::int64_t *rbegin /* use fixed-size integers */) 
	{
		//Here are various functions that allow to carry on the computations:
		compute_betti_numbers
		(number_of_cells, 
		dimension,
		xbegin, 
		xbegin + number_of_cells, 
		rbegin);
	}
	
	
	
	void compute_persistence_intervals(const std::uint64_t number_of_cells /* use fixed-size integers */,
	const std::uint64_t dimension /* use fixed-size integers */,
	const std::uint64_t *xbegin /* use fixed-size integers */,
	std::int64_t *rbegin /* use fixed-size integers */) 
	{
		//Here are various functions that allow to carry on the computations:
		compute_persistence_intervals
		(number_of_cells, 
		dimension,
		xbegin, 
		xbegin + number_of_cells, 
		rbegin);
	}
}











//Here is an example of a complex with a filtration:
//
//      1---(8)---2
//      |         |
//     (5)       (7)
//      |         |
//      4---(6)---3
//
//And here is the array to encode it:
//1 0 -1 2 0 -1 3 0 -1 4 0 -1 5 1 1 4 -1 6 1 4 3 -1 7 1 2 3 -1 8 1 1 2 -1
//
//Here is a bit more complicated stuff:
//      1---(8)---2---(12)---6
//      |         |          |
//     (7)       (10)       (13)
//      |         |          |
//      4---(9)---3---(11)---5
//And here is the array to encode it:
//1 0 -1 2 0 -1 3 0 -1 4 0 -1 5 0 -1 6 0 -1 7 1 1 4 -1 8 1 1 2 -1 9 1 4 3 -1 10 1 2 3 -1 11 1 3 5 -1 12 1 2 6 -1 13 1 5 6 -1




/**
THIS IS HOW TO RUN IT VIA INTERFACE THAT USE PURE C

Compile this C++ file to a shared object:
g++ -std=c++11  -shared -O3 -fPIC -o interface_with_julia.so interface_with_julia.cpp

Then open Julia and run: 

const libhandle = Libdl.dlopen(joinpath(pwd(), "interface_with_julia.so"))

//here are the options for functions to use to compute Betti numbers or persistence:
const funhandle = Libdl.dlsym(libhandle, :compute_betti_numbers_not_optimal) //non optymality is here because the boundary matrix is stored twice (once as vector, and the other time as phat boundary matrix).
const funhandle = Libdl.dlsym(libhandle, :compute_betti_numbers)
const funhandle = Libdl.dlsym(libhandle, :compute_persistence_intervals)  //the output of this function is a vector of creating cell and killing cell, creatinc cell and killing cell etc etc.
																		  //The infinite pairs comes at the end, they are separated from the other pairs by -1. In this case, each id of a cell that comes after -1 is a creator of infinite homology class
																		  //REMEMBER that when computing persistence intervals, the size of the array have to be the number of cells + 1 (to store the -1 to separate finite from infinite pairs!!)
																		  


function compute(number_of_cells::UInt64,dimension::UInt64,A::Array{Int64,2})
    result = Vector{Int64}(dimension)
    ccall(funhandle, Void, (UInt64, UInt64, Ref{Int64}, Ref{Int64}),
          number_of_cells, dimension, A, result)
    return result
end





cells = [1 0 -1 2 0 -1 3 0 -1 4 0 -1 5 1 1 4 -1 6 1 4 3 -1 7 1 2 3 -1 8 1 1 2 -1]

result = compute(convert(UInt64,8),convert(UInt64,2),cells)
* 
* 
THIS IS HOW TO DO IT VIA INTERFACE THAT USES C++ VIA Cxx:
 
Compile this C++ file to a shared object:
g++ -std=c++11  -shared -O3 -fPIC -o interface_with_julia.so interface_with_julia.cpp

Then open Julia and run: 

using Cxx
const path_to_lib = pwd()
addHeaderDir(path_to_lib, kind=C_System)
Libdl.dlopen(joinpath(path_to_lib, "interface_with_julia"), Libdl.RTLD_GLOBAL)
cxxinclude("interface_with_julia.cpp")


cells = [ [],[],[],[],[0,3],[2,3],[1,2],[0,1] ]       // NOTE THAT ENUMERATION STAT FRoM ZERO HERE
cpp_cells = convert(cxxt"std::vector< std::vector< unsigned > >", cells)
dimensions = [0,0,0,0,1,1,1,1]
cpp_dimensions = convert(cxxt"std::vector<unsigned>", dimensions)


To compute Betti numbers use:
cxx_Betti_numbers = icxx"compute_Z2_Betti_numbers_from_sparse_matrix($cpp_cells,$cpp_dimensions);"
println("Cxx sums: $(collect(cxx_Betti_numbers))") 

To compute persistence intervals use:
cxx_persistence = icxx"compute_persistence_from_sparse_matrix_return_vector($cpp_cells,$cpp_dimensions);"
println("Cxx sums: $(collect(cxx_persistence))")   //THE OUTPUT FORMAT IS A COLLECTION OF PAIRS PACKED IN A VECTOR (JULIA DO NOT LIKE VECTOR<PAIR>. THEY HAVE TO BE UNPACKED IN CPP LEVEL
 
 



For the further reference, here are my questions on forums that helped me to get this done:
https://discourse.julialang.org/t/calling-c-function-having-std-vectors-as-input-and-output-parameters-from-julia/9224
https://stackoverflow.com/questions/48851198/calling-c-function-having-stdvectors-as-input-and-output-parameters-from-jul/48956351#48956351
**/
