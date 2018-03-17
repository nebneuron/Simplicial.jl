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


#include <cstdint>
#include <iterator>


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
	phat::persistence_pairs pairs = reduce_boundary_matrix_with_phat(sparse_boundary_matrix, dimensions_of_cells);

	
	
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
 * For simple test if this all make sense.
**/ 
void simple_test()
{
	std::cout << "Working \n";
}

void pring_number( int i )
{
	std::cout << "Here is the number : " << i << "\n";
}

void print_vector( const std::vector<int>& vect )
{
	std::cout << "Here is the vector \n";
	for ( int i = 0 ; i != vect.size() ; ++i )
	{
		std::cout << vect[i] << " ";
	}
	std::cout << std::endl;
}


template <class RandomIt, class OutputIt>
OutputIt compute_sum(const std::uint64_t nrows, RandomIt xbegin, RandomIt xend,
                     OutputIt rbegin) {
  const std::size_t ncols{std::distance(xbegin, xend) / nrows};
  typename std::iterator_traits<OutputIt>::value_type sum{0};
  
  for (std::size_t row = 0; row < nrows; row++) {
    for (std::size_t col = 0; col < ncols; col++)
      sum += xbegin[col * nrows + row];

    *rbegin++ = sum;
    sum = 0;
  }
  return rbegin;
}



extern "C" {
void compute_sum(const std::uint64_t m /* use fixed-size integers */,
                 const std::uint64_t n /* use fixed-size integers */,
                 const std::int64_t *xbegin /* use fixed-size integers */,
                 std::int64_t *rbegin /* use fixed-size integers */) {
  compute_sum(m, xbegin, xbegin + m * n, rbegin);
}
}
