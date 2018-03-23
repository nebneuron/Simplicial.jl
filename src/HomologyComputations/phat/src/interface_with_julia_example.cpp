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
 
//This file contains a simple example that shows how to compte Betti numbers or
//persistent homology using phat via interface with Julia.
#include "interface_with_julia.cpp"

int main( int argc, char** argv ) 
{
//We will create a boundary matrix of an empty 2d cube. The numbers below
//indicate the ids of the cells. The filtration values are gieven in the 
//brackets.
//        (2)
// (0)1----5----2(0)
//    |         | 
//    |         | 
// (1)4         6(3)
//    |         |
//    |         |
// (0)0----7----3(0)
//        (4)

	std::vector< unsigned > dimensions_of_cells( 8 );
	dimensions_of_cells[0] = dimensions_of_cells[1] = dimensions_of_cells[2] = dimensions_of_cells[3] = 0;
	dimensions_of_cells[4] = dimensions_of_cells[5] = dimensions_of_cells[6] = dimensions_of_cells[7] = 1;


	std::vector< std::vector< unsigned > > sparse_boundary_matrix(8);
	 
	//sparse_boundary_matrix[0] - empty
	//sparse_boundary_matrix[1] - empty
	//sparse_boundary_matrix[2] - empty
	//sparse_boundary_matrix[3] - empty
	
	sparse_boundary_matrix[4].push_back(0);
	sparse_boundary_matrix[4].push_back(1);
	
	sparse_boundary_matrix[5].push_back(1);
	sparse_boundary_matrix[5].push_back(2);
	
	sparse_boundary_matrix[6].push_back(2);
	sparse_boundary_matrix[6].push_back(3);
	
	sparse_boundary_matrix[7].push_back(0);
	sparse_boundary_matrix[7].push_back(3);
			
	std::vector<  unsigned > betti = compute_Z2_Betti_numbers_from_sparse_matrix(sparse_boundary_matrix,dimensions_of_cells);
	
	std::cout << "Here are Betti numbers : \n";
	for ( size_t i = 0 ; i != betti.size() ; ++i )std::cout << betti[i] << " ";
	std::cout << std::endl;
	
	std::vector< std::pair< unsigned,unsigned > >  intervals = compute_persistence_from_sparse_matrix( sparse_boundary_matrix, dimensions_of_cells );
	
	std::cout << "And here are the persistence intervals : \n";
	for ( size_t i = 0 ; i != intervals.size() ; ++i )
	{
		std::cout << intervals[i].first << " " << intervals[i].second << std::endl;
	}
	std::cout << "That's all folks! \n";
												  
	return 0;
}
