#include "config.h"
#include <iostream>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#include <iostream>
#include <array>
#include "flowFunctions.hh"

/*
for every triangle checks if corners are above, below or inbetween the edges --> if corners are on different sides border within trinagle

*/

//using namespace Dune;



int main(int argc, char *argv[])
{
    // Set up MPI if available
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;

    

    // Start with a structured grid
    const std::array<unsigned, dim> n = {1, 1};
    const Dune::FieldVector<double, dim> lower = {0, 0};
    const Dune::FieldVector<double, dim> upper = {1, 1};

    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);
    const GridView gridView = grid->leafGridView();

    flowFragment f1 = {{0.5, 0}, {0.5, 1}, 0.5};
    //flowFragment f1 = {{0.5, 0}, {0.5, 1}, 0.5};
    flowWithFragments(grid, f1);

    std::vector<double> height(gridView.indexSet().size(2));
    height = applyFlowHeightFragments(grid, f1, height);



    //for(auto& v : vertices(gridView)){
    //    height[gridView.indexSet().index(v)] = v.geometry().corner(0)[0];
    //}
   

    // Write grid to file
    Dune::VTKWriter<GridView> vtkWriter(gridView);
    

    vtkWriter.addVertexData(height, "test");
    vtkWriter.write("refined_grid_single");
    
}
