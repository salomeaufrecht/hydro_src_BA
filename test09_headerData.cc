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

using namespace Dune;



int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;

    

    // Start with a structured grid
    const std::array<unsigned, dim> n = {6, 6};
    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {6, 6};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    const GridView gridView = grid->leafGridView();



    

    flowFragment f1 = {{3.5, 3.5}, {6, 3.5}, 9};
    flowFragment f2 = {{0, 0}, {3.5, 3.5}, 0.1, 0.1};
    flowFragment f3 = {{1.5, 1.5}, {1.5, 2}, 0.1};
    flowFragment f4 = {{1.5, 2}, {1.5, 3}, 0.2};
    flowFragment f5 = {{1.5, 3}, {1.5, 4}, 0.3};
    flowFragment f6 = {{1.5, 4}, {1.5, 5}, 0.4};
    flowFragment f7 = {{1.5, 5}, {1.5, 6}, 0.5};

    flowWithFragments(grid, f1);
    
    std::vector<flowFragment> fragments = {f1, f2, f3, f4, f5, f6, f7};

    for(auto f : fragments){
        flowWithFragments(grid, f);
    }

    std::vector<double> height(gridView.indexSet().size(2));
    for(auto f : fragments){
        height = applyFlowHeightFragments(grid, f, height);
    }


    

    // Write grid to file
    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.addVertexData(height, "height");
    vtkWriter.write("refined_grid_with_height");
    
}
