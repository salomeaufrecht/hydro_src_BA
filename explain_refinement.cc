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


using namespace Dune;


#include "flowFunctions.hh"


// { grid_setup_begin }
int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    
    typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;

        // Start with a structured grid
    const std::array<unsigned, 2> n = {1, 1};
    const FieldVector<double, 2> lower = {0, 0};
    const FieldVector<double, 2> upper = {1, 1};

    
            
    std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    const GridView gridView01 = grid01->leafGridView();
    flowFragment f01 = {Dune::FieldVector<double, 2>{0, 0.5}, Dune::FieldVector<double, 2>{1, 0.5}, 0.3};
    std::vector<flowFragment> fragments01 = {f01};


    refineGridwithFragments(grid01, fragments01, 0.4, {1.0, 1.0});


    VTKWriter<GridView> vtkWriter(gridView01);
    vtkWriter.write("explain_refinement");
   

}