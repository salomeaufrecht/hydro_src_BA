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
    const std::array<unsigned, dim> n = {1, 1};
    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {1, 1};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);
    const GridView gridView = grid->leafGridView();

    std::vector<double> w;
    std::vector<double> factors;
    std::vector<int> triangles;

    //for(double width : std::vector<double> {0.5, 0.4, 0.3, 0.2, 0.1, 0.05}){
    for(double width : std::vector<double> {20, 15, 10, 7.5, 5, 2.5, 1}){
        for(double fac : std::vector<double> {1, 0.75, 0.5, 0.25, 0.1}){ //factor should be at least 0.5
            std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);
            const GridView gridView = grid->leafGridView();

            flowFragment f1 = {{0, 0.5}, {1, 0.5}, width};
            flowWithFragments(grid, f1, fac);
            int c1 = 0;
            for(const auto& element : elements(gridView)){
                c1++;
            }
            std::cout << "factor " << fac << ": " << c1 << " triangles" << std::endl;
            factors.push_back(fac);
            triangles.push_back(c1);
            w.push_back(width);
            std::string s = "test_" + std::to_string(fac);

            VTKWriter<GridView> vtkWriter(gridView);
            vtkWriter.write(s);
        }
    }

    for(int i = 0; i < factors.size(); i++){
        std::cout <<"width: " << w[i] << ", Factor: " << factors[i] << ", triangles: " << triangles[i] << std::endl;
    }
    

   

    // Write grid to file
    //VTKWriter<GridView> vtkWriter(gridView);
    //vtkWriter.write("refined_grid_borderInTriangle");
    
}
