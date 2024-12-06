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
constexpr int dim = 2; // Grid and world dimension

double signedDistance(FieldVector<double, dim> p1, FieldVector<double, dim> p2){
    return (p1 - p2).two_norm();
}

// { grid_setup_begin }
int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    
    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;


    // Start with a structured grid
    const std::array<unsigned, dim> n = {3, 3};
    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {3, 3};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

  

    double epsilon = 0.01;
    double width = 0.2;

    constexpr size_t num_points = 11; 
    std::array<FieldVector<double, dim>, num_points> points;
    double y = 1.5;
    for (size_t i = 0; i < num_points; ++i) {
        points[i] = {1.0 + i * 0.1, y}; // x-Werte von 1.0 bis 2.0
    }
    

    for (auto point : points)
    {
        double minDistBorderU = 1;
        double minDistBorderL = 1;
        int c=0;
        while((minDistBorderL > epsilon || minDistBorderU > epsilon) && c<50){
            std::cout << "\n Upper: " << minDistBorderU << std::endl;
            std::cout << "Lower: " << minDistBorderL << std::endl;
            c++;

            minDistBorderU = 1;
            minDistBorderL = 1;
            
            for (const auto& element : elements(gridView)){
                if(std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) < minDistBorderU){
                    minDistBorderU = std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, element.geometry().center()));
                }
                if(std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) < minDistBorderL){
                    minDistBorderL = std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, element.geometry().center()));
                }
            }

            std::cout << "\n Upper1: " << minDistBorderU << std::endl;
            std::cout << "Lower1: " << minDistBorderL << std::endl;
            std::cout << (minDistBorderL > epsilon) << (minDistBorderU > epsilon) << (c<50) << std::endl;

            if(minDistBorderL < epsilon && minDistBorderU < epsilon) break;


            for (const auto& element : elements(gridView)){
                if(std::abs(std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) - minDistBorderU) < 0.000001 ||
                    std::abs(std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) - minDistBorderL) < 0.00001){
                    std::cout << "mark" << std::endl;
                    grid->mark(1, element);
                }
            }
            
            grid->preAdapt();
            grid->adapt();
            grid->postAdapt();

            std::cout << "\n Upper2: " << minDistBorderU << std::endl;
            std::cout << "Lower2: " << minDistBorderL << std::endl;
            std::cout << (minDistBorderL > epsilon) << (minDistBorderU > epsilon) << (c<50) << "\n" << std::endl;
        
        
        }
   
    }
    

    // Write grid to file

    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("refined_grid");
    
}