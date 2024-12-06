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
#include <typeinfo>

/*
for every triangle checks if corners are above, below or inbetween the edges --> if corners are on different sides border within trinagle

*/

using namespace Dune;
constexpr int dim = 2; // Grid and world dimension


double distance(FieldVector<double, dim> p1, FieldVector<double, dim> p2){
    return std::abs((p1 - p2).two_norm());
}






void horizontalFlow(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, double width, FieldVector<double, dim> pos ){

    //TODO check if pos is at edge of map --> special case

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    double epsilon = 0.001;

    double yBorderU = pos[1] + width/2;
    double yBorderD = pos[1] - width/2;

    double xRangeL = pos[0] - 0.5;
    double xRangeR = pos[0] + 0.5;
    
    
    int c=0;
    while( c<100){

        c++;
        
        for (const auto& element : elements(gridView)){


            bool between = false;
            bool upper = false;
            bool under = false;

            bool continueFor = false;
            int xOut = 0;

            for (int i = 0; i < 3; i++){
                if(element.geometry().corner(i)[0] < xRangeL || element.geometry().corner(i)[0] > xRangeR){
                    xOut ++;
                }
                if (element.geometry().corner(i)[1] > yBorderU) upper = true;
                else if (element.geometry().corner(i)[1] < yBorderD) under = true;
                //else if (element.geometry().corner(i)[1] == yBorderD || element.geometry().corner(i)[1] == yBorderU){ //border on edge of triangle
                //    grid->mark(-1, element); 
                //    continueFor = true;
                //    break;
                //}
                else between = true;
            }

            if(xOut == 3) continue;

            //if(continueFor) continue;

            if(between + upper + under <= 1) continue;  //border not within triangle

            double epsilon2 = 3e-2;
            double dist01 = distance(element.geometry().corner(0), element.geometry().corner(1));
            double dist02 = distance(element.geometry().corner(0), element.geometry().corner(2));
            double dist12 = distance(element.geometry().corner(1), element.geometry().corner(2));
            
            double minDist = std::min({dist01, dist02, dist12});

            if (minDist < epsilon2) continue;

            int hypo = 0;
            Dune::FieldVector<double, 2> v;
            double yCenterHeight;
            double maxDist = std::max({dist01, dist02, dist12});
            if (maxDist == dist01){ 
                v = element.geometry().corner(0) - element.geometry().corner(1);
                if(element.geometry().corner(2)[1] == element.geometry().corner(0)[1]) {
                    yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(0)[1]) / 2;}
                else yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(1)[1]) / 2;
                hypo = 0;}
            if (maxDist == dist02) {
                v = element.geometry().corner(0) - element.geometry().corner(2);
                if (element.geometry().corner(1)[1] == element.geometry().corner(2)[1]){
                    yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(2)[1]) / 2;}
                else yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(0)[1]) / 2;
                hypo = 1;}
            if (maxDist == dist12) {
                v = element.geometry().corner(2) - element.geometry().corner(1);
                if (element.geometry().corner(0)[1] == element.geometry().corner(1)[1]){
                    yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(1)[1]) / 2;
                }
                else yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(2)[1]) / 2;
                hypo = 2;}

            Dune::FieldVector<double, 2> u1{1.0, 0.0};
            Dune::FieldVector<double, 2> u2{0, 1.0};
            double dotProduct1 = v[0] * u1[0] + v[1] * u1[1];
            double dotProduct2 = v[0] * u2[0] + v[1] * u2[1];
            if (std::abs(dotProduct1) < 1e-8 || std::abs(dotProduct2) < 1e-8){
                grid->mark(1, element);//hypothenuse vertical or horizontal
                continue;
            }

            //double yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(2)[1]) / 2 + (element.geometry().corner(0)[1] - element.geometry().corner(1)[1]) / 2;

            if(yBorderU - yCenterHeight < epsilon || yBorderD - yCenterHeight < 3e-3) continue; // border close to middle of triangle (in right angle)

            grid->mark(1,element);
        }
       
        grid->preAdapt();
        grid->adapt();
        grid->postAdapt();            
    }  
}



// { grid_setup_begin }
int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    

    // Start with a structured grid
    const std::array<unsigned, dim> n = {6, 6};
    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {6, 6};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

  

    horizontalFlow(grid, 0.1, FieldVector<double, dim> {{4, 1.5}});
    horizontalFlow(grid, 0.3, FieldVector<double, dim> {{3, 1.5}});
    

    // Write grid to file
    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("refined_grid_borderInTriangle");
    
}