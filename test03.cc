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



using namespace Dune;
constexpr int dim = 2; // Grid and world dimension


double signedDistance(FieldVector<double, dim> p1, FieldVector<double, dim> p2){
    return (p1 - p2).two_norm();
}

void horizontalFlow(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, double width, FieldVector<double, dim> pos ){

    //TODO check if pos is at border --> special case

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    double epsilon = 0.01;

    constexpr size_t num_points = 11; 
    double dist = 1/num_points;
    std::array<FieldVector<double, dim>, num_points> points;
    double y = pos[1];
    for (size_t i = 0; i < num_points; ++i) {
        points[i] = {pos[0]-0.5 + i * 0.1, y}; // x-Values: center left quadrat to center right quardrat
    }
    
    for (auto point : points)
    {
        double minDistBorderU = 1;
        double minDistBorderL = 1;
        int c=0;
        while((minDistBorderL > epsilon || minDistBorderU > epsilon) && c<50){
    
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
        
            if(minDistBorderL < epsilon && minDistBorderU < epsilon) break;

            for (const auto& element : elements(gridView)){
                if(std::abs(std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) - minDistBorderU) < 0.000001 ||
                    std::abs(std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) - minDistBorderL) < 0.00001){
                    grid->mark(1, element);
                }
            }    
            grid->preAdapt();
            grid->adapt();
            grid->postAdapt();         
        }
    }
}

void horizontalFlow2(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, double width, FieldVector<double, dim> pos ){

    //TODO check if pos is at border --> special case

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    double epsilon = 0.001;

    constexpr size_t num_points = 11; 
    double dist = 1/num_points;
    std::array<FieldVector<double, dim>, num_points> points;
    double y = pos[1];
    for (size_t i = 0; i < num_points; ++i) {
        points[i] = {pos[0]-0.5 + i * 0.1, y}; // x-Values: center left quadrat to center right quardrat
    }
    
    for (auto point : points)
    {
        double minDistBorderU = 1;
        double minDistBorderL = 1;
        int c=0;
        while((minDistBorderL > epsilon || minDistBorderU > epsilon) && c<50){
    
            c++;

            double minDistBorderU = 1;
            double minDistBorderL = 1;

            //std::vector<unsigned int> markedIndicesU;
            //std::vector<unsigned int> markedIndicesL;

            auto minElementU1 = gridView.end<0>();
            auto minElementU2 = gridView.end<0>();
            auto minElementU3 = gridView.end<0>();
            auto minElementU4 = gridView.end<0>();
            auto minElementL1 = gridView.end<0>();
            auto minElementL2 = gridView.end<0>();
            auto minElementL3 = gridView.end<0>();
            auto minElementL4 = gridView.end<0>();
     


            
            for (auto it = gridView.begin<0>(); it != gridView.end<0>(); ++it){
                //auto* elemPtr = &element;
                if(std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, it->geometry().center())) < minDistBorderU){
                    minDistBorderU = std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, it->geometry().center()));
                    minElementU1 = it;
                    minElementU2 = gridView.end<0>();
                    minElementU3 = gridView.end<0>();
                    minElementU4 = gridView.end<0>();
                    //minElementU = it;
                    //markedIndicesU = {(gridView.indexSet().index(*it))};
                }
                if(std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, it->geometry().center())) == minDistBorderU){
                    //markedIndicesU.push_back(gridView.indexSet().index(*it));
                    if (gridView.end<0>() == minElementU1)  minElementU1 = it;
                    else if (gridView.end<0>() == minElementU2)   minElementU2 = it;
                    else if (gridView.end<0>() == minElementU3)   minElementU3 = it;
                    else if (gridView.end<0>() == minElementU4)   minElementU4 = it;
                    else  std::cout << "5th min element" << std::endl;
                }


                if(std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, it->geometry().center())) < minDistBorderL){
                    minDistBorderL = std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, it->geometry().center()));
                    //minElementL = it;
                    minElementL1 = it;
                    minElementL2 = gridView.end<0>();
                    minElementL3 = gridView.end<0>();
                    minElementL4 = gridView.end<0>();
                    //markedIndicesL = {(gridView.indexSet().index(*it))};
                }
                if(std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, it->geometry().center())) == minDistBorderL){
                    //markedIndicesL.push_back(gridView.indexSet().index(*it));
                    if (gridView.end<0>() == minElementL1)  minElementL1 = it;
                    else if (gridView.end<0>() == minElementL2)  minElementL2 = it;
                    else if (gridView.end<0>() == minElementL3)   minElementL3 = it;
                    else if (gridView.end<0>() == minElementL4)   minElementL4 = it;
                    else  std::cout << "5th min element" << std::endl;
                }
            }
        
            if(minDistBorderL < epsilon && minDistBorderU < epsilon) break;

            if(minElementU1 != gridView.end<0>()) grid->mark(1, *minElementU1);
            if(minElementU2 != gridView.end<0>()) grid->mark(1, *minElementU2);
            if(minElementU3 != gridView.end<0>()) grid->mark(1, *minElementU3);
            if(minElementU4 != gridView.end<0>()) grid->mark(1, *minElementU4);
            if(minElementL1 != gridView.end<0>()) grid->mark(1, *minElementL1);
            if(minElementL2 != gridView.end<0>()) grid->mark(1, *minElementL2);
            if(minElementL3 != gridView.end<0>()) grid->mark(1, *minElementL3);
            if(minElementL4 != gridView.end<0>()) grid->mark(1, *minElementL4);

            //for (int index : markedIndicesU) {
            //    //auto element = gridView.indexSet().entity<0>(index); // Element abrufen
            //    //std::cout << "Markiertes Element mit Index " << index
            //    //        << ", Volumen: " << element.geometry().volume() << std::endl;
            //    //Dune::MarkingStrategy::BaseMarkingStrategy<Grid> marker;
            //    //marker.mark(*index, Dune::MarkingStrategy::refine);
            //    grid->mark(1, element);
            //}
            //for (int index : markedIndicesL) {
            //    auto element = gridView.indexSet().entity<0>(index); // Element abrufen
            //    //std::cout << "Markiertes Element mit Index " << index
            //    //        << ", Volumen: " << element.geometry().volume() << std::endl;
            //    //Dune::MarkingStrategy::BaseMarkingStrategy<Grid> marker;
            //    //marker.mark(*index, Dune::MarkingStrategy::refine);
            //    grid->mark(1, element);
            //}

            //for (const auto& element : elements(gridView)){
            //    if(std::abs(std::abs(signedDistance(point+FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) - minDistBorderU) < 0.000001 ||
            //        std::abs(std::abs(signedDistance(point-FieldVector<double, dim>{{0, width/2}}, element.geometry().center())) - minDistBorderL) < 0.00001){
            //        grid->mark(1, element);
            //    }
            //}    
            grid->preAdapt();
            grid->adapt();
            grid->postAdapt();         
        }
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

  

    
    horizontalFlow(grid, 0.3, FieldVector<double, dim> {{3, 1.5}});
    horizontalFlow(grid, 0.1, FieldVector<double, dim> {{4, 1.5}});
    horizontalFlow2(grid, 0.1, FieldVector<double, dim> {{4, 2.5}});

    // Write grid to file

    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("refined_grid");
    
}