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
const double pixelSize = 90;

struct flowFragment
    {
        FieldVector<double, dim> start;
        FieldVector<double, dim> end;
        double widthStart;
        double widthEnd = -1;
    };


double distance(FieldVector<double, dim> p1, FieldVector<double, dim> p2){
    return std::abs((p1 - p2).two_norm());
}

double squaredDistance(FieldVector<double, dim> p1, FieldVector<double, dim> p2){
    FieldVector<double, dim> diff = (p1 - p2);
    return diff[0]*diff[0] + diff[1]*diff[1];
}



void flow(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, FieldVector<double, dim> start, FieldVector<double, dim> end, double widthStart, double widthEnd = -1, double minSizeFactor = 0.5 ){



    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    FieldVector<double, dim> flowVector = end - start;
    FieldVector<double, dim> normal_flowVector = { -flowVector[1]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1]), 
                                                    flowVector[0]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1])} ;

    if (widthEnd < 0) widthEnd = widthStart;

    FieldVector<double, dim> start1 = start + widthStart * normal_flowVector;
    FieldVector<double, dim> start2 = start - widthStart * normal_flowVector;
    FieldVector<double, dim> end1 = end + widthEnd * normal_flowVector;
    FieldVector<double, dim> end2 = end - widthEnd * normal_flowVector;

    FieldVector<double, dim> b1 = end1 - start1;
    FieldVector<double, dim> b2 = end2 - start2;



    double minSize =  minSizeFactor*std::min(widthStart, widthEnd);
   
    
    int c=0;
    bool change = true;
    while( c<100 && change){
        change = false;
        c++;  
        for (const auto& element : elements(gridView)){

            bool between = false;
            bool out1 = false;
            bool out2 = false;

            bool continueFor = false;
            int xOut = 0;
            int yOut = 0;

            for (int i = 0; i < 3; i++){
                auto cornerI = element.geometry().corner(i);
                if(cornerI[0] < std::min({start1[0], start2[0], end1[0], end2[0]})) xOut--; //test if trangle is outside of desired range
                else if (cornerI[0] > std::max({start1[0], start2[0], end1[0], end2[0]}))xOut++;
                if(cornerI[1] < std::min({start1[1], start2[1], end1[1], end2[1]})) yOut--;
                else if (cornerI[1] > std::max({start1[1], start2[1], end1[1], end2[1]})) yOut ++;
    
                FieldVector<double, dim> vec1 = cornerI - start1;
                FieldVector<double, dim> vec2 = cornerI - start2;
                double cross1 = b1[0] * vec1[1] - b1[1] * vec1[0];
                double cross2 = b2[0] * vec2[1] - b2[1] * vec2[0];
            
                if (std::abs(cross1) < 1e-5 ||std::abs( cross2) < 1e-5){ //border on edge of triangle
                    grid->mark(-1, element); 
                    continueFor = true;
                    change = true;
                    break;
                }
                
                if(cross1 > 0) out1 = true;
                else if(cross2 < 0) out2 = true;
                else between = true;
            }

            
            if(std::abs(xOut) == 3 || std::abs(yOut) == 3 || continueFor || between + out1 + out2 <= 1) continue; // all 3 corners outside x/y boundry/border on corner/border not within triangle
    

            double dist01 = squaredDistance(element.geometry().corner(0), element.geometry().corner(1));
            double dist02 = squaredDistance(element.geometry().corner(0), element.geometry().corner(2));

            
            if (std::min({dist01, dist02}) < minSize*minSize) continue; //small enough

            Dune::FieldVector<double, 2> hypo;
            double yCenterHeight;

            if (dist01 == dist02){ //dist12
                hypo = element.geometry().corner(2) - element.geometry().corner(1);
                if (element.geometry().corner(0)[1] == element.geometry().corner(1)[1]){
                    yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(1)[1]) / 2;
                }
                else yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(2)[1]) / 2;
                }
            if ( dist01 > dist02){  //calc hypothenuse v; y value middle of height of triangle
                hypo = element.geometry().corner(0) - element.geometry().corner(1);
                if(element.geometry().corner(2)[1] == element.geometry().corner(0)[1]) {
                    yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(0)[1]) / 2;}
                else yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(1)[1]) / 2;
                }
            else {
                hypo = element.geometry().corner(0) - element.geometry().corner(2);
                if (element.geometry().corner(1)[1] == element.geometry().corner(2)[1]){
                    yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(2)[1]) / 2;}
                else yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(0)[1]) / 2;
                }
           

            double epsilonAngle = 1e-8;


            //TODO: till which angle? 
            if((std::abs(b1[0]) < epsilonAngle || std::abs(b1[1]) < epsilonAngle) && ((std::abs(b2[0]) < epsilonAngle || std::abs(b2[1]) < epsilonAngle))){ //flow (and border) horizontal or vertical (dot product with {1,0}/{0,1})
                if (std::abs(hypo[0]) < 1e-8 || std::abs(hypo[1]) < 1e-8){ //hypothenuse vertical or horizontal
                    grid->mark(1, element);//hypothenuse vertical or horizontal
                    change = true;
                    continue;
                } 
                double epsilon3 = 1e-2;
                if(std::abs(start1[1] - yCenterHeight) < epsilon3 || std::abs(start2[1] - yCenterHeight) < epsilon3) continue; // border close to middle of triangle (in right angle)

            }
            //else if(std::abs(flowVector[0] + flowVector[1]) < 1e-8 || std::abs(flowVector[0] - flowVector[1]) < 1e-8){ //flow diagonal
            //    if (std::abs(hypo[0] + hypo[1]) < 1e-8 || std::abs(hypo[0] - hypo[1]) < 1e-8){ //hypothenuse diagonal
            //        grid->mark(1, element);//hypothenuse vertical or horizontal
            //        change = true;
            //        continue;
            //    }
            //    // TODO: continue if right position (analog to middle in other case)
            //}
        
            change = true;
            grid->mark(1,element);
        }
       
        grid->preAdapt();
        grid->adapt();  
        grid->postAdapt();            
    }  
}

void flowWithFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f, double minSizeFactor=0.6){
    if(f.widthStart>1) flow(grid, f.start, f.end, f.widthStart/pixelSize, f.widthEnd/pixelSize);
    else flow(grid, f.start, f.end, f.widthStart, f.widthEnd);
}


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
    flowFragment f3 = {{1.5, 1.5}, {1.5, 6}, 0.05, 0.5};

    flowWithFragments(grid, f1);
    
    std::vector<flowFragment> fragments = {f1, f2, f3};

    for(auto f : fragments){
        flowWithFragments(grid, f);
    }


    

    // Write grid to file
    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.write("refined_grid_borderInTriangle");
    
}
