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


    std::ofstream outFile("test_log_diag.txt");
    std::ofstream outTable("test_table_diag.txt");

    
    outTable << "\\begin{tabular}{ c |" ;
    std::vector<double> minSizeFactors;
    for( int i = 1; i < 8; i++){
        double msf = i/10.0;
        minSizeFactors.push_back(msf);
        outTable << " c";
    }
    outTable << " }" << std::endl;
    for(double msf :minSizeFactors){
        outTable << " & " << msf;
    }
     outTable << " \\\\" << std::endl << "\\hline" << std::endl;

            
    for(int w = 1; w<6; w++){
        double width = w/10.0;

        outTable << width;
        for(double msf : minSizeFactors){
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

            const GridView gridView01 = grid01->leafGridView();
            
            flowFragment f01 = {Dune::FieldVector<double, 2>{0, 0}, Dune::FieldVector<double, 2>{1, 1}, width};
            std::vector<flowFragment> fragments01 = {f01};

            refineGridwithFragments(grid01, fragments01, msf, {1.0, 1.0});
        
                
            std::string name = "1x1grid_diag_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }

            outFile << "Breite: " << width << ", minSIzeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable << " & " << counter;
            // Write grid to file

            VTKWriter<GridView> vtkWriter(gridView01);
            vtkWriter.write(name);
        }
        outFile << std::endl;
        outTable << " \\\\" << std::endl;
    }

     
    outTable.close();
    outFile.close();
}