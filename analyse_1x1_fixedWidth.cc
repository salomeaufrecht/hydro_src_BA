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

//     realCellSize = {90.0, 90.0}; 
// einfügen in code sonst geht es nicht
 


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
    const FieldVector<double, 2> upper = {90, 90};

    std::ofstream outFile("test_log_fixed_width.txt");
    std::ofstream outTable("test_table_fixed_width.txt");

    
    outTable << "\\begin{tabular}{ c |" ;
    std::vector<double> minSizeFactors = {2, 4, 5, 10};
    std::vector<double> widths = {9.6, 46.3, 65.5, 90};

    for( int i = 1; i < minSizeFactors.size(); i++){
        outTable << " c";
    }
    outTable << " }" << std::endl;
    outTable << "Breite:";
    for(double w : widths){
        outTable << " & " <<   w;
    }
     outTable << " \\\\" << std::endl << "\\hline" << std::endl<< "Anzahl Elemente:";

    
    for(int i = 0; i<4; i++){
        double width = widths[i];
        
            
        std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

        const GridView gridView01 = grid01->leafGridView();
        flowFragment f01 = {Dune::FieldVector<double, 2>{45, 0}, Dune::FieldVector<double, 2>{45, 90}, width, width, 1, 1000, 2};
        std::vector<flowFragment> fragments01 = {f01};

        std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};


        refineGridwithFragments(grid01, rivers,{1, 1}, {90, 90}, 50, minSizeFactors[i], true);
    
            
        std::string name = "1x1grid_hor_fixed_w_"+ std::to_string(width) + "_msf_" + std::to_string( minSizeFactors[i]);

        int counter=0;
        for(auto element : elements(gridView01)){
            counter++;
        }

        outFile << "Breite: " << width << ", minSIzeFactor: " <<  minSizeFactors[i] << ": " << counter << " elements" << std::endl;
        outTable << " & " << counter;
        // Write grid to file

        VTKWriter<GridView> vtkWriter(gridView01);
        vtkWriter.write(name);
        
        outFile << std::endl;
        //outTable << " \\\\" << std::endl;
    }

    outTable << std::endl << std::endl;

    std::cout << "Diagonal" << std::endl;


    for( int i = 1; i < minSizeFactors.size(); i++){
        outTable << " c";
    }
    outTable << " }" << std::endl;
    outTable << "Breite:";
    for(double w : widths){
        outTable << " & " <<   w;
    }
     outTable << " \\\\" << std::endl << "\\hline" << std::endl<< "Anzahl Elemente:";

    
    for(int i = 0; i<4; i++){
        double width = widths[i];
        
            
        std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

        const GridView gridView01 = grid01->leafGridView();
        flowFragment f01 = {Dune::FieldVector<double, 2>{0, 0}, Dune::FieldVector<double, 2>{90, 90}, width, width, 1, 1000, 2};
        std::vector<flowFragment> fragments01 = {f01};

        std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};


        refineGridwithFragments(grid01, rivers,{1, 1}, {90, 90}, 50, minSizeFactors[i], false);
    
            
        std::string name = "1x1grid_diag_fixed_w_"+ std::to_string(width) + "_msf_" + std::to_string( minSizeFactors[i]);

        int counter=0;
        for(auto element : elements(gridView01)){
            counter++;
        }

        outFile << "Breite: " << width << ", minSizeFactor: " <<  minSizeFactors[i] << ": " << counter << " elements" << std::endl;
        outTable << " & " << counter;
        // Write grid to file

        VTKWriter<GridView> vtkWriter(gridView01);
        vtkWriter.write(name);
        
        outFile << std::endl;
        //outTable << " \\\\" << std::endl;
    }

    outTable << std::endl << std::endl;

     
    outTable.close();
    outFile.close();
}