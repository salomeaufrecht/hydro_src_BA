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

    std::vector<int> amountElements(300, 0);    


    std::ofstream outFile("test_log.txt");
    std::ofstream outTable("test_table.txt");

    
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

    int pos=0;
    for(int w = 1; w<90; w+=5){
        double width = w; ///40.0;

        outTable << width;
        for(double msf : minSizeFactors){
            pos++;
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

            const GridView gridView01 = grid01->leafGridView();
            flowFragment f01 = {Dune::FieldVector<double, 2>{45, 0}, Dune::FieldVector<double, 2>{45, 90}, width, width, 1, 1000, 2};
            std::vector<flowFragment> fragments01 = {f01};

            std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};


            refineGridwithFragments(grid01, rivers,{1, 1}, {90, 90}, 50, msf, true);
        
                
            std::string name = "1x1grid_hor_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }

            amountElements[pos] = counter;


            outFile << "Breite: " << width << ", minSIzeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable << " & " << counter;
            // Write grid to file

            VTKWriter<GridView> vtkWriter(gridView01);
            vtkWriter.write(name);
        }
        outFile << std::endl;
        outTable << " \\\\" << std::endl;
    }
    
    outTable << std::endl << std::endl << std::endl;

    const std::array<unsigned, 2> n_ = {3, 1};
    const FieldVector<double, 2> lower_ = {0, 0};
    const FieldVector<double, 2> upper_ = {270, 90};


    outTable << "\\begin{tabular}{ c |" ;
    for( int i = 0; i < minSizeFactors.size(); i++){
        outTable << " c";
    }
    outTable << " }" << std::endl;
    for(double msf :minSizeFactors){
        outTable << " & " << msf;
    }
     outTable << " \\\\" << std::endl << "\\hline" << std::endl;

    pos = 0;
    for(int w = 1; w<90; w+=5){
        double width = w; ///40.0;

        outTable << width;
        for(double msf : minSizeFactors){
            pos++;
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower_, upper_, n_);
            const GridView gridView01 = grid01->leafGridView();
            flowFragment f01 = {Dune::FieldVector<double, 2>{135, 0}, Dune::FieldVector<double, 2>{135, 90}, width, width, 1, 1000, 2};
            std::vector<flowFragment> fragments01 = {f01};

            std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};


            refineGridwithFragments(grid01, rivers,{3, 1}, {90, 90}, 50, msf, true);
        
                
            std::string name = "3x1grid_hor_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }

    

            outFile << "Breite: " << width << ", minSIzeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable << " & " <<amountElements[pos] << "+" << counter - amountElements[pos];
            // Write grid to file

            //VTKWriter<GridView> vtkWriter(gridView01);
            //vtkWriter.write(name);
        }
        outFile << std::endl;
        outTable << " \\\\" << std::endl;
    }



     
    outTable.close();
    outFile.close();
}