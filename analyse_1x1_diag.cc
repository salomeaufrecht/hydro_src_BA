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
    const FieldVector<double, 2> upper = {90, 90};


    std::ofstream outFile("test_log_diag.txt");
    std::ofstream outTable("test_table_diag.txt");

    std::vector<int> amountElements1(300, 0);

    
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

    int pos = 0;
    for(int w = 1; w<90; w+=5){
        double width = w; ///10.0;

        outTable << width;
        for(double msf : minSizeFactors){
            pos++;
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

            const GridView gridView01 = grid01->leafGridView();
            
            flowFragment f01 = {Dune::FieldVector<double, 2>{-90, 180}, Dune::FieldVector<double, 2>{180, -90}, width, width, 1, 1000, 3};
            std::vector<flowFragment> fragments01 = {f01};

            std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};

            refineGridwithFragments(grid01, rivers,{1, 1}, {90, 90}, 50, msf, true);
        
                
            std::string name = "1x1grid_diag_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }

            amountElements1[pos] = counter;

            outFile << "Breite: " << width << ", minSizeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable << " & " << counter;
            // Write grid to file

            VTKWriter<GridView> vtkWriter(gridView01);
            vtkWriter.write(name);
        }
        outFile << std::endl;
        outTable << " \\\\" << std::endl;
    }

    outTable << std::endl << std::endl << std::endl;

        const std::array<unsigned, 2> n2 = {2, 1};
    const FieldVector<double, 2> lower2 = {0, 0};
    const FieldVector<double, 2> upper2 = {180, 90};

    std::vector<int> amountElements2(300, 0);

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
        double width = w; ///10.0;

        //outTable << width;
        for(double msf : minSizeFactors){
            pos++;
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower2, upper2, n2);

            const GridView gridView01 = grid01->leafGridView();
            
            flowFragment f01 = {Dune::FieldVector<double, 2>{-90, 180}, Dune::FieldVector<double, 2>{180, -90}, width, width, 1, 1000, 3};
            std::vector<flowFragment> fragments01 = {f01};

            std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};

            refineGridwithFragments(grid01, rivers,{2, 1}, {90, 90}, 50, msf, true);
        
                
            //std::string name = "1x1grid_diag_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }

            amountElements2[pos] = counter - amountElements1[pos];

            //outFile << "Breite: " << width << ", minSizeFactor: " << msf << ": " << counter << " elements" << std::endl;
            //outTable << "+" << amountElements2[pos];
            // Write grid to file

            //VTKWriter<GridView> vtkWriter(gridView01);
            //vtkWriter.write(name);
        }
        //outFile << std::endl;
        //outTable << " \\\\" << std::endl;
    }

    const std::array<unsigned, 2> n4 = {2, 2};
    const FieldVector<double, 2> lower4 = {0, 0};
    const FieldVector<double, 2> upper4 = {180, 180};

    pos = 0;

    for(int w = 1; w<90; w+=5){
        double width = w; ///10.0;

        outTable << width;
        for(double msf : minSizeFactors){
            pos++;
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower4, upper4, n4);

            const GridView gridView01 = grid01->leafGridView();
            
            flowFragment f01 = {Dune::FieldVector<double, 2>{-90, 180}, Dune::FieldVector<double, 2>{180, -90}, width, width, 1, 1000, 3};
            std::vector<flowFragment> fragments01 = {f01};

            std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};

            refineGridwithFragments(grid01, rivers,{2, 2}, {90, 90}, 50, msf, true);
        
                
            std::string name = "4x4_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }

            //amountElements2[pos] = counter - amountElements1[pos];

            //outFile << "Breite: " << width << ", minSizeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable <<" & " << amountElements1[pos] << "+" << 2*amountElements2[pos] << "/" << 2*(counter - 2*amountElements2[pos] - amountElements1[pos]);
            // Write grid to file

            VTKWriter<GridView> vtkWriter(gridView01);
            vtkWriter.write(name);
        }
        //outFile << std::endl;
        outTable << " \\\\" << std::endl;
    }

     
    outTable.close();
    outFile.close();
}