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
    const std::array<unsigned, 2> n = {4, 3};
    const FieldVector<double, 2> lower = {0, 0};
    const FieldVector<double, 2> upper = {360, 270};


    std::ofstream outFile("test_log_diag.txt");
    std::ofstream outTable("test_table_diag.txt");

    std::vector<int> amountElements1(300, 0);
    std::vector<int> amountElements2(300, 0);
    std::vector<int> amountElements3(300, 0);

    
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
            
            flowFragment f01 = {Dune::FieldVector<double, 2>{0, 360}, Dune::FieldVector<double, 2>{360, 0}, width, width, 1, 1000, 3};
            std::vector<flowFragment> fragments01 = {f01};

            std::vector<std::vector<flowFragment>> rivers = {fragments01, fragments01, fragments01, fragments01};

            refineGridwithFragments(grid01, rivers,{4, 3}, {90, 90}, 50, msf, true);
        
                
            std::string name = "4x3grid_diag_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter1=0;
            int counter2 = 0;
            int counter3 = 0;
            bool cont1 = false;
            bool cont2 = false;
            bool cont3 = false;
            for(auto element : elements(gridView01)){
                cont1 = false;
                cont2 = false;
                cont3 = false;
                for(int i = 0; i < 3; i++){
                    if(element.geometry().corner(i)[0] < 180 || element.geometry().corner(i)[0] >270 || element.geometry().corner(i)[1] < 90 ||element.geometry().corner(i)[1] >180){
                        cont1 = true;
                    }
                    if(element.geometry().corner(i)[0] < 90 || element.geometry().corner(i)[0] >180 || element.geometry().corner(i)[1] < 90 ||element.geometry().corner(i)[1] >180){
                        cont2 = true;
                    }
                    if(element.geometry().corner(i)[0] < 90 || element.geometry().corner(i)[0] >180 || element.geometry().corner(i)[1] < 0 ||element.geometry().corner(i)[1] >90){
                        cont3 = true;
                    }
                }
                //std::cout << cont1 << cont2 << cont3;
                //std::cout << element.geometry().corner(0) << std::endl;
                if(!cont1) {
                    counter1++;
                }
                if(!cont2) counter2++;
                if(!cont3) {
                    //std::cout << element.geometry().corner(0) << std::endl;
                    counter3++;}
            }



            //outFile << "Breite: " << width << ", minSizeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable << " & " << counter1 << "+" << 2*counter2 << "/" << 2*counter3;
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