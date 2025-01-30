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
    const std::array<unsigned, 2> n = {5, 1};
    const FieldVector<double, 2> lower = {0, 0};
    const FieldVector<double, 2> upper = {5, 1};

    std::ofstream outFile("cells_around_log.txt");
    std::ofstream outTable("cells_around_table.txt");

    
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

    std::vector<int> amountElements(300, 0);    
    for(int w = 1; w<21; w++){
        double width = w/40.0;

        //outTable << width;
        for(double msf : minSizeFactors){

            const std::array<unsigned, 2> n = {1, 1};
            const FieldVector<double, 2> lower = {0, 0};
            const FieldVector<double, 2> upper = {1, 1};
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

            const GridView gridView01 = grid01->leafGridView();
            flowFragment f01 = {Dune::FieldVector<double, 2>{0, 0.5}, Dune::FieldVector<double, 2>{1, 0.5}, width};
            std::vector<flowFragment> fragments01 = {f01};


            refineGridwithFragments(grid01, fragments01, msf, {1.0, 1.0});
        

            //std::string name = "1x1grid_hor_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }
            std::cout << "here " << w-1 << " " << int(msf*10)-1 << std::endl;
            amountElements[w*10+int(msf*10)] = counter;
            //amountElements[0][0] = counter;

            //outFile << "Breite: " << width << ", minSIzeFactor: " << msf << ": " << counter << " elements" << std::endl;
            //outTable << " & " << counter;
            // Write grid to file

            //VTKWriter<GridView> vtkWriter(gridView01);
            //vtkWriter.write(name);
        }
    }

            
    for(int w = 1; w<21; w++){
        double width = w/40.0;

        outTable << width;
        for(double msf : minSizeFactors){
            
            std::shared_ptr<Grid> grid01 = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

            const GridView gridView01 = grid01->leafGridView();
            flowFragment f01 = {Dune::FieldVector<double, 2>{2.5, 0}, Dune::FieldVector<double, 2>{2.5, 1}, width};
            std::vector<flowFragment> fragments01 = {f01};


            refineGridwithFragments(grid01, fragments01, msf, {1.0, 1.0});
        
                
            std::string name = "1x3grid_hor_w_"+ std::to_string(width) + "_msf_" + std::to_string(msf);

            int counter=0;
            for(auto element : elements(gridView01)){
                counter++;
            }
            counter -= amountElements[w*10+int(msf*10)] + 8;

            outFile << "Breite: " << width << ", minSIzeFactor: " << msf << ": " << counter << " elements" << std::endl;
            outTable << " & " << counter;
            // Write grid to file

            VTKWriter<GridView> vtkWriter(gridView01);
            vtkWriter.write(name);
        }
        outFile << std::endl;
        outTable << " \\\\" << std::endl;
    }
    outTable << "\\end{tabular}";

     
    outTable.close();
    outFile.close();
}