#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include <algorithm>
// C++ includes
#include <math.h>
#include <iostream>
// dune-common includes
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/timer.hh>
// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
// dune-grid includes
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/alugrid/grid.hh>


#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

using namespace Dune;

int main(int argc, char **argv){

   
    int n = 4;
    int m = 3;

    constexpr int dim = 2;
    using Grid = UGGrid<dim>;

    GridFactory<UGGrid<2> > factory;



    for (int i = 0; i < n; i++){
        for (int j = 0; j < m ; j++){
            factory.insertVertex({i, j});
        }
    }

    for (int i = 0; i < n-1; i++){
        for (int j = 0; j < m-1 ; j++){
            factory.insertElement(GeometryTypes::simplex(2), {i*m+j+1, i*m+j, (i+1)*m+j});
            factory.insertElement(GeometryTypes::simplex(2), {i*m+j+1, (i+1)*m+j+1, (i+1)*m+j});
        }
    }

    std::unique_ptr<UGGrid<2> > grid = factory.createGrid();

    // VTK-Writer initialisieren und Datei schreiben
    //Dune::VTKWriter<GridType::LeafGridView> vtkWriter(grid->leafGridView());
    //vtkWriter.write("test_grid_output");

    using GridView = Grid::LeafGridView;
    GridView gridView = grid->leafGridView();

    VTKWriter<GridView> vtkWriter( gridView);
    //vtkWriter.addVertexData(x, "solution");
    std::cout << (vtkWriter.write("test_grid_output"));


}