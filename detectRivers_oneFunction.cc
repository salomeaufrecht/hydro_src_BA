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
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-pdelab includes
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/variablefactories.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/solver/newton.hh>

// dune-vtk includes
#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/writers/vtkrectilineargridwriter.hh>
#include <dune/vtk/writers/vtkstructuredgridwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>
#include <dune/vtk/datacollectors/spdatacollector.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>


// include stuff from dune-hydro
#include <dune/hydro/nonlineardiffusionfv.hh>
#include <dune/hydro/nonlineardiffusionfv_velocity.hh>
#include <dune/hydro/geotiffreader.hh>
#include <dune/hydro/netcdfreader.hh>
#include <dune/hydro/rasterdataset.hh>
#include <dune/vtk/pvdwriter.hh>

#include "flowFunctions.hh"



int main(int argc, char **argv)
{
    try
    {
        // Maybe initialize MPI
        Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
        std::cout << "This is dune-hydro." << std::endl;
        if (Dune::MPIHelper::isFake)
        std::cout << "This is a sequential program." << std::endl;
        else
        std::cout << "I am rank " << helper.rank() << " of " << helper.size()
                    << " processes!" << std::endl;

        // register GDAL drivers once at the beginning
        GDALAllRegister();


        // open ini file
        Dune::ParameterTree ptree;
        Dune::ParameterTreeParser ptreeparser;
        ptreeparser.readINITree("nakhon.ini", ptree);
        ptreeparser.readOptions(argc, argv, ptree);

        // Nakhonsitamarat
        if (true)
        {
        // read data files
        GeoTIFFImage<GInt16> image("n00e090_con.tif", "", 1, 2);
        GeoTIFFImage<GUInt32> aimage("N00E090_acc.tif", "", 1, 2);
        GeoTIFFImage<GByte> dimage("n00e090_dir.tif", "", 1, 2);

        double dx = std::abs(image.dLong());
        double dy = std::abs(image.dLat());

        std::array<int, 2> N;
        N[0] = 50; //1800
        N[1] = 50; //1200
        std::array<double, 2> H;
        H[0] = 90; 
        H[1] = 90; 
        Dune::FieldVector<double, 2> L;
        L[0] = N[0] * H[0];
        L[1] = N[1] * H[1];

      auto elevation_raster    = RasterDataSet<float>           (99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); //original: 99.0, 8.0 (99 breite, 8 höhe)
      auto accumulation_raster = RasterDataSet<float>           (99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      auto direction_raster    = RasterDataSet<unsigned char>   (99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
        //auto elevation_raster    = RasterDataSet<float>(99.2 + 0.5 * dx, 8.4 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); //original: 99.0, 8.0 (99 breite, 8 höhe)
        //auto accumulation_raster = RasterDataSet<float>(99.2 + 0.5 * dx, 8.4 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
        //auto direction_raster    = RasterDataSet<unsigned char>(99.2 + 0.5 * dx, 8.4 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
        elevation_raster.paste(image);
        accumulation_raster.paste(aimage);
        direction_raster.paste(dimage);

        // create grid
        typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
        using GridView = Grid::LeafGridView;

        const std::array<unsigned, 2> n = {N[0], N[1]};
        const Dune::FieldVector<double, 2> lower = {0.5 * H[0], 0.5 * H[1]};
        const Dune::FieldVector<double, 2> upper = {L[0] - 0.5 * H[0], L[1] - 0.5 * H[1]};

        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);
        GridView gridView = grid->leafGridView();

        // compute rivers and add them to the map
        std::vector<double> height = addRiversToMap(grid, H, N, accumulation_raster, direction_raster, elevation_raster);

        

        // write result to file
        Dune::VTKWriter<GridView> vtkWriter(gridView);
        vtkWriter.addVertexData(height, "height");
        vtkWriter.write("addRiversToMap");

    }
      return 0;
  }
  catch (Dune::Exception &e)
  {
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (GeoTIFFReaderException &e)
  {
    std::cerr << "geotiff error: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
