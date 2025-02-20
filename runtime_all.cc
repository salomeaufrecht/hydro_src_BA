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
#include <dune/common/timer.hh>


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
      std::array<double, 2> H;
      H[0] = 90; 
      H[1] = 90;
      Dune::FieldVector<double, 2> L;

      std::array<int, 2> N;
      std::vector<double> size;
      std::vector<double> runtime_detectFragments;
      std::vector<double> runtime_refineGrid;
      std::vector<double> runtime_overallHeight;
      std::vector<double> runtime_applyFlowHeightFragments;

      for(int i = 5; i < 600; i+=8){
            
        N[0] = i;
        N[1] = i;
        L[0] = N[0] * H[0];
        L[1] = N[1] * H[1];
        std::cout << i << std::endl;

        //was 9.04, 8.11

      auto elevation_raster    = RasterDataSet<float>         (99.0 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); //original: 99.0, 8.0 (99 breite, 8 höhe)
      auto accumulation_raster = RasterDataSet<float>         (99.0 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      auto direction_raster    = RasterDataSet<unsigned char> (99.0 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      //auto elevation_raster    = RasterDataSet<float>         (99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); //original: 99.0, 8.0 (99 breite, 8 höhe)
      //auto accumulation_raster = RasterDataSet<float>         (99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      //auto direction_raster    = RasterDataSet<unsigned char> (99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      elevation_raster.paste(image);
      accumulation_raster.paste(aimage);
      direction_raster.paste(dimage);

        // create grid
        typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid_;
        using GridView = Grid_::LeafGridView;
          
        const std::array<unsigned, 2> n_ = {N[0], N[1]};
        const Dune::FieldVector<double, 2> lower = {0.5 * H[0], 0.5 * H[1]};
        const Dune::FieldVector<double, 2> upper = {L[0] - 0.5 * H[0], L[1] - 0.5 * H[1]};
        std::shared_ptr<Grid_> grid = Dune::StructuredGridFactory<Grid_>::createSimplexGrid(lower, upper, n_);
        const GridView gridView = grid->leafGridView();



        Dune::Timer timeDetectFragments;
        std::vector<std::vector<flowFragment>> rivers = detectFragments(accumulation_raster, direction_raster, H, N);
        double elapsedTimeDetectFragments = timeDetectFragments.stop();


        Dune::Timer timerRefine;  // Timer starten
        refineGridwithFragments(grid, rivers,N,H); 
        double elapsedTimeRefine = timerRefine.stop();  // Stoppe den Timer und bekomme die Zeit in Sekunden
        
        Dune::Timer timerOverallHeight;  // Timer starten
        std::vector<double> height = overallHeight(gridView, elevation_raster, H, N);
        double elapsedTimeOverallHeight = timerOverallHeight.stop();  // Stoppe den Timer und bekomme die Zeit in Sekunden

        Dune::Timer timerApplyFlowHeightFragments;  // Timer starten
        height= applyFlowHeightFragments(gridView, rivers, height, elevation_raster, H, N);
        double elapsedTimeApplyFlowHeightFragments = timerApplyFlowHeightFragments.stop();  // Stoppe den Timer und bekomme die Zeit in Sekunden

        size.push_back(i*i);
        runtime_detectFragments.push_back(elapsedTimeDetectFragments);
        runtime_refineGrid.push_back(elapsedTimeRefine);
        runtime_overallHeight.push_back(elapsedTimeOverallHeight);
        runtime_applyFlowHeightFragments.push_back(elapsedTimeApplyFlowHeightFragments);

        }

      std::cout << std::endl << "size = [ " << size[0] ;
      for (auto s : size){
          std::cout << ", " << s;
      }
      std::cout << "]" << std::endl;

      std::cout << "detectFragments = [" << runtime_detectFragments[0] ;
      for (auto r : runtime_detectFragments){
          std::cout << ", "<< r ;
      }
      std::cout << "]" << std::endl;

      std::cout << "refineGrid = [" << runtime_refineGrid[0] ;
      for (auto r : runtime_refineGrid){
          std::cout << ", "<< r ;
      }
      std::cout << "]" << std::endl;

      std::cout << "overallHeight = [" << runtime_overallHeight[0] ;
      for (auto r : runtime_overallHeight){
          std::cout << ", "<< r ;
      }
      std::cout << "]" << std::endl;

      std::cout << "applyFlowHeightFragments = [" << runtime_applyFlowHeightFragments[0] ;
      for (auto r : runtime_applyFlowHeightFragments){
          std::cout << ", "<< r ;
      }
      std::cout << "]" << std::endl;

      
      std::cout << std::endl << "plt.plot(size, detectFragmetents, \"-\", label = \"time detectFragments\")" << std::endl;
      std::cout << "plt.plot(size, refineGrid, \"-\", label=\"time refineGrid\")" << std::endl;
      std::cout << "plt.plot(size, overallHeight, \"-\", label = \"time overallHeight\" )" << std::endl;
      std::cout << "plt.plot(size, applyFlowHeightFragments, \"-\", label = \"time applyFlowHeightFragments\")" << std::endl;
      std::cout << "plt.ylabel(\"sec\")" << std::endl;
      std::cout << "plt.xlabel(\"size of map\")" << std::endl;
      std::cout << "plt.legend()" <<std::endl;
      return 0;
    }
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
