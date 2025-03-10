
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>

// dune-common includes
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
// dune-grid includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>

// dune-pdelab includes
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/backend/istl.hh>

// dune-vtk includes
#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>


// include stuff from dune-hydro
#include <dune/hydro/geotiffreader.hh>
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


      // create YaspGrid 
      double dx = std::abs(image.dLong());
      double dy = std::abs(image.dLat());
      double ox = image.originLong();
      double oy = image.originLat();
      std::array<int, 2> N;
      N[0] = 80;
      N[1] = 100;
      std::array<double, 2> H;
      H[0] = 90; 
      H[1] = 90; 
      Dune::FieldVector<double, 2> L;
      L[0] = N[0] * H[0];
      L[1] = N[1] * H[1];

      typedef Dune::YaspGrid<2> Grid;
      typedef Grid::ctype DF;
      typedef double RF;
      auto gridp = std::make_shared<Grid>(L, N, std::bitset<2>(0ULL), 1);

      // now make raster canvas in cell-centered mode 
      auto elevation_raster    = RasterDataSet<float>         (99.4 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); 
      auto accumulation_raster = RasterDataSet<float>         (99.4 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      auto direction_raster    = RasterDataSet<unsigned char> (99.4 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      
      elevation_raster.paste(image);
      accumulation_raster.paste(aimage);
      direction_raster.paste(dimage); 

      // create grid with original map  
      typedef Grid::LeafGridView GV;
      GV gv = gridp->leafGridView();

      using FEM = Dune::PDELab::P0LocalFiniteElementMap<RF, RF, 2>;
      using CON = Dune::PDELab::P0ParallelConstraints;
      using VBE = Dune::PDELab::ISTL::VectorBackend<>;
      using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
      using Z = Dune::PDELab::Backend::Vector<GFS, RF>;
      using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS, Z>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

      FEM fem(Dune::GeometryTypes::cube(2));
      CON con;
      GFS gfs(gv, fem, con);
      gfs.name("Vh");
      Z z(gfs);
      ZDGF zdgf(gfs, z);
      Z az(gfs);
      ZDGF azdgf(gfs, az);
      Z dz(gfs);
      ZDGF dzdgf(gfs, dz);

      auto bathymmetrylambda = [&](const auto &e, const auto &xlocal){
        auto x = e.geometry().center();
        int cellx = std::floor(x[0] / H[0]);
        int celly = std::floor(x[1] / H[1]);
        return elevation_raster(cellx, celly);
      };
      auto bathymmetrygf = Dune::PDELab::makeGridFunctionFromCallable(gv, bathymmetrylambda);
      Dune::PDELab::interpolate(bathymmetrygf, gfs, z);

      auto accumulationlambda = [&](const auto &e, const auto &xlocal){
        auto x = e.geometry().center();
        int cellx = std::floor(x[0] / H[0]);
        int celly = std::floor(x[1] / H[1]);
        return accumulation_raster(cellx, celly);
      };
      auto accumulationgf = Dune::PDELab::makeGridFunctionFromCallable(gv, accumulationlambda);

      auto directionlambda = [&](const auto &e, const auto &xlocal){
        auto x = e.geometry().center();
        int cellx = std::floor(x[0] / H[0]);
        int celly = std::floor(x[1] / H[1]);
        int value = static_cast<int>(direction_raster(cellx, celly));
        return static_cast<float>(value);
      };
      auto directiongf = Dune::PDELab::makeGridFunctionFromCallable(gv, directionlambda);


      Dune::PDELab::interpolate(bathymmetrygf, gfs, z);
      Dune::PDELab::interpolate(accumulationgf, gfs, az);
      Dune::PDELab::interpolate(directiongf, gfs, dz);
      using Writer = Dune::VtkImageDataWriter<GV>;
      Dune::PvdWriter<Writer> pvdWriter(gv, Dune::Vtk::FormatTypes::COMPRESSED, Dune::Vtk::DataTypes::FLOAT32);
      pvdWriter.addCellData(std::make_shared<VTKF>(zdgf, "bathymmetry"));
      pvdWriter.addCellData(std::make_shared<VTKF>(azdgf, "accumulation"));
      pvdWriter.addCellData(std::make_shared<VTKF>(dzdgf, "direction"));
      std::string fullfilename = "rivers.vti";
      pvdWriter.writeTimestep(0.0, fullfilename);

   
      
      // create grid for refinement
      typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid_;
      using GridView = Grid_::LeafGridView;

      const std::array<unsigned, 2> n = {N[0], N[1]};
      const Dune::FieldVector<double, 2> lower = {0, 0};
      const Dune::FieldVector<double, 2> upper = {L[0], L[1]};


      std::shared_ptr<Grid_> grid = Dune::StructuredGridFactory<Grid_>::createSimplexGrid(lower, upper, n);
      const GridView gridView = grid->leafGridView();

      std::vector<std::vector<flowFragment>> fragments = detectFragments(accumulation_raster, direction_raster, H, N, 50, 200);

      std::cout << "Start refining (might take a while)" << std::endl;
      refineGridwithFragments(grid, fragments, N, H,50,  0.5, false  );

      std::vector<double> height = overallHeight(gridView, elevation_raster, H, N);
      height= applyFlowHeightFragments(gridView, fragments, height, elevation_raster, H, N);

      Dune::VTKWriter<GridView> vtkWriter(gridView);
      vtkWriter.addVertexData(height, "height");
      
      vtkWriter.write("mapWithRefinedRivers02");
    }

    // close GDAL
    GDALDestroyDriverManager();
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
