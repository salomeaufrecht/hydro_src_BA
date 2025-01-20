// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*************************************/
/* This example solves surface flow in the vicinity of Heidelberg
   with 30m resolution and constant rainfall. Uses structured mesh
   and finite volume method.
 */
/*************************************/

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

template<typename GV>
class MyVTKFunction
    : public Dune::VTKFunction<GV>
{
  typedef typename GV::Grid::ctype DF;
  enum
  {
    n = GV::dimension
  };
  using Entity = typename GV::Grid::template Codim<0>::Entity;

public:
  MyVTKFunction(const GV& gv_, std::vector<double>& data_, std::string s_)
      : gv(gv_), data(data_), s(s_)
  {
  }

  virtual int ncomps() const override
  {
    return 1;
  }

  virtual double evaluate(int comp, const Entity &e, const Dune::FieldVector<DF, n> &local) const override
  {
    auto x = e.geometry().global(local);
    double minnorm = 1e10;
    int mincorner;
    for (int i=0; i<e.geometry().corners(); i++)
    {
      auto xc = e.geometry().corner(i);
      xc -= x;
      auto norm = xc.two_norm();
      if (norm<minnorm) { minnorm=norm; mincorner=i; }
    }
    return data[gv.indexSet().subIndex(e,mincorner,2)];
  }

  virtual std::string name() const override
  {
    return s;
  }

private:
  const GV& gv;
  std::vector<double>& data;
  std::string s;
};




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
      // first tile
      GeoTIFFImage<GInt16> image("n00e090_con.tif", "", 1, 2);
      GeoTIFFImage<GUInt32> aimage("N00E090_acc.tif", "", 1, 2);
      GeoTIFFImage<GByte> dimage("n00e090_dir.tif", "", 1, 2);
      // second tile
      //GeoTIFFImage<GInt16> image2("n00e100_con.tif", "", 1, 2);
      //GeoTIFFImage<GUInt32> aimage2("N00E100_acc.tif", "", 1, 2);
      //GeoTIFFImage<GByte> dimage2("n00e100_dir.tif", "", 1, 2);

      // create YaspGrid 
      double dx = std::abs(image.dLong());
      double dy = std::abs(image.dLat());
      double ox = image.originLong();
      double oy = image.originLat();
      const int dim = 2;
      std::array<int, dim> N;
      N[0] = 16; //1800
      N[1] = 16; //1200
      std::array<double, dim> H;
      H[0] = 6; 
      H[1] = 6; 
      Dune::FieldVector<double, dim> L;
      L[0] = N[0] * H[0];
      L[1] = N[1] * H[1];

      typedef Dune::YaspGrid<dim> Grid;
      typedef Grid::ctype DF;
      typedef double RF;
      auto gridp = std::make_shared<Grid>(L, N, std::bitset<dim>(0ULL), 1);

      // now make raster canvas in cell-centered mode
      auto elevation_raster = RasterDataSet<float>(99.411 + 0.5 * dx, 8.1 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); //original: 99.0, 8.0 (99 breite, 8 höhe)

      elevation_raster.paste(image);
      //elevation.paste(image2);
      auto accumulation_raster = RasterDataSet<float>(99.411 + 0.5 * dx, 8.1 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      accumulation_raster.paste(aimage);
      //accumulation_raster.paste(aimage2);
      auto direction_raster = RasterDataSet<unsigned char>(99.411 + 0.5 * dx, 8.1 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      direction_raster.paste(dimage);
      //direction_raster.paste(dimage2);

      // write a grid file    
      typedef Grid::LeafGridView GV;
      GV gv = gridp->leafGridView();

      using FEM = Dune::PDELab::P0LocalFiniteElementMap<RF, RF, dim>;
      using CON = Dune::PDELab::P0ParallelConstraints;
      using VBE = Dune::PDELab::ISTL::VectorBackend<>;
      using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
      using Z = Dune::PDELab::Backend::Vector<GFS, RF>;
      using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS, Z>;
      using VTKF = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

      FEM fem(Dune::GeometryTypes::cube(dim));
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



  std::vector<flowFragment> rivers = detectFragments(accumulation_raster, direction_raster, H, N);


        
    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid_;
    using GridView = Grid_::LeafGridView;


    // create grid
    const std::array<unsigned, dim> n_ = {N[0], N[1]};
    const Dune::FieldVector<double, dim> lower = {0.5 * H[0], 0.5 * H[1]};
    const Dune::FieldVector<double, dim> upper = {L[0] - 0.5 * H[0], L[1] - 0.5 * H[1]};


  std::cout << "make grid" << std::endl;
  std::shared_ptr<Grid_> grid = Dune::StructuredGridFactory<Grid_>::createSimplexGrid(lower, upper, n_);
  std::cout << "grid done" << std::endl;

  std::cout << "make gridView" << std::endl;
  const GridView gridView = grid->leafGridView();
  std::cout << "gridview done" << std::endl;

  std::cout << "refine grid" << std::endl;
  refineGridwithFragments(grid, rivers, 0.5, H); 
  std::cout << "refine grid done " << std::endl;

  std::vector<double> height(gridView.indexSet().size(2), 0);
  std::cout << "set river height grid" << std::endl;
  height= applyFlowHeightFragments(grid, rivers, height);
  std::cout << "river heigth done" << std::endl;

  std::cout << "set overall height" << std::endl;
  height = overallHeight(grid, height, elevation_raster, H, N);
  std::cout << "overall heigth done" << std::endl;


    //for(auto& f : fragments){
    //  int endI=round(f.start[0]/H[0]);
    //  int endJ = round(f.start[1]/H[1]);
    //  if(skip[endJ*N[0]+endI]==-1) continue;
    //  //std::cout << "height" << std::endl;
    //  height = adjustFlowHeightFragments(grid, f, height); 
    //}
     
     //Write grid to file

    int subsampling = 2;
    //Dune::SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementIntervals(subsampling));
    //auto f = std::make_shared<MyVTKFunction<GridView>>(gridView,height,"height");
    //vtkWriter.addCellData(f);
    Dune::VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.addVertexData(height, "height");
   
    vtkWriter.write("nakhon_rivers");


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
