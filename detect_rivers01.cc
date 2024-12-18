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
    std::cout << "Hello World! This is dune-hydro." << std::endl;
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

      // create YaspGrid with 1800 x 1200 cells
      double dx = std::abs(image.dLong());
      double dy = std::abs(image.dLat());
      double ox = image.originLong();
      double oy = image.originLat();
      const int dim = 2;
      std::array<int, dim> N;
      N[0] = 60; //1800
      N[1] = 60; //1200
      std::array<double, dim> H;
      H[0] = 1; //90
      H[1] = 1; //90
      Dune::FieldVector<double, dim> L;
      L[0] = N[0] * H[0];
      L[1] = N[1] * H[1];

    

      typedef Dune::YaspGrid<dim> Grid;
      typedef Grid::ctype DF;
      typedef double RF;
      auto gridp = std::make_shared<Grid>(L, N, std::bitset<dim>(0ULL), 1);

      // now make raster canvas in cell-centered mode
      //auto elevation_raster = RasterDataSet<float>(500, 100, dx, dy, N[0], N[1], 0, 1);
      auto elevation_raster = RasterDataSet<float>(99.0 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);

      elevation_raster.paste(image);
      //elevation.paste(image2);
      auto accumulation_raster = RasterDataSet<float>(99.0 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
      //auto accumulation_raster = RasterDataSet<float>(500, 100, dx, dy, N[0], N[1], 0, 1);
      accumulation_raster.paste(aimage);
      //accumulation_raster.paste(aimage2);
      //auto direction_raster = RasterDataSet<unsigned char>(500, 100, dx, dy, N[0], N[1], 0, 1);
      auto direction_raster = RasterDataSet<unsigned char>(99.0 + 0.5 * dx, 8.0 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
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




  std::vector<flowFragment> rivers;
  int size = int(N[0]*N[1]);
  std::vector<int> skip(size, 0);

      for (int i = 1; i< N[0]-1; i++){ //TODO: start from 0 and check if end is outside
        for (int j=1; j<N[1]-1; j++){
          if(skip[j*N[0]+i]==-1) {std::cout << "."; continue;}
          if(accumulation_raster(i, j)>50 && accumulation_raster(i, j)< 4294967295)  {

            Dune::FieldVector<double, 2> start = {i, j};
             Dune::FieldVector<double, 2> end = start;
             

             int dir = int(direction_raster(i, j));
             //std::cout <<"corr dir: " << dir << ", end: " << end << std::endl;

             int endI=round(end[0]);
             int endJ = round(end[1]);
             int endDir = int(direction_raster(endI, endJ));

             double accStart = accumulation_raster(i, j);
             double accEnd = accStart;
             

            while (endDir == dir){
              
              //std::cout << " enddir: " << endDir << ", dir:" << dir << ", end:" << end << std::endl;
              skip[endJ*N[0]+endI]=-1;
              
              Dune::FieldVector<double, 2> currPoint = end;
              switch (dir)
              {
              case 1: end= currPoint+ Dune::FieldVector<double, 2>{1, 0};
                break;
              case 2: end= currPoint+ Dune::FieldVector<double, 2>{1,-1};
                break;
                case 4: end= currPoint+ Dune::FieldVector<double, 2>{0, -1};
                break;
                case 8: end= currPoint+ Dune::FieldVector<double, 2>{-1, -1};
                break;
                case 16: end= currPoint+ Dune::FieldVector<double, 2>{-1, 0};
                break;
                case 32: end= currPoint+ Dune::FieldVector<double, 2>{-1, 1};
                break;
                case 64: end= currPoint+ Dune::FieldVector<double, 2>{0, 1};
                break;
                case 128: end= currPoint+ Dune::FieldVector<double, 2>{1, 1};
                break;
                case 255: end= currPoint+ Dune::FieldVector<double, 2>{0, 0};
                break;
              default:
                break;
              }

              endI = round(end[0]); 
              endJ = round(end[1]);
              endDir = int(direction_raster(endI, endJ));
              accEnd = accumulation_raster(endI, endJ);

              if(end[0] > N[0]-1 || end[1]>N[1]-1 || end[0] <1 || end[1]<1 || std::abs(accStart-accEnd)>200) {
                if(currPoint==start &&  std::abs(accStart-accEnd)>100){break;}
                end=currPoint;
                endI = round(end[0]); 
                endJ = round(end[1]);
                skip[endJ*N[0]+endI]=0; //end should not be skipped
                accEnd = accumulation_raster(endI, endJ); 
                break;
                }

              //std::cout << "new enddir: " << endDir<< ", new end: " << end << std::endl;

            }

            skip[j*N[0]+i] =0; //start should not be skipped
            double volume = (accStart+accEnd)/2;
            double width = std::sqrt(volume/3000) * 3; //width:depht 3:1
            double depht = std::sqrt(volume/3000);
            width = std::min(0.6, width);
            depht = std::min(0.4, depht);

            start[0] = start[0]*H[0];
            start[1] = start[1]*H[1];
            end[0] = end[0]*H[0];
            end[1] = end[1]*H[1];

            flowFragment f = {start, end, width*H[1]};
            f.depht = depht;
            //std::cout << "depht: " << depht << ", width: " << width << ", start: " << start << ", end: " << end << std::endl;

            rivers.push_back(f);
        }
      }
      }

    


        
    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid_;
    using GridView = Grid_::LeafGridView;


    // Start with a structured grid
    const std::array<unsigned, dim> n_ = {N[0], N[1]};
    const Dune::FieldVector<double, dim> lower = {0, 0};
    const Dune::FieldVector<double, dim> upper = {L[0], L[1]};


  std::cout << "make grid" << std::endl;
    std::shared_ptr<Grid_> grid = Dune::StructuredGridFactory<Grid_>::createSimplexGrid(lower, upper, n_);
    //grid->globalRefine(1); for conforming
  std::cout << "grid done" << std::endl;
  std::cout << "make gridView" << std::endl;

    const GridView gridView = grid->leafGridView();
    std::cout << "gridview done" << std::endl;

  std::vector<flowFragment> fragments;

  fragments = rivers;

    for(auto& f : fragments){
     int endI=round(f.start[0]/H[0]);
      int endJ = round(f.start[1]/H[1]);
      if(skip[endJ*N[0]+endI]==-1) continue;
      flowWithFragments(grid, f, 0.5, H[0]);
    }

    std::vector<double> height(gridView.indexSet().size(2));
    // height = initializeHeight(height, elevation_raster, N[0], N[1]);

  int c = 0;
    for(auto& f : fragments){
     int endI=round(f.start[0]/H[0]);
      int endJ = round(f.start[1]/H[1]);
      if(skip[endJ*N[0]+endI]==-1){std::cout << "."; continue;}
      c++;
      height = applyFlowHeightFragments(grid, f, height);
    }
    std::cout << c << std::endl;

    height = overallHeigth(grid, height, elevation_raster);
     
    //for(auto& f : fragments){
    //  int endI=round(f.start[0]/H[0]);
    //  int endJ = round(f.start[1]/H[1]);
    //  if(skip[endJ*N[0]+endI]==-1) continue;
    //  //std::cout << "height" << std::endl;
    //  height = adjustFlowHeightFragments(grid, f, height); //TODO: make more efficient, better
    //}
     
     //Write grid to file

    int subsampling = 2;
    //Dune::SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementIntervals(subsampling));
    //auto f = std::make_shared<MyVTKFunction<GridView>>(gridView,height,"height");
    //vtkWriter.addCellData(f);
    Dune::VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.addVertexData(height, "height");
   
    vtkWriter.write("nakthon_rivers");


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
