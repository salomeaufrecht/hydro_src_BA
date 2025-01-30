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
        GeoTIFFImage<GInt16> image("n00e090_con.tif", "", 1, 2);
        GeoTIFFImage<GUInt32> aimage("N00E090_acc.tif", "", 1, 2);
        GeoTIFFImage<GByte> dimage("n00e090_dir.tif", "", 1, 2);
    
        double dx = std::abs(image.dLong());
        double dy = std::abs(image.dLat());
        const int dim = 2;
        std::array<double, dim> H;
        H[0] = 1; 
        H[1] = 1;
        Dune::FieldVector<double, dim> L;

        std::array<int, dim> N;
        std::vector<std::array<double, 2>> runtime;
        for(int i = 5; i < 1001; i+=5){
            
            N[0] = i;
            N[1] = i;
            std::cout << i << std::endl;

            //auto elevation_raster = RasterDataSet<float>(99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1); //original: 99.0, 8.0 (99 breite, 8 höhe)
            //elevation_raster.paste(image);
            auto accumulation_raster = RasterDataSet<float>(99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
            accumulation_raster.paste(aimage);
            auto direction_raster = RasterDataSet<unsigned char>(99.405 + 0.5 * dx, 8.11 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
            direction_raster.paste(dimage);

            Dune::Timer timer;  // Timer starten
            std::vector<flowFragment> rivers = detectFragments(accumulation_raster, direction_raster, H, N);
            double elapsedTime = timer.stop();  // Stoppe den Timer und bekomme die Zeit in Sekunden
            runtime.push_back({i*i, elapsedTime});
        }

        std::cout << std::endl << "x1 = [ " << runtime[0][0] ;
        for (auto r : runtime){
            std::cout << ", " << r[0];
        }
        std::cout << "]" << std::endl << "y1 = [" << runtime[0][1] ;
        for (auto r : runtime){
            std::cout << ", "<< r[1] ;
        }
        std::cout << "]" << std::endl;
        std::cout << std::endl << "plt.plot(x1, y1, \"-\") \nplt.ylabel(\"sec\") \nplt.xlabel(\"size of map\")" << std::endl;
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
