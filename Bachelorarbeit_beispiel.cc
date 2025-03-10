#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>

// dune-common includes
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/hydro/geotiffreader.hh>

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

        // read data files
        GeoTIFFImage<GInt16> image("n00e090_con.tif", "", 1, 2);
        GeoTIFFImage<GUInt32> aimage("N00E090_acc.tif", "", 1, 2);
        GeoTIFFImage<GByte> dimage("n00e090_dir.tif", "", 1, 2);

        double dx = std::abs(image.dLong());
        double dy = std::abs(image.dLat());

        //size of the map
        std::array<int, 2> N;
        N[0] = 30;
        N[1] = 30;
        std::array<double, 2> H;
        H[0] = 90; 
        H[1] = 90; 
        Dune::FieldVector<double, 2> L;
        L[0] = N[0] * H[0];
        L[1] = N[1] * H[1];

        //paste data into raster
        auto elevation_raster    = RasterDataSet<float>(99.21 + 0.5 * dx, 8.4 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
        auto accumulation_raster = RasterDataSet<float>(99.21 + 0.5 * dx, 8.4 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
        auto direction_raster    = RasterDataSet<unsigned char>(99.21 + 0.5 * dx, 8.4 + 0.5 * dy, dx, dy, N[0], N[1], 0, 1);
        elevation_raster.paste(image);
        accumulation_raster.paste(aimage);
        direction_raster.paste(dimage);

        // create grid
        typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
        using GridView = Grid::LeafGridView;

        const std::array<unsigned, 2> n = {N[0], N[1]};
        const Dune::FieldVector<double, 2> lower = {0, 0};
        const Dune::FieldVector<double, 2> upper = {L[0] , L[1] };

        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);
        GridView gridView = grid->leafGridView();

        //compute rivers and add them to the map
        std::vector<std::vector<flowFragment>> rivers = detectFragments(accumulation_raster, direction_raster, H, N);
        refineGridwithFragments(grid, rivers, N, H);
        std::vector<double> height = overallHeight(gridView, elevation_raster, H, N);
        height= applyFlowHeightFragments(gridView, rivers, height, elevation_raster, H, N);
       
        // or use one function to do all of the above
        //std::vector<double> height = addRiversToMap(grid, H, N, accumulation_raster, direction_raster, elevation_raster);


        // write result to file
        Dune::VTKWriter<GridView> vtkWriter(gridView);
        vtkWriter.addVertexData(height, "height");
        vtkWriter.write("mapWithRefinedRivers01");

    
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
