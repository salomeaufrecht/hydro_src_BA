#ifndef FLOWFUNCTIONS_HH
#define FLOWFUNCTIONS_HH


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
#include <dune/hydro/rasterdataset.hh>

#include <array>


constexpr int dim = 2; // Dimension des Grids

/**
 * @brief stores data for each fragment (widthEnd and depht being optional)
 * 
 */
struct flowFragment
{
    Dune::FieldVector<double, dim> start;
    Dune::FieldVector<double, dim> end;
    double widthStart;
    double widthEnd = -1;
    double depht = -2;
};

/**
 * @brief for each fragment stores more detailed information than flowFragment.
 * All information can be computed from flowFragment but is stored here to avoid recomputation
 * 
 */
struct fragmentBoundaries{
    Dune::FieldVector<double, dim> start1;
    Dune::FieldVector<double, dim> start2;
    Dune::FieldVector<double, dim> b1;
    Dune::FieldVector<double, dim> b2;
    double minX;
    double maxX;
    double minY;
    double maxY;
    int direction = 0; //0: any; 1: horizontal; 2: vertikal; 3: diagonal
    double minSize = 0.1;
    double depht = 0;
};


/**
 * @brief applies height for all river fragments by setting corresponding height value to -depht
 * 
 * @param grid with refined areas where rivers are
 * @param fragments vector with all flow fragments
 * @param height vector with entry for each vertex where the corresponding entry will be set to -depht
 * @return std::vector<double> height
 */
std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize);


/**
     * @brief applies height to all vertices in grid
     * @param grid alugrid representing an map
     * @param height vector with entry for each vertex where the height value will be added to and stored (changes of height values should be done in here (e.g. -1 for rivers))
     * @param map rasterdataset with height values for every cell in original data
     * @param cellSize size of each cell
     */
std::vector<double> overallHeight(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize);


/**
 * @brief given a list of flow fragmants it refines the grid to a level where the flow can be portrailed properly
 * 
 * @param grid grid in which rivers should be inserted
 * @param fragments vector with all flow fragments that should be inserted
 * @param minSizeFactor the minimal size to which a triangle should be refined relative to the width of each fragment 
 * (e.g. f.width=2, minSizeFactor=0.5 --> stop refinement when triangles are smaller than 1 )
 * @param cellSize pixelSize
 */
void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor = 0.4, std::array<double, 2> cellSize = {1.0, 1.0});


/**
 * @brief detects rivers by looking at high accumulation values. Rivers are splitted in straight fragments
 * 
 * @param accumulation_raster raster with accumulation value for every cell
 * @param direction_raster raster with direction value for every cell
 * @param cellSize size of each cell in original grid
 * @param gridSize size of grid (amount of cells)
 * @return std::vector<flowFragment> 
 */
std::vector<flowFragment> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize);

std::vector<double> adjustFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f,
                                             std::vector<double> height);

std::vector<double> adjustFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                    Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, 
                                    double depht, std::vector<double> height);

#endif // FLOWFUNCTIONS_HH
