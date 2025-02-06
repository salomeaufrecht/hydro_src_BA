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


/**
 * @brief Stores basic data for each flow fragment, including start and end points, 
 *        as well as optional parameters for width and depth.
 * 
 * The struct defines a flow fragment in terms of its geometric properties. 
 * `widthEnd` and `depht` are optional parameters with default values, used to 
 * define the fragment's width and depth where applicable.
 */
struct flowFragment
{
    Dune::FieldVector<double, 2> start; ///< Starting point of the flow fragment.
    Dune::FieldVector<double, 2> end;   ///< Endpoint of the flow fragment.
    double widthStart;                    ///< Width of the flow at the starting point.
    double widthEnd = -1;                 ///< Width of the flow at the endpoint (default: -1 means same as start).
    double depht = -2;                    ///< Depth of the flow (default: -2 means no value).
};

/**
 * @brief Stores detailed boundary information for a flow fragment to optimize computations.
 * 
 * The struct holds precomputed details about a flow fragment's boundaries, including 
 * its bounding box, direction, and additional derived properties. This avoids redundant 
 * recomputation of these properties during processing.
 */
struct fragmentBoundaries{
    Dune::FieldVector<double, 2> start1; ///< One endpoint of the flow fragment's first boundary.
    Dune::FieldVector<double, 2> start2; ///< One endpoint of the flow fragment's second boundary.
    Dune::FieldVector<double, 2> b1;     ///< Vector defining the first boundary.
    Dune::FieldVector<double, 2> b2;     ///< Vector defining the second boundary.
    double minX;                           ///< Minimum x-coordinate of the bounding box.
    double maxX;                           ///< Maximum x-coordinate of the bounding box.
    double minY;                           ///< Minimum y-coordinate of the bounding box.
    double maxY;                           ///< Maximum y-coordinate of the bounding box.
    int direction = 0;                     ///< Direction of the fragment (0: any, 1: horizontal, 2: vertical, 3: diagonal).
    double minSize = 0.1;                  ///< Minimum size of the fragment for refinement purposes.
    double depht = 0;                      ///< Depth of the flow fragment.
    Dune::FieldVector<double, 2> normal;  ///< The vector pointing away from the flow vector.
};

/**
 * @brief Applies depth to all river fragments by setting the corresponding height values.
 * 
 * Updates the height values for vertices in the grid based on the depth of the flow fragments. 
 * The height value for each affected vertex is substracted by `depht`.
 * 
 * @param grid The ALUGrid representing the map.
 * @param fragments A vector containing all flow fragments.
 * @param height A vector of height values for each vertex, modified in this function.
 * @param elevation_raster The elevation raster dataset for reference.
 * @param cellSize The size of each cell in the grid.
 * @param gridSize The dimensions of the grid.
 * @return std::vector<double> The updated height values.
 */
std::vector<double> applyFlowHeightFragments(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
        std::vector<flowFragment> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, 2> gridSize);

/**
 * @brief Applies height values to all vertices in the grid.
 * 
 * Modifies the height vector for the grid vertices by incorporating elevation data from 
 * a raster dataset.
 * 
 * @param grid The ALUGrid representing the map.
 * @param height A vector of height values for each vertex (updated in this function).
 * @param elevation_raster The elevation raster dataset containing height values.
 * @param cellSize The size of each cell in the grid.
 * @return std::vector<double> The updated height values.
 */
std::vector<double> overallHeight(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
                                RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, 2> gridSize);

/**
 * @brief Refines the grid to ensure accurate representation of flow fragments.
 * 
 * Adjusts the grid's resolution to properly capture the geometric details of the flow 
 * fragments. Refinement stops when the triangle size is smaller than the specified 
 * minimum size relative to the fragment width.
 * 
 * @param grid The grid where the refinement will be applied.
 * @param fragments A vector containing all flow fragments to be refined.
 * @param minSizeFactor The minimum size factor for refinement, relative to fragment width (default: 0.4).
 * @param cellSize The size of each grid cell (default: {1.0, 1.0}).
 * @param maxIterations The maximal amount of refinement iterations (default: 50).
 * @param minMinSize The minimum for minSize in world coordinates (default: 0.2).
 */
void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor = 0.4, std::array<double, 2> cellSize = {1.0, 1.0}, int maxIterations=50, double minMinSize=0.2);

/**
 * @brief Detects flow fragments based on accumulation and direction raster data.
 * 
 * Identifies flow fragments by analyzing accumulation and direction data from raster datasets. 
 * The results include start and end points, as well as additional geometric details for 
 * each fragment.
 * 
 * @param accumulation_raster The raster dataset containing accumulation values.
 * @param direction_raster The raster dataset containing flow direction values.
 * @param cellSize The size of each grid cell in physical units.
 * @param gridSize The dimensions of the grid.
 * @param minAcc The minimum accumulation value for flow detection (default: 50).
 * @param maxAccDiff The maximum accumulation difference for flow detection (default: 200).
 * @param scaleDephtFactor The scaling factor to calculate flow depth from accumalation values(default: 5000).
 * @param scaleWidthFactor The scaling factor to calculate flow width from accumalation values(default: 500).
 * @param minWidth The minimum width of a flow fragment in worldcoordinates(default: 1.0).
 * @return std::vector<flowFragment> A vector of detected flow fragments.
 */
std::vector<flowFragment> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster,
                                             std::array<double, 2> cellSize, std::array<int, 2> gridSize, double minAcc = 50, double maxAccDiff=200,
                                             double scaleDephtFactor=3000, double scaleWidthFactor=500, double minWidth = 1.0);


RasterDataSet<float> removeUpwardsRivers(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                            std::array<int, 2> gridSize, double minAcc = 50);

/**
 * @brief Detects rivers in an map, refines the grid to represent them accurately and adjusts the height values.
 * 
 * This function combines the detection of flow fragments, grid refinement, and height adjustment.
 * 
 * @param grid The grid where the refinement will be applied.
 * @param cellSize The size of each grid cell.
 * @param gridSize The dimensions of the grid.
 * @param accumulation_raster The raster dataset containing accumulation values.
 * @param direction_raster The raster dataset containing flow direction values.
 * @param elevation_raster The raster dataset containing elevation values.
 * @param minSizeFactor The minimum size factor for refinement (default: 0.4).
 * @param minAcc The minimum accumulation value for flow detection (default: 50).
 * @param maxAccDiff The maximum accumulation difference for flow detection (default: 200).
 * @param scaleDephtFactor The scaling factor to calculate flow depth from accumulation values (default: 400).
 * @param scaleWidthFactor The scaling factor to calculate flow width from accumulation values (default: 5000).
 * @param maxIterations The maximal amount of refinement iterations (default: 50).
 * @param minWidth The minimum width of a flow fragment in world coordinates (default: 1.0).
 * @param minMinSize The minimum for minSize in world coordinates (default: 0.2).
 * @return std::vector<double> 
 */
std::vector<double> addRiversToMap(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, std::array<double, 2> cellSize, std::array<int, 2> gridSize,
                                    RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                                    double minSizeFactor = 0.4, double minAcc = 50, double maxAccDiff = 200, double scaleDephtFactor = 400, 
                                    double scaleWidthFactor = 5000, int maxIterations = 50, double minWidth = 1.0, double minMinSize = 0.2);


/**
 * @brief Calculates the real cell size based on the grid size and cell size.
 * 
 * @param cellSize The selected size of each grid cell.
 * @param gridSize The dimensions of the grid.
 * @return std::array<double, 2> The actual size of each grid cell.
 */
std::array<double, 2> calcRealCellSize(std::array<double, 2> cellSize, std::array<int, 2> gridSize);

#endif // FLOWFUNCTIONS_HH
