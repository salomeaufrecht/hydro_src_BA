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


constexpr int dim = 2; // dimension of grid


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
    Dune::FieldVector<double, dim> start; ///< Starting point of the flow fragment.
    Dune::FieldVector<double, dim> end;   ///< Endpoint of the flow fragment.
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
    Dune::FieldVector<double, dim> start1; ///< One endpoint of the flow fragment's first boundary.
    Dune::FieldVector<double, dim> start2; ///< One endpoint of the flow fragment's second boundary.
    Dune::FieldVector<double, dim> b1;     ///< Vector defining the first boundary.
    Dune::FieldVector<double, dim> b2;     ///< Vector defining the second boundary.
    double minX;                           ///< Minimum x-coordinate of the bounding box.
    double maxX;                           ///< Maximum x-coordinate of the bounding box.
    double minY;                           ///< Minimum y-coordinate of the bounding box.
    double maxY;                           ///< Maximum y-coordinate of the bounding box.
    int direction = 0;                     ///< Direction of the fragment (0: any, 1: horizontal, 2: vertical, 3: diagonal).
    double minSize = 0.1;                  ///< Minimum size of the fragment for refinement purposes.
    double depht = 0;                      ///< Depth of the flow fragment.
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
std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize);

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
std::vector<double> overallHeight(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize);

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
 */
void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor = 0.4, std::array<double, 2> cellSize = {1.0, 1.0});

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
 * @return std::vector<flowFragment> A vector of detected flow fragments.
 */
std::vector<flowFragment> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize);

/**
 * @brief Adjusts height values for a specific flow fragment.
 * 
 * Modifies the height values of the vertices influenced by a given flow fragment, 
 * considering its depth and geometry.
 * 
 * @param grid The ALUGrid representing the map.
 * @param f The flow fragment whose height values are to be adjusted.
 * @param height A vector of height values for the grid vertices.
 * @return std::vector<double> The updated height values.
 */
std::vector<double> adjustFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f,
                                             std::vector<double> height);

/**
 * @brief Adjusts height values between two points based on flow parameters.
 * 
 * Modifies the height values of the vertices along a line defined by the 
 * start and end points, taking into account the width, depth, and geometry of the flow.
 * 
 * @param grid The ALUGrid representing the map.
 * @param start The starting point of the flow.
 * @param end The endpoint of the flow.
 * @param widthStart The width of the flow at the starting point.
 * @param widthEnd The width of the flow at the endpoint.
 * @param depht The depth of the flow.
 * @param height A vector of height values for the grid vertices.
 * @return std::vector<double> The updated height values.
 */
std::vector<double> adjustFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                    Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, 
                                    double depht, std::vector<double> height);

#endif // FLOWFUNCTIONS_HH
