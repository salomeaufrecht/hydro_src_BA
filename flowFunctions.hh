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
#include <tuple>

#include <array>


/**
 * @brief Stores basic data for each flow fragment, including start and end points, 
 *        as well as optional parameters for width and depth.
 * 
 * The struct defines a flow fragment in terms of its geometric properties. 
 */
struct flowFragment
{
    Dune::FieldVector<double, 2> start; ///< Starting point of the flow fragment.
    Dune::FieldVector<double, 2> end;   ///< Endpoint of the flow fragment.
    double widthStart;                    ///< Width of the flow at the starting point.
    double widthEnd = -1;                 ///< Width of the flow at the endpoint (default: -1 means same as start).
    double depht = -2;                    ///< Depth of the flow (default: -2 means no value).
    int maxIterationsFragment = 1000;     ///< Maximum number of iterations for refinement.
    int direction = 0;                     ///< Direction of the fragment (0: any, 1: horizontal, 2: vertical, 3: diagonal).	
    double minSize = -1;               ///< Minimum size of the fragment for refinement purposes.
    bool operator<(const flowFragment& f) const{
        return start[0] < f.start[0] || (start[0] == f.start[0] && start[1] < f.start[1]);
    }
    bool operator==(const flowFragment& f) const{
        return start[0] == f.start[0] && start[1] == f.start[1] && end[0] == f.end[0] && end[1] == f.end[1] && widthStart == f.widthStart && widthEnd == f.widthEnd && depht == f.depht;
    }
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
    Dune::FieldVector<double, 2> normal;  ///< The vector pointing away from the flow vector.
    int direction = 0;                     ///< Direction of the fragment (0: any, 1: horizontal, 2: vertical, 3: diagonal).
    double minSize = 0.1;                  ///< Minimum size of the fragment for refinement purposes.
    double depht = 1;                      ///< Depth of the flow fragment.
    int maxIterationsFragment = 50;               ///< Maximum number of iterations for refinement.
};


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
 * @param minAcc The minimum accumulation threshold for detecting fragments (default: 50).
 * @param maxAccDiff The maximum allowed difference in accumulation within a fragment (default: 200).
 * @param fixedWidth Whether to use a fixed width for all fragments (default: true).
 * @param scaleDephtFactor The scaling factor to calculate flow depth from accumalation values(default: 600).
 * @param scaleWidthFactor The scaling factor to calculate flow width from accumalation values(default: 90).
 * @param minWidth The minimum allowed width of a fragment (default: 1.0).
 * @param maxWidth The maximum allowed width of a fragment (default: -1 which means cellSize).
 * @return std::vector<flowFragment> A vector of detected flow fragments grouped in bounding boxes.
 */
std::vector<std::vector<flowFragment>> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, 
                            std::array<double, 2> cellSize, std::array<int, 2> gridSize, double minAcc=50, double maxAccDiff=200, bool fixedWidth=true, 
                            double scaleDephtFactor=600, double scaleWidthFactor=90, double minWidth=1.0, double maxWidth=-1);

/**
 * @brief Refines the grid to ensure accurate representation of flow fragments.
 * 
 * Adjusts the grid's resolution to properly capture the geometric details of the flow 
 * fragments. Refinement stops when the triangle size is smaller than the specified 
 * minimum size relative to the fragment width.
 * 
 * @param grid The grid where the refinement will be applied.
 * @param fragments A vector containing all flow fragments to be refined.
 * @param gridSize The dimensions of the grid.
 * @param cellSize The size of each grid cell (default: {90.0, 90.0}).
 * @param maxIterations The maximal amount of refinement iterations (default: 50).
 * @param minSizeFactor The minimum size factor for refinement, relative to fragment width (default: 0.4).
 * @param exactCalc Whether to use exact calculations for the fragment width (default: false).
 */
void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, 
        std::vector<std::vector<flowFragment>> fragments,std::array<int, 2> gridSize, std::array<double, 2> cellSize = {90.0, 90.0}, int maxIterations=50, 
        double minSizeFactor=0.4, bool exactCalc = false);

/**
 * @brief Applies height values to all vertices in the grid.
 * 
 * Modifies the height vector for the grid vertices by incorporating elevation data from 
 * a raster dataset.
 * 
 * @param gridView The gridview of the ALUGrid representing the map.
 * @param elevation_raster The elevation raster dataset containing height values.
 * @param cellSize The size of each cell in the grid.
 * @param gridSize The dimensions of the grid.
 * @return std::vector<double> The updated height values.
 */
std::vector<double> overallHeight(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
                                RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, 2> gridSize);

/**
 * @brief Applies depth to all river fragments by setting the corresponding height values.
 * 
 * Updates the height values for vertices in the grid based on the depth of the flow fragments. 
 * The height value for each affected vertex is substracted by `depht`.
 * 
 * @param grid The gridview of the ALUGrid representing the map.
 * @param fragments A vector containing all flow fragments grouped by bounding box.
 * @param height A vector of height values for each vertex, modified in this function.
 * @param elevation_raster The elevation raster dataset for reference.
 * @param cellSize The size of each cell in the grid.
 * @param gridSize The dimensions of the grid.
 * @return std::vector<double> The updated height values.
 */
std::vector<double> applyFlowHeightFragments(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
       std::vector<std::vector<flowFragment>> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, 
        std::array<int, 2> gridSize);

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
 * @param minAcc The minimum accumulation value for flow detection (default: 50).
 * @param maxAccDiff The maximum accumulation difference for flow detection (default: 200).
 * @param fixedWidth Whether to use a fixed width for all fragments (default: true).
 * @param scaleDephtFactor The scaling factor to calculate flow depth from accumulation values (default: 400).
 * @param scaleWidthFactor The scaling factor to calculate flow width from accumulation values (default: 5000).
 * @param minSizeFactor The minimum size factor for refinement (default: 0.4).
 * @param maxIterations The maximal amount of refinement iterations (default: 50).
 * @param minWidth The minimum width of a flow fragment in world coordinates (default: 1.0).
 * @param exactCalc Whether to use exact calculations for the fragment width (default: false).
 * @return std::vector<double> 
 */
std::vector<double> addRiversToMap(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid,
                                    std::array<double, 2> cellSize, std::array<int, 2> gridSize,
                                    RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                                    double minAcc = 50, double maxAccDiff = 200, bool fixedWidth = true, double scaleDephtFactor = 600, 
                                    double scaleWidthFactor = 90, double minSizeFactor = 0.4, int maxIterations = 50, double minWidth = 1.0, double maxWidth = -1,
                                    bool exactCalc = false);

#endif // FLOWFUNCTIONS_HH
