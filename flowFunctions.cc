#include "flowFunctions.hh"

/**
 * @brief Calculates the normal vector of a given 2D vector.
 * 
 * The function computes a unit vector normal to the input vector. 
 * The normal vector is rotated 90 degrees counter-clockwise.
 * 
 * @param vector The input vector for which the normal is calculated.
 * @return Dune::FieldVector<double, 2> The normalized normal vector.
 * @throws std::invalid_argument if the input vector has zero length.

 */
Dune::FieldVector<double, 2> calculateNormal(const Dune::FieldVector<double, 2>& vector) {
    double length = std::sqrt(vector[0] * vector[0] + vector[1] * vector[1]);
    if (length < 1e-8) {
        throw std::invalid_argument("Zero length vector provided");
    }
    return { -vector[1] / length, vector[0] / length };
}

/**
 * @brief Calculates the neighbors of a coordinate in a given grid.
 * 
 * The function determines the indices of the nearest grid neighbors (in raster data) for a given 
 * coordinate if this coordinate lies on the border between two grid cells in this coordinate direction. It ensures that the indices are within the bounds of the grid.
 * 
 * @param coord The coordinate whose neighbors are calculated (in raster coordinates).
 * @param gridSize The size of the grid.
 * @return std::pair<int, int> The indices of the two neighboring cells. 
 */
std::pair<int, int> calculateNeighborsOneDirection(double coord, int gridSize) {
    int neighbor1 = std::round(coord);
    int neighbor2 = -2;

    if (std::abs(coord - neighbor1) < 1e-8) {
        neighbor2 = int(coord);
        if (neighbor1 == neighbor2) neighbor2--;
    } else {
        neighbor1 = int(coord);
    }

    if (neighbor1 > gridSize - 1) neighbor1 = gridSize - 1;
    if (neighbor2 < 0 ) neighbor2 = -2;

    return {neighbor1, neighbor2};
}

/**
 * @brief Determines neighbors of a grid cell for a given point in raster coordinates.
 * 
 * The function identifies neighboring grid cells for a point in raster date.
 * It returns the indices of the adjacent cells in both x and y directions if the point lies on the 
 * border between grids. Otherwise, neighbor x/y 2 will be -2.
 * 
 * @param point The point in raster coordinates.
 * @param gridSize The dimensions of the grid (number of cells in x and y directions).
 * @return std::vector<int> The indices of the neighboring cells in x and y directions.
 */
std::vector<int> calcNeighbors(const Dune::FieldVector<double, 2>& point, std::array<int, 2> gridSize){

    auto [neighborX1, neighborX2] = calculateNeighborsOneDirection(point[0], gridSize[0]);
    auto [neighborY1, neighborY2] = calculateNeighborsOneDirection(point[1], gridSize[1]);

    return {neighborX1, neighborX2, neighborY1, neighborY2};
}

/**
 * @brief Calculates the squared distance between two points.
 * 
 * The function computes the Euclidean squared distance between two points.
 * 
 * @param p1 The first point.
 * @param p2 The second point.
 * @return double The squared distance between the two points.
 */
double squaredDistance(const Dune::FieldVector<double, 2>& p1, const Dune::FieldVector<double, 2>& p2){
    Dune::FieldVector<double, 2> diff = (p1 - p2);
    return diff[0]*diff[0] + diff[1]*diff[1];
}

/**
 * @brief Converts raster coordinates to global coordinates.
 * 
 * The function maps raster grid coordinates to global coordinates using 
 * specified scaling factors and offsets.
 * 
 * @param rasterCoord The raster coordinates.
 * @param resCellSize Scaling factors for each dimension.
 * @param half Offsets for each dimension.
 * @return Dune::FieldVector<double, 2> The corresponding global coordinates.
 */
Dune::FieldVector<double, 2> convertToGlobalCoordinates(const Dune::FieldVector<double, 2>& rasterCoord, const std::array<double, 2>& resCellSize, const std::array<double, 2>& half) {
    Dune::FieldVector<double, 2> globalCoord;
    globalCoord[0] = rasterCoord[0] * resCellSize[0] + half[0];
    globalCoord[1] = rasterCoord[1] * resCellSize[1] + half[1];
    return globalCoord;
}

/**
 * @brief Converts global coordinates to raster coordinates.
 * 
 * The function maps global coordinates to the raster coordinate system.
 * 
 * @param globalCoord The global coordinates.
 * @param cellSize The size of the grid cell in each dimension.
 * @param resCellSize The actual cell size.
 * @return Dune::FieldVector<double, 2> The corresponding raster coordinates.
 */
Dune::FieldVector<double, 2> convertToRasterCoordinates(const Dune::FieldVector<double, 2>& globalCoord, const std::array<double, 2>& cellSize, 
                                                            const std::array<double, 2>& resCellSize) {
    Dune::FieldVector<double, 2> rasterCoord;
    rasterCoord[0] = (globalCoord[0] - cellSize[0]/2)/(resCellSize[0]);
    rasterCoord[1] = (globalCoord[1] - cellSize[1]/2)/(resCellSize[1]);
    return rasterCoord;
}

/**
 * @brief Calculates boundaries for flow fragments.
 * 
 * The function computes boundary information for each flow fragment, including 
 * axis-aligned bounding boxes, directional properties, and adjusted start/end points.
 * 
 * @param fragments A vector of flow fragments to process.
 * @param minSizeFactor A scaling factor for the minimum boundary size (default is 0.4).
 * @return std::vector<fragmentBoundaries> The computed boundaries for all fragments.
 */
std::vector<fragmentBoundaries> calcFragmentBoundaries(std::vector<flowFragment> fragments, double minSizeFactor=0.4){
    
    std::vector<fragmentBoundaries> fragmentsBoundaries;

        for (auto f : fragments){ //store relevant information for each fragment so that it needs to be computed only once
        if(f.end==f.start){
            std::cerr << "Flow without end" << std::endl; //debug message
            f.end[0] -=0.001;
            f.end[1] -=0.001;
        }

        Dune::FieldVector<double, 2> flowVector = f.end - f.start;
        Dune::FieldVector<double, 2> normalFlowVector = calculateNormal(flowVector);

        if (f.widthEnd < 0) f.widthEnd = f.widthStart;

        Dune::FieldVector<double, 2> start1 = f.start + (f.widthStart/2) * normalFlowVector;
        Dune::FieldVector<double, 2> start2 = f.start - (f.widthStart/2) * normalFlowVector;
        Dune::FieldVector<double, 2> end1 = f.end + (f.widthEnd/2) * normalFlowVector;
        Dune::FieldVector<double, 2> end2 = f.end - (f.widthEnd/2) * normalFlowVector;
        
        if(abs(f.widthStart - std::sqrt(squaredDistance(start1, start2)))>1e-8){
            std::cerr << "width start: " << f.widthStart << " != " << std::sqrt(squaredDistance(start1, start2))  << std::endl;
        }

        Dune::FieldVector<double, 2> b1 = end1 - start1;
        Dune::FieldVector<double, 2> b2 = end2 - start2;


        int direction = 0;
        if(std::abs(flowVector[1]) < 1e-8 && f.widthStart==f.widthEnd) direction = 1;
        else if(std::abs(flowVector[0]) < 1e-8 && f.widthStart==f.widthEnd) direction = 2;
        else if((std::abs(flowVector[0] + flowVector[1]) < 1e-8 || std::abs(flowVector[0] - flowVector[1]) < 1e-8) && f.widthStart==f.widthEnd) direction=3;


        double minX = std::min({start1[0], start2[0], end1[0], end2[0]}); //axis aligned bounding box for each fragment
        double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
        double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
        double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});
        
        
        double minSize =  minSizeFactor*std::min(f.widthStart, f.widthEnd);

        fragmentBoundaries fB = fragmentBoundaries{start1, start2, b1, b2, minX, maxX, minY, maxY, direction, minSize};
        fB.depht = f.depht;
        fB.normal = normalFlowVector;
        fragmentsBoundaries.push_back(fB);
    }

    return fragmentsBoundaries;
}



/**
 * @brief Interpolates the elevation at a given point in the grid.
 * 
 * The function uses bilinear interpolation to estimate the elevation value 
 * at a specified point if this point is located between cells.
 * 
 * @param point The point in raster coordinates.
 * @param gridSize The size of the grid in x and y directions.
 * @param elevation_raster The raster dataset containing elevation values.
 * @return double The interpolated elevation at the specified point.
 */
double interpolateElevation(const Dune::FieldVector<double, 2>& point, const std::array<int, 2>& gridSize, const RasterDataSet<float>& elevation_raster) {
    auto [neighborX1, neighborX2] = calculateNeighborsOneDirection(point[0], gridSize[0]);
    auto [neighborY1, neighborY2] = calculateNeighborsOneDirection(point[1], gridSize[1]);

    if (neighborX2 == -2 && neighborY2 == -2) {
        return elevation_raster(neighborX1, neighborY1);
    } else if (neighborX2 == -2) {
        return (elevation_raster(neighborX1, neighborY1) + elevation_raster(neighborX1, neighborY2)) / 2;
    } else if (neighborY2 == -2) {
        return (elevation_raster(neighborX1, neighborY1) + elevation_raster(neighborX2, neighborY1)) / 2;
    } else {
        return (elevation_raster(neighborX1, neighborY1) + elevation_raster(neighborX1, neighborY2) +
                elevation_raster(neighborX2, neighborY1) + elevation_raster(neighborX2, neighborY2)) / 4;
    }
}

/**
 * @brief Determines whether a point lies within a flow fragment.
 * 
 * The function checks whether a given point lies within the boundaries 
 * of a flow fragment using axis-aligned bounding boxes and directional checks.
 * 
 * @param point The point to check.
 * @param fB The fragment boundaries containing the flow fragment's details.
 * @return bool true if the point lies within the flow fragment, false otherwise.
 */
bool isPointInFlow(const Dune::FieldVector<double, 2>& point, const fragmentBoundaries& fB) {
    if (point[0] < fB.minX || point[0] > fB.maxX || point[1] < fB.minY || point[1] > fB.maxY) {
        return false;
    }

    if (fB.direction == 0 || fB.direction == 3) {
        Dune::FieldVector<double, 2> vec1 = point - fB.start1;
        Dune::FieldVector<double, 2> vec2 = point - fB.start2;
        double cross1 = fB.b1[0] * vec1[1] - fB.b1[1] * vec1[0];
        double cross2 = fB.b2[0] * vec2[1] - fB.b2[1] * vec2[0];
        return !(cross1 > 0 || cross2 < 0);
    }

    return true;
}


std::vector<double> overallHeight(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
                                RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, 2> gridSize){

    std::array<double, 2> widthGrid = {cellSize[0] * (gridSize[0] - 1), cellSize[1] * (gridSize[1] - 1)}; //frame is gridSize*cellSize but on eighter side half a cell is missing
    std::array<double, 2> resCellSize = {widthGrid[0]/gridSize[0], widthGrid[1]/gridSize[1]}; //acctual cell size is slightly smaller than callSize

    std::vector<double> height(gridView.indexSet().size(2), 0);

    for(auto& v : vertices(gridView)){
        auto corner = v.geometry().corner(0);
        Dune::FieldVector<double, 2> refCorner = convertToRasterCoordinates(corner, cellSize, resCellSize);
        
        double elevation_val = interpolateElevation(refCorner, gridSize, elevation_raster);

        if(elevation_val < -5000 || elevation_val > 10000) elevation_val = 0; //wrong values

        height[gridView.indexSet().index(v)] = elevation_val;
    }
    return height;
}



std::vector<double> applyFlowHeightFragments(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
        std::vector<flowFragment> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, std::array<int, 2> gridSize){

    std::vector<double> originalHeight = height;

    Dune::FieldVector<double, 2> widthGrid = {cellSize[0] * (gridSize[0] - 1), cellSize[1] * (gridSize[1] - 1)};
    std::array<double, 2> resCellSize = {widthGrid[0]/gridSize[0], widthGrid[1]/gridSize[1]};
    std::array<double, 2> half = {cellSize[0]/2, cellSize[1]/2};

    std::vector<fragmentBoundaries> fragmentsBoundaries = calcFragmentBoundaries(fragments);

    
    for(auto& v : vertices(gridView)){ //check for each vertex if its part of a fragment
        auto corner = v.geometry().corner(0); //avoid recomputation
 
        for (const auto& fB : fragmentsBoundaries){

            if (!isPointInFlow(corner, fB)) continue;

            if(fB.direction == 3){ //diagonal //following not essential to work but smoothens results
                Dune::FieldVector<double, 2> refCorner = convertToRasterCoordinates(corner, cellSize, resCellSize);

                std::vector<int> neighbors= calcNeighbors(refCorner, gridSize);
                int neighborX1 = neighbors[0];
                int neighborX2 = neighbors[1];
                int neighborY1 = neighbors[2];
                int neighborY2 = neighbors[3];


                if(neighborX2==-2 && neighborY2 == -2){ // vertex in cell (of raster data)
                    
                    if(std::abs(fB.b1[0]-fB.b1[1])<1e-8){ //flow from bottom left to top right or vice versa
                        Dune::FieldVector<double, 2> bottomLeftCellCorner = convertToGlobalCoordinates({neighborX1, neighborY1}, resCellSize, half);
                        
                        if(!isPointInFlow(bottomLeftCellCorner, fB)){ //flow through corner of cell (bottom left corner not in flow)
                            double x_inCell = refCorner[0] - neighborX1;
                            double y_inCell = refCorner[1] - neighborY1;

                            int x1 = neighborX1;
                            int x2 = neighborX1;
                            int y1 = neighborY1;
                            int y2 = neighborY1;
                            if(x_inCell < y_inCell){ //flow through top left corner
                                x1 = neighborX1 - 1;
                                y2 = neighborY1 + 1;
                            } else { //flow through bottom right corner
                                x1 = neighborX1 + 1;
                                y2 = neighborY1 - 1;
                            }
                            x1 = std::clamp(x1, 0, gridSize[0]-1);
                            y2 = std::clamp(y2, 0, gridSize[1]-1);
                            height[gridView.indexSet().index(v)] = std::min(double(originalHeight[gridView.indexSet().index(v)]), 
                                                                                double((elevation_raster(x1, y1)+elevation_raster(x2, y2))/2)) - fB.depht;
                            continue;
                        }
                    }
                    else{//flow from bottom right to top left or vice versa
                        Dune::FieldVector<double, 2> bottomRightCellCorner = convertToGlobalCoordinates({neighborX1+1, neighborY1}, resCellSize, half);

                        if(!isPointInFlow(bottomRightCellCorner, fB)){ //flow through corner of cell (bottom right corner not in flow)
                            double x_inCell = refCorner[0]- neighborX1;
                            double y_inCell = refCorner[1]- neighborY1;

                            int x1 = neighborX1;
                            int x2 = neighborX1;
                            int y1 = neighborY1;
                            int y2 = neighborY1;
                            if(x_inCell + y_inCell > 1){ //flow through top right corner
                                x1 = neighborX1 + 1;
                                y2 = neighborY1 + 1;
                            } else { //flow through bottom left corner
                                x1 = neighborX1 - 1;
                                y2 = neighborY1 - 1;
                            }
                            x1 = std::clamp(x1, 0, gridSize[0]-1);
                            y2 = std::clamp(y2, 0, gridSize[1]-1);
                            height[gridView.indexSet().index(v)] = std::min(double(originalHeight[gridView.indexSet().index(v)]), 
                                                                                double((elevation_raster(x1, y1)+elevation_raster(x2, y2))/2)) - fB.depht;
                            continue;
                        }
                    }         
                } else if(neighborX2==-2 && neighborY2 != -2){ //vertex between two cells in y direction
                    height[gridView.indexSet().index(v)] = std::min(elevation_raster(neighborX1, neighborY1), elevation_raster(neighborX1, neighborY2)) - fB.depht;
                    continue;
                } else if(neighborX2 != -2 && neighborY2 == -2){ //vertex between two cells in x direction
                    height[gridView.indexSet().index(v)] = std::min(elevation_raster(neighborX1, neighborY1), elevation_raster(neighborX2, neighborY1)) - fB.depht;
                    continue;
                } else{ //vertex between four cells
                    height[gridView.indexSet().index(v)] = std::min((elevation_raster(neighborX1, neighborY1) + elevation_raster(neighborX2, neighborY2))/2, 
                        	                                                (elevation_raster(neighborX2, neighborY1) + elevation_raster(neighborX1, neighborY2))/2) - fB.depht;
                    continue;
                }
            }   
            height[gridView.indexSet().index(v)] = std::min(originalHeight[gridView.indexSet().index(v)] - fB.depht, height[gridView.indexSet().index(v)]);
        }
    }
    return height;
}

void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor, std::array<double, 2> cellSize, int maxIterations){
    
    typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    std::vector<fragmentBoundaries> fragmentsBoundaries = calcFragmentBoundaries(fragments, minSizeFactor);
    
    int c=0;
    bool change = true;
    while( c<maxIterations && change){ 
        std::cout << c << std::endl;
        c++; 
        change = false;
         
        for (const auto& element : elements(gridView)){ //check for each element if its part of the border of a flow --> refine
            
            auto corner0 = element.geometry().corner(0);
            auto corner1 = element.geometry().corner(1);
            auto corner2 = element.geometry().corner(2); 

            std::vector<Dune::FieldVector<double, 2>> corners = {corner0, corner1, corner2};
    
            for (const auto& fB : fragmentsBoundaries){

                bool between = false;
                bool out1 = false;
                bool out2 = false;

                bool continueFragment = false;
                int xOut = 0;
                int yOut = 0;
                for (int i = 0; i < 3; i++){

                    auto cornerI =corners[i];

                    if (cornerI[0] < fB.minX-cellSize[0] || cornerI[0] > fB.maxX+cellSize[0] || cornerI[1] < fB.minY-cellSize[1] || cornerI[1] > fB.maxY+cellSize[1]){  //element too far away
                        continueFragment = true; 
                        break;
                    }
                    if(cornerI[0] < fB.minX) xOut--; //test if (corner of) trangle is outside of desired range
                    else if (cornerI[0] > fB.maxX)xOut++;
                    if(cornerI[1] < fB.minY) yOut--;
                    else if (cornerI[1] > fB.maxY) yOut ++;

                    if(fB.direction == 1){ //horizontal
                        if(cornerI[1] > fB.maxY) out1 = true;
                        else if (cornerI[1] < fB.minY) out2 = true;
                        else between = true;
                        continue;
                    }

                    if(fB.direction == 2){ //vertical
                        if(cornerI[0] > fB.maxX) out2 = true;
                        else if (cornerI[0] < fB.minX) out1 = true;
                        else between = true;
                        continue;
                    }
        
                    Dune::FieldVector<double, 2> vec1 = cornerI - fB.start1;
                    Dune::FieldVector<double, 2> vec2 = cornerI - fB.start2;
                    double cross1 = fB.b1[0] * vec1[1] - fB.b1[1] * vec1[0];
                    double cross2 = fB.b2[0] * vec2[1] - fB.b2[1] * vec2[0];
                    
                    if(cross1 > 0) out1 = true;
                    else if(cross2 < 0) out2 = true;
                    else between = true;
                }


                
                if(continueFragment || std::abs(xOut) == 3 || std::abs(yOut) == 3 || between + out1 + out2 <= 1){continue;} // all 3 corners outside x/y boundary/one corner too far away
                if(out1+out2==2) { //whole flow inside triangle
                    grid->mark(1, element); 
                    change=true; 
                    continue;}

                double dist01 = squaredDistance(corner0, corner1);
                double dist02 = squaredDistance(corner0, corner2);

                if (std::min({dist01, dist02}) < fB.minSize*fB.minSize){ continue;} //small enough

                Dune::FieldVector<double, 2> hypo;
                Dune::FieldVector<double, 2> hypoCenter;

                //calc hypothenuse v; point on center of hypothenuse
                if (std::abs(dist01-dist02)< 1e-6){ //90deg at 0
                    hypo = corner2 - corner1;
                    hypoCenter = corner1 + hypo/2;
                }
                else if ( dist01 > dist02){  //90deg at 2
                    hypo = corner0 - corner1;
                    hypoCenter = corner1 + hypo/2;
                }
                else { //90deg at 1
                    hypo = corner0 - corner2;
                    hypoCenter = corner2 + hypo/2;
                }


                if(fB.direction == 1 || fB.direction == 2){ //flow (and border) horizontal or vertical (dot product with {1,0}/{0,1})
                    if (std::abs(hypo[0]) < 1e-8 || std::abs(hypo[1]) < 1e-8){ //hypotenuse vertical or horizontal
                        grid->mark(1, element);//hypothenuse vertical or horizontal
                        change = true;
                        continue;
                    }
                    if(fB.direction == 1 && (std::abs(fB.start1[1] - hypoCenter[1]) < fB.minSize/2 || std::abs(fB.start2[1] - hypoCenter[1]) < fB.minSize/2)){continue;} // border close to middle of triangle (in right angle)
                    else if(fB.direction == 2 && (std::abs(fB.start1[0] - hypoCenter[0]) < fB.minSize/2 || std::abs(fB.start2[0] - hypoCenter[0]) < fB.minSize/2)){ continue;} // border close to middle of triangle (in right angle)
                }
                else if(fB.direction == 3){ //flow diagonal
                    if (std::abs(hypo[0] + hypo[1]) < 1e-8 || std::abs(hypo[0] - hypo[1]) < 1e-8){ //hypothenuse diagonal
                        grid->mark(1, element);
                        change = true;
                        continue;
                    }
                    //checks if hypoCenter shifted by minSize/2 away from border in each both directions is on time one one side and one time on another
                    Dune::FieldVector<double, 2> minSizeHalf = {std::sqrt(2) * fB.minSize/2, std::sqrt(2) * fB.minSize/2};
                    if(abs(fB.normal[0]-fB.normal[1])<1e-8){
                        minSizeHalf = {std::sqrt(2) * -fB.minSize/2, std::sqrt(2)* fB.minSize/2 };
                    }
                                 
                    Dune::FieldVector<double, 2> vec1_U = hypoCenter + minSizeHalf - fB.start1;
                    Dune::FieldVector<double, 2> vec1_D = hypoCenter - minSizeHalf - fB.start1;
                    Dune::FieldVector<double, 2> vec2_U = hypoCenter + minSizeHalf - fB.start2;
                    Dune::FieldVector<double, 2> vec2_D = hypoCenter - minSizeHalf - fB.start2;

                    double cross1_U = fB.b1[0] * vec1_U[1] - fB.b1[1] * vec1_U[0];
                    double cross1_D = fB.b1[0] * vec1_D[1] - fB.b1[1] * vec1_D[0];
                    double cross2_U = fB.b2[0] * vec2_U[1] - fB.b2[1] * vec2_U[0];
                    double cross2_D = fB.b2[0] * vec2_D[1] - fB.b2[1] * vec2_D[0];
    
                    if((cross1_U > 0 && cross1_D < 0) ||(cross2_U > 0 && cross2_D < 0))continue; // points on different sides of border --> hypoCenter closer than minSize/2
                }
            
                change = true;
                grid->mark(1,element);
            }
        }
        
        grid->preAdapt();
        grid->adapt();  
        grid->postAdapt();            
         
    }
}


std::vector<flowFragment> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, 
                            std::array<double, 2> cellSize, std::array<int, 2> gridSize, double minAcc, double maxAccDiff, double scaleDephtFactor, 
                            double scaleWidthFactor, double minWidth){
    std::vector<flowFragment> rivers;
    int size = int(gridSize[0] * gridSize[1]);
    std::vector<int> skip(size, 0);


    const std::array<double, 2> resCellSize = {((gridSize[0]-1.0)/gridSize[0]) * cellSize[0], ((gridSize[1]-1.0)/gridSize[1]) * cellSize[1]};
    const std::array<double, 2> half ={cellSize[0]/2, cellSize[1]/2};

    std::cout << "Start finding flow fragments" << std::endl;
    for (int i = 0; i < gridSize[0]; i++){ 
        for (int j = 0; j < gridSize[1]; j++){
            if(skip[j*gridSize[0]+i]) continue;
            if(accumulation_raster(i, j) > minAcc){
                Dune::FieldVector<double, 2> start = {i, j};
                Dune::FieldVector<double, 2> end = start;

                int dir = int(direction_raster(i, j));
                int endI=  i;
                int endJ = j;
                int endDir = dir;
                if(dir>128 || dir==0) continue;
                double accStart = accumulation_raster(i, j);
                double accEnd = accStart;
                

                while (endDir == dir){ //straight flow
                    skip[endJ*gridSize[0]+endI]=1; //elements in middle of fragment don't need to be checked again
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
                    default:
                        break;
                
                    }
                    
                    if(end[0] >= gridSize[0] || end[1] >= gridSize[1] || end[0] < 0 || end[1] < 0) { //end is outside of grid
                        endI = round(currPoint[0]); 
                        endJ = round(currPoint[1]);
                        accEnd = accumulation_raster(endI, endJ); 
                        end=(currPoint + end)/2; //undo changes
                        break;
                    }
                    
                    endI = round(end[0]); 
                    endJ = round(end[1]);
                    endDir = int(direction_raster(endI, endJ));
                    accEnd = accumulation_raster(endI, endJ);

                    if(std::abs(accStart-accEnd)>maxAccDiff){
                        if(currPoint==start) break;
                        endI = round(currPoint[0]); 
                        endJ = round(currPoint[1]);
                        skip[endJ*gridSize[0]+endI]=0; //end should not be skipped
                        accEnd = accumulation_raster(endI, endJ); 
                        end = currPoint;
                        break;
                    }

                }

                skip[j*gridSize[0]+i] = 0; //start should not be skipped (skipped in first loop iteration)
                double volume = (accStart+accEnd)/2;
                double width = std::sqrt(volume/(scaleWidthFactor*(cellSize[0]+cellSize[1])/2)); //width of flow in cells
                double depht = std::sqrt(volume/(scaleDephtFactor*(cellSize[0]+cellSize[1])/2));
                width = std::min(1.0, width); //max width cellsize
                depht = std::min(15/((cellSize[0]+cellSize[1])/2), depht); //max depht 15m
                width = std::max(minWidth/((cellSize[0]+cellSize[1])/2), width); //min width
                depht = std::max(0.2/((cellSize[0]+cellSize[1])/2), depht); //min depht 0.2m
                //std::cout << "width: " << width << " --> "   << width*cellSize[0] << std::endl;
                //std::cout << "depht: " <<depht << " --> "    << depht*cellSize[0] << std::endl << std::endl;

                //std::cout << "final: " << start << " - " << end << std::endl;
                Dune::FieldVector<double, 2>shift = {0.5, 0.5}; //fragmetns should start in middle of cell
                start = convertToGlobalCoordinates(start+shift, resCellSize, half);
                end = convertToGlobalCoordinates(end+shift, resCellSize, half);

                
                flowFragment f = {start, end, width * cellSize[0]}; 
                if(dir==4 || dir==64) f.widthStart= width * cellSize[1];
                f.depht = depht * cellSize[0];
                rivers.push_back(f);
            }
        }
    }

    for (int i = rivers.size()-1; i>=0; i--){ //delete all fragments that can be skipped but were not skipped due to order
        auto start = rivers[i].start;
        int start_index_x = round(start[0]/(cellSize[0] * ((gridSize[0]-1.0)/gridSize[0])) -1);
        int start_index_y = round(start[1]/(cellSize[1] * ((gridSize[1]-1.0)/gridSize[1])) -1);
        if (skip[start_index_y*gridSize[0]+start_index_x]) rivers.erase(rivers.begin()+i);
    }

    std::cout << "\n Find fragments done. " << std::endl;

    return rivers;
}

RasterDataSet<float> removeUpwardsRivers(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                            std::array<double, 2> cellSize, std::array<int, 2> gridSize, double minAcc){



    std::cout << "\nStart removing upward rivers" << std::endl;
    bool change = true;
    while(change){
        change=false;
    
        for (int i = 0; i < gridSize[0]; i++){ 
            for (int j = 0; j < gridSize[1]; j++){
                if(accumulation_raster(i, j) > minAcc){
                    Dune::FieldVector<int, 2> start = {i, j};
                    Dune::FieldVector<int, 2> end = start;

                    int dir = int(direction_raster(i, j));
                    if(dir>128 || dir==0) continue;

                    switch (dir)
                    {
                    case 1:     end = end + Dune::FieldVector<int, 2>{1, 0};
                        break;
                    case 2:     end = end + Dune::FieldVector<int, 2>{1,-1};
                        break;
                    case 4:     end = end + Dune::FieldVector<int, 2>{0, -1};
                        break;
                    case 8:     end = end + Dune::FieldVector<int, 2>{-1, -1};
                        break;
                    case 16:    end = end + Dune::FieldVector<int, 2>{-1, 0};
                        break;
                    case 32:    end = end + Dune::FieldVector<int, 2>{-1, 1};
                        break;
                    case 64:    end = end + Dune::FieldVector<int, 2>{0, 1};
                        break;
                    case 128:   end = end + Dune::FieldVector<int, 2>{1, 1};
                        break;
                    default:
                        break;
                    }
                    //std::cout << currPoint <<", " << elevation_raster(currPoint[0], currPoint[1]) << std::endl;

                    if(end[0] >= gridSize[0] || end[1] >= gridSize[1] ||  end[0] < 0 || end[1] < 0 || accumulation_raster(end[0], end[1]) < minAcc){
                        break;}
                    
                    if(elevation_raster(end[0], end[1]) > elevation_raster(start[0], start[1])){ //end is higher than start
                        elevation_raster(end[0], end[1]) = elevation_raster(start[0], start[1]);
                        change = true;
                        std::cout << "change: "  << end << std::endl;
                    }
                    
                }
            }
        }
    }

    std::cout << "\nRemove upward rivers done. " << std::endl;

    return elevation_raster;
}


std::vector<double> addRiversToMap(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, std::array<double, 2> cellSize, std::array<int, 2> gridSize,
                                    RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                                    double minSizeFactor, double minAcc, double maxAccDiff, double scaleDephtFactor, double scaleWidthFactor, int maxIterations, 
                                    double minWidth){
    typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    elevation_raster = removeUpwardsRivers(accumulation_raster, direction_raster, elevation_raster, cellSize, gridSize, minAcc);
    std::vector<flowFragment> rivers = detectFragments(accumulation_raster, direction_raster, cellSize, gridSize, minAcc, maxAccDiff, scaleDephtFactor, scaleWidthFactor, minWidth);

    refineGridwithFragments(grid, rivers, minSizeFactor, cellSize, maxIterations);

    std::vector<double> height = overallHeight(gridView, elevation_raster, cellSize, gridSize);
    height = applyFlowHeightFragments(gridView, rivers, height, elevation_raster, cellSize, gridSize);
    return height;
}
