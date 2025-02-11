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

double diagDistance(const Dune::FieldVector<double, 2> start, const Dune::FieldVector<double, 2> vec, const Dune::FieldVector<double, 2> point){
    const Dune::FieldVector<double, 2> startToPoint = point - start;
    double d = std::copysign(1, vec[0]) * startToPoint[0]- std::copysign(1, vec[1]) * startToPoint[1];
    double dist = d / sqrt(2);
    return abs(dist);
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
 * @param realCellSize The size of the grid cell in each dimension.
 * @param cellSize The size of the grid cell in each dimension.
 * 
 * @return Dune::FieldVector<double, 2> The corresponding global coordinates.
 */
Dune::FieldVector<double, 2> convertToGlobalCoordinates(const Dune::FieldVector<double, 2>& rasterCoord, const std::array<double, 2>& cellSize, 
                                                        const std::array<double, 2>& realCellSize) {
    Dune::FieldVector<double, 2> globalCoord;
    globalCoord[0] = rasterCoord[0] * realCellSize[0] + cellSize[0]/2;
    globalCoord[1] = rasterCoord[1] * realCellSize[1] + cellSize[1]/2;
    return globalCoord;
}


/**
 * @brief Converts global coordinates to raster coordinates.
 * 
 * The function maps global coordinates to the raster coordinate system.
 * 
 * @param globalCoord The global coordinates.
 * @param cellSize The size of the grid cell in each dimension.
 * @param realCellSize The actual cell size.
 * @return Dune::FieldVector<double, 2> The corresponding raster coordinates.
 */
Dune::FieldVector<double, 2> convertToRasterCoordinates(const Dune::FieldVector<double, 2>& globalCoord, const std::array<double, 2>& cellSize, 
                                                            const std::array<double, 2>& realCellSize) {
    Dune::FieldVector<double, 2> rasterCoord;
    rasterCoord[0] = (globalCoord[0] - cellSize[0]/2)/(realCellSize[0]);
    rasterCoord[1] = (globalCoord[1] - cellSize[1]/2)/(realCellSize[1]);
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
            std::cerr << "Wrong boderders: width start: " << f.widthStart << " != distance boundaries" << std::sqrt(squaredDistance(start1, start2))  << std::endl;
        }

        Dune::FieldVector<double, 2> b1 = end1 - start1;
        Dune::FieldVector<double, 2> b2 = end2 - start2;


        int direction = 0;
        if(std::abs(flowVector[1]) < 1e-8 && f.widthStart==f.widthEnd) direction = 1; //horizontal
        else if(std::abs(flowVector[0]) < 1e-8 && f.widthStart==f.widthEnd) direction = 2; //vetical
        else if((std::abs(flowVector[0] + flowVector[1]) < 1e-8 || std::abs(flowVector[0] - flowVector[1]) < 1e-8) && f.widthStart==f.widthEnd) direction=3; //diagonal


        double minX = std::min({start1[0], start2[0], end1[0], end2[0]}); //axis aligned bounding box for each fragment
        double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
        double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
        double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});

        if(direction==30){ //TODO: delete
            minX = std::min({f.start[0], f.end[0]}); //axis aligned bounding box for each fragment
            maxX = std::max({f.start[0], f.end[0]});
            minY = std::min({f.start[1], f.end[1]});
            maxY = std::max({f.start[1], f.end[1]});
        }
        
        
        double minSize =  minSizeFactor*std::min(f.widthStart, f.widthEnd);

        fragmentBoundaries fB = fragmentBoundaries{start1, start2, b1, b2, minX, maxX, minY, maxY, direction, minSize};
        fB.depht = f.depht;
        fB.normal = normalFlowVector;
        fragmentsBoundaries.push_back(fB);
        //std::cout << sqrt(squaredDistance(start1, start2)) << " " << sqrt(squaredDistance(end1, end2)) << std::endl;
        //std::cout <<"start: " << start1 << " ; " << start2 << " | end: " << end1 << " ; " << end2 << std::endl;

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

    std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize);

    std::vector<double> height(gridView.indexSet().size(2), 0);

    for(auto& v : vertices(gridView)){
        auto corner = v.geometry().corner(0);
        Dune::FieldVector<double, 2> rasterCorner = convertToRasterCoordinates(corner, cellSize, realCellSize);
        
        double elevation_val = interpolateElevation(rasterCorner, gridSize, elevation_raster);

        if(elevation_val < -5000 || elevation_val > 10000) elevation_val = 0; //wrong values

        height[gridView.indexSet().index(v)] = elevation_val;
    }
    return height;
}



std::vector<double> applyFlowHeightFragments(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
       std::vector<std::vector<flowFragment>> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, 
        std::array<int, 2> gridSize){

    std::vector<double> originalHeight = height;

    std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize);

    Dune::FieldVector<double, 2> bisectors = convertToGlobalCoordinates({std::floor(gridSize[0]/2), std::floor(gridSize[0]/2)}, cellSize, realCellSize);
    double xBisector = bisectors[0];
    double yBisector = bisectors[1];

    std::vector<fragmentBoundaries> fragmentsBoundaries00 = calcFragmentBoundaries(fragments[0]);
    std::vector<fragmentBoundaries> fragmentsBoundaries01 = calcFragmentBoundaries(fragments[1]);
    std::vector<fragmentBoundaries> fragmentsBoundaries10 = calcFragmentBoundaries(fragments[2]);
    std::vector<fragmentBoundaries> fragmentsBoundaries11 = calcFragmentBoundaries(fragments[3]);

    std::vector<std::vector<fragmentBoundaries>> fragmentsBoundaries = {fragmentsBoundaries00, fragmentsBoundaries01, fragmentsBoundaries10, fragmentsBoundaries11};


    
    for(auto& v : vertices(gridView)){ //check for each vertex if its part of a fragment
        auto corner = v.geometry().corner(0); //avoid recomputation

        int selectFB = 0;
        if(corner[0] >= xBisector) selectFB += 2;
        if(corner[1] >= yBisector) selectFB += 1;
        
 
        for (const auto& fB : fragmentsBoundaries[selectFB]){

            if (!isPointInFlow(corner, fB)) continue;

            if(fB.direction == 3){ //diagonal //following not essential to work but smoothens results
                Dune::FieldVector<double, 2> refCorner = convertToRasterCoordinates(corner, cellSize, realCellSize);

                std::vector<int> neighbors= calcNeighbors(refCorner, gridSize);
                int neighborX1 = neighbors[0];
                int neighborX2 = neighbors[1];
                int neighborY1 = neighbors[2];
                int neighborY2 = neighbors[3];


                if(neighborX2==-2 && neighborY2 == -2){ // vertex in cell (of raster data)
                    
                    if(std::abs(fB.b1[0]-fB.b1[1])<1e-8){ //flow from bottom left to top right or vice versa
                        Dune::FieldVector<double, 2> bottomLeftCellCorner = convertToGlobalCoordinates({neighborX1, neighborY1}, cellSize, realCellSize);
                        
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
                        Dune::FieldVector<double, 2> bottomRightCellCorner = convertToGlobalCoordinates({neighborX1+1, neighborY1}, cellSize, realCellSize);

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
        std::vector<std::vector<flowFragment>> fragments, std::array<int, 2> gridSize, double minSizeFactor, std::array<double, 2> cellSize, int maxIterations, double minMinSize){
    
    std::cout << "Refining";
    typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize); //TODO check if gridView.Sze=gridsize

    Dune::FieldVector<double, 2> bisectors = convertToGlobalCoordinates({std::floor(gridSize[0]/2), std::floor(gridSize[0]/2)}, cellSize, realCellSize);
    double xBisector = bisectors[0];
    double yBisector = bisectors[1];

    std::vector<fragmentBoundaries> fragmentsBoundaries00 = calcFragmentBoundaries(fragments[0], minSizeFactor);
    std::vector<fragmentBoundaries> fragmentsBoundaries01 = calcFragmentBoundaries(fragments[1], minSizeFactor);
    std::vector<fragmentBoundaries> fragmentsBoundaries10 = calcFragmentBoundaries(fragments[2], minSizeFactor);
    std::vector<fragmentBoundaries> fragmentsBoundaries11 = calcFragmentBoundaries(fragments[3], minSizeFactor);

    std::vector<std::vector<fragmentBoundaries>> fragmentsBoundaries = {fragmentsBoundaries00, fragmentsBoundaries01, fragmentsBoundaries10, fragmentsBoundaries11};

    if(minMinSize >0){
        for (auto& fragmentBoundarieI : fragmentsBoundaries){
            for(auto& fB : fragmentBoundarieI){
                if(fB.minSize*(realCellSize[0]+realCellSize[1])/2 < minMinSize ) fB.minSize = minMinSize/(realCellSize[0]+realCellSize[1])*2; //minimum discripancy 0.2m
                //fB.minSize = 7;
            }
        }
    }
    
    int c=0;
    bool change = true;
    while( c<8 && change){ 
        std::cout << ".";
        c++; 
        change = false;
         
        for (const auto& element : elements(gridView)){ //check for each element if its part of the border of a flow --> refine
            //std::cout << "\nelem: ";
            auto corner0 = element.geometry().corner(0);
            auto corner1 = element.geometry().corner(1);
            auto corner2 = element.geometry().corner(2); 

            std::vector<Dune::FieldVector<double, 2>> corners = {corner0, corner1, corner2};

            int selectFB = 0;
            for (auto c : corners){
                if(c[0] == xBisector) continue;
                if(c[0] > xBisector){
                    selectFB += 2;
                }
                break;
            }
            for (auto c : corners){
                if(c[1] == yBisector) continue;
                if(c[1] > yBisector){
                    selectFB += 1;
                }
                break;
            }

            bool marked = false;
            //std::cout << "selectFB: " << selectFB << " ";

            for (const auto& fB : fragmentsBoundaries[selectFB]){
                //std::cout << "\nfragment boundaries: \n" << "start: " << fB.start1 << " ; " << fB.start2 << "\nb: " <<fB.b1 << ", " << fB.b2 <<"\n"  << std::endl;
                if(marked) { break;}
                //if(passt) std::cout << "_";

                bool between = false;
                bool out1 = false;
                bool out2 = false;

                bool continueFragment = false;
                int xOut = 0;
                int yOut = 0;
                
                for (int i = 0; i < 3; i++){

                    auto cornerI =corners[i];

                    if (cornerI[0] < fB.minX-realCellSize[0] || cornerI[0] > fB.maxX+realCellSize[0] || cornerI[1] < fB.minY-realCellSize[1] || cornerI[1] > fB.maxY+realCellSize[1]){  //element too far away
                        continueFragment = true; 
                        break;
                    }
                    
                    //TODO: change this
                    //if (abs(cornerI[0]-  fB.minX) < 1e-4 || abs(cornerI[0] - fB.maxX) < 1e-4){ //corner on bounding box
                    //    continueFragment = true; 
                    //    break;
                    //}  
                    //else 
                    if(cornerI[0] < fB.minX) xOut--; //test if (corner of) triangle is outside of desired range
                    else if (cornerI[0] > fB.maxX)xOut++;

                    //TODO: change this
                    //if (abs(cornerI[1]- fB.minY) < 1e-4 || abs(cornerI[1] - fB.maxY) < 1e-4){ //corner on bounding box
                    //    continueFragment = true; 
                    //    break;
                    //}
                    //else 
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

                    //std::cout << "corner: " << cornerI  <<", " << out1 << out2 << between << " cross1: " << cross1 << " cross2: " << cross2 << std::endl;
                }


                
                if(continueFragment || std::abs(xOut) == 3 || std::abs(yOut) == 3 || between + out1 + out2 <= 1){continue;} // all 3 corners outside x/y boundary/one corner too far away
                //std::cout << "corners: " << corners[0] << " \n       " << corners[1] << " \n       " << corners[2] << std::endl;
                //std::cout << "out1: " << out1 << " out2: " << out2 << " between: " << between << std::endl;
                if(out1+out2==2) { //whole flow inside triangle
                    //std::cout <<"out1+out2==2" << std::endl;
                    grid->mark(1, element); 
                    marked = true;
                    change=true; 
                    continue;}

                double dist01 = squaredDistance(corner0, corner1);
                double dist02 = squaredDistance(corner0, corner2);

                if (std::min({dist01, dist02}) < fB.minSize*fB.minSize){ continue;} //small enough



                //following makes only small difference
                if(fB.direction == 1){
                    Dune::FieldVector<double, 2> pointOfHorSide = {-100, -100};
                    if(abs((corner2-corner0)[1]) < 1e-8) pointOfHorSide = corner0;
                    else if (abs((corner2-corner1)[1]) < 1e-8) pointOfHorSide = corner1;
                    else if (abs((corner1-corner0)[1]) < 1e-8) pointOfHorSide = corner0;
                    if(abs((pointOfHorSide -fB.start1)[1]) < 1e-1){
                        //std::cout << ":"; 
                        continue;} //TODO: <?; restrict size of triangle?
                    else if(abs((pointOfHorSide -fB.start2)[1]) < 1e-1) {
                        //std::cout << ":"; 
                        continue;} //TODO: <?; restrict size of triangle?

                }
                if(fB.direction == 2){
                    Dune::FieldVector<double, 2> pointOfVertSide = {-100, -100};
                    if(abs((corner2-corner0)[0]) < 1e-8) pointOfVertSide = corner0;
                    else if (abs((corner2-corner1)[0]) < 1e-8) pointOfVertSide = corner1;
                    else if (abs((corner1-corner0)[0]) < 1e-8) pointOfVertSide = corner0;
                    if(abs((pointOfVertSide -fB.start1)[0]) < 1e-1) {
                        //std::cout << ":"; 
                        continue;} //TODO: <?; restrict size of triangle?
                    else if(abs((pointOfVertSide -fB.start2)[0]) < 1e-1) {
                        //std::cout << ":"; 
                        continue;} //TODO: <?; restrict size of triangle?

                }
                if(fB.direction == 3){
                    Dune::FieldVector<double, 2> pointOfDiagSide = {-100, -100};
                    Dune::FieldVector<double, 2> side1 = (corner2 - corner0);
                    Dune::FieldVector<double, 2> side2 = (corner2 - corner1);
                    Dune::FieldVector<double, 2> side3 = (corner1 - corner0);             

                    side1[0]=side1[0]*fB.b1[0];
                    side2[0]=side2[0]*fB.b1[0];
                    side3[0]=side3[0]*fB.b1[0];
                    side1[1]=side1[1]*fB.b1[1];
                    side2[1]=side2[1]*fB.b1[1];
                    side3[1]=side3[1]*fB.b1[1];

                    if(abs(side1[0]-side1[1]) < 1e-8) pointOfDiagSide = corner0;
                    if(abs(side2[0]-side2[1]) < 1e-8) pointOfDiagSide = corner1;
                    if(abs(side3[0]-side3[1]) < 1e-8) pointOfDiagSide = corner0;


                    if(diagDistance(fB.start1, fB.b1, pointOfDiagSide) < 1e-1 || diagDistance(fB.start2, fB.b1, pointOfDiagSide) < 1e-1) {
                        //std::cout << ";"; 
                        continue;} //TODO: <?; restrict size of triangle?
                }

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


                if(fB.direction == 1 || fB.direction == 2){ //flow (and border) horizontal or vertical
                    if (std::abs(hypo[0]) < 1e-8 || std::abs(hypo[1]) < 1e-8){ //hypotenuse vertical or horizontal

                        //grid->mark(1, element);//hypothenuse vertical or horizontal
                        //marked = true;
                        //change = true;
                        ////std::cout << "wrong direction" << std::endl;
                        //continue;
                    }
                    if(fB.direction == 1 && (std::abs(fB.start1[1] - hypoCenter[1]) < fB.minSize/2 || std::abs(fB.start2[1] - hypoCenter[1]) < fB.minSize/2)){continue;} //horizontal: border close to middle of triangle (in right angle)
                    else if(fB.direction == 2 && (std::abs(fB.start1[0] - hypoCenter[0]) < fB.minSize/2 || std::abs(fB.start2[0] - hypoCenter[0]) < fB.minSize/2)){continue;} //vertical: border close to middle of triangle (in right angle)
                }
                else if(fB.direction == 3){ //flow diagonal
                    if (std::abs(hypo[0] + hypo[1]) < 1e-8 || std::abs(hypo[0] - hypo[1]) < 1e-8){ //hypothenuse diagonal
                        //grid->mark(1, element);
                        ////std::cout << "hypo diagonal" << std::endl;
                        //marked = true;
                        //change = true;
                        //continue;
                    }
                    
                    //checks if hypoCenter shifted by minSize/2 away from border in each both directions is on time one one side and one time on another
                    //Dune::FieldVector<double, 2> minSizeHalf = {std::sqrt(2) * fB.minSize/2, std::sqrt(2) * fB.minSize/2};
                    //if(abs(fB.normal[0]+fB.normal[1])<1e-8){
                    //    minSizeHalf = {std::sqrt(2) * -fB.minSize/2, std::sqrt(2)* fB.minSize/2 };
                    //}
                    //if(hypoCenter[0]<90 && hypoCenter[1] > 195){
                    //    std::cout << "\nfragment boundaries: \n" << "start: " << fB.start1 << " ; " << fB.start2 << "\nb: " <<fB.b1 << ", " << fB.b2  << std::endl;
//
                    //    std::cout << "\n hypoCenter= " << hypoCenter << "\nhc-msh = " << hypoCenter-minSizeHalf << "\nhc+msh = " << hypoCenter + minSizeHalf << std::endl;             
                    //    std::cout <<"start-hc: x: " << fB.start1[0] - hypoCenter[0] << " y: " << fB.start1[1] - hypoCenter[1] << "\nstart-hc: x: " << fB.start2[0] - hypoCenter[0] << " y: " << fB.start2[1] - hypoCenter[1] << std::endl;
                    //}
                    Dune::FieldVector<double, 2> vec = (hypoCenter) - fB.start1;

                    //Dune::FieldVector<double, 2> vec1_U = (hypoCenter + minSizeHalf) - fB.start1;
                    //Dune::FieldVector<double, 2> vec1_D = (hypoCenter - minSizeHalf) - fB.start1;
                    //Dune::FieldVector<double, 2> vec2_U = (hypoCenter + minSizeHalf) - fB.start2;
                    //Dune::FieldVector<double, 2> vec2_D = (hypoCenter - minSizeHalf) - fB.start2;

                    //if(hypoCenter[0]<90 && hypoCenter[1] > 195) std::cout << "vec1_U: " << vec1_U << " vec1_D: " << vec1_D << " vec2_U: " << vec2_U << " vec2_D: " << vec2_D << std::endl;
                    
                    //double cross1_U = fB.b1[0] * vec1_U[1] - fB.b1[1] * vec1_U[0];
                    //double cross1_D = fB.b1[0] * vec1_D[1] - fB.b1[1] * vec1_D[0];
                    //double cross2_U = fB.b2[0] * vec2_U[1] - fB.b2[1] * vec2_U[0];
                    //double cross2_D = fB.b2[0] * vec2_D[1] - fB.b2[1] * vec2_D[0];
                    //if(hypoCenter[0]<90 && hypoCenter[1] > 195) std::cout << fB.b1[0] <<"*"<< vec1_U[1] <<"-"<< fB.b1[1] <<"*" << vec1_U[0] << "=" << cross1_U << std::endl;
                    //if(hypoCenter[0]<90 && hypoCenter[1] > 195) std::cout << "cross1_U: " << cross1_U << " cross1_D: " << cross1_D << " cross2_U: " << cross2_U << " cross2_D: " << cross2_D << std::endl;
                    //if(passt)std::cout << ".";
                    //std::cout << fB.start1 <<", " << fB.b1 << "\n" << hypoCenter << ": "<< hypoCenter- fB.start1 << std::endl;
                    //std::cout << diagDistance(fB.start1, fB.b1, hypoCenter) << ", " << diagDistance(fB.start2, fB.b1, hypoCenter) << std::endl;
                    if(abs(diagDistance(fB.start1, fB.b1, hypoCenter))<fB.minSize/2 || abs(diagDistance(fB.start2, fB.b1, hypoCenter))<fB.minSize/2){continue;}
                    //if(abs(abs(vec[0])-abs(vec[1]))<fB.minSize*sqrt(2)) std::cout << "a";
                    //if((cross1_U * cross1_D < 0) ||(cross2_U * cross2_D < 0) ){std::cout << "b\n";continue;} // points on different sides of border --> hypoCenter closer than minSize/2
                    //if(abs(cross1_U) < 1e-8 || abs(cross1_D) < 1e-8 || abs(cross2_U) < 1e-8 || abs(cross2_D) < 1e-8){std::cout << "close enough" << std::endl; passt=true;continue;} //points on border
                }
            
                //std::cout << "end" << std::endl;    
                if(squaredDistance(fB.start1, hypoCenter) < 625 || squaredDistance(fB.start1, hypoCenter) < 625) continue; //TODO check effect of this
                grid->mark(1,element);
                marked = true;
                change = true;
            }
        }
        
        grid->preAdapt();
        grid->adapt();  
        grid->postAdapt();            
         
    }
}


std::vector<std::vector<flowFragment>> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, 
                            std::array<double, 2> cellSize, std::array<int, 2> gridSize, double minAcc, double maxAccDiff, bool fixedWidth, double scaleDephtFactor, 
                            double scaleWidthFactor, double minWidth){
    int size = int(gridSize[0] * gridSize[1]);
    std::vector<int> skip(size, 0);

    int xBisector = std::floor(gridSize[0]/2);
    int yBisector = std::floor(gridSize[1]/2);


    std::vector<flowFragment> rivers00;
    std::vector<flowFragment> rivers01;
    std::vector<flowFragment> rivers10;
    std::vector<flowFragment> rivers11;

    const std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize);

    std::cout << "Start finding flow fragments" << std::endl;
    for (int i = 0; i < gridSize[0]; i++){ 
        for (int j = 0; j < gridSize[1]; j++){
            if(skip[j*gridSize[0]+i]) continue;
            if(accumulation_raster(i, j) > minAcc ){
                Dune::FieldVector<double, 2> start = {i, j};
                Dune::FieldVector<double, 2> end = start;

                int dir = int(direction_raster(i, j));
                if(dir>128 || dir==0) continue;
                
                int endI=  i;
                int endJ = j;
                int endDir = dir;
                
                double accStart = accumulation_raster(i, j);
                double accEnd = accStart;

                //if(dir != 2 && dir != 8 && dir != 32 && dir != 128) continue; //only horizontal or vertical flow
                

                while (endDir == dir ){ //straight flow
                    skip[endJ*gridSize[0]+endI]=1; //elements in middle of fragment don't need to be checked again
                    Dune::FieldVector<double, 2> currPoint = end;
                    switch (dir)
                    {
                    case 1:     end= end+ Dune::FieldVector<double, 2>{1, 0};
                        break;
                    case 2:     end= end+ Dune::FieldVector<double, 2>{1,-1};
                        break;
                    case 4:     end= end+ Dune::FieldVector<double, 2>{0, -1};
                        break;
                    case 8:     end= end+ Dune::FieldVector<double, 2>{-1, -1};
                        break;
                    case 16:    end= end+ Dune::FieldVector<double, 2>{-1, 0};
                        break;
                    case 32:    end= end+ Dune::FieldVector<double, 2>{-1, 1};
                        break;
                    case 64:    end= end+ Dune::FieldVector<double, 2>{0, 1};
                        break;
                    case 128:   end= end+ Dune::FieldVector<double, 2>{1, 1};
                        break;
                    default:
                        break;
                
                    }
                    
                    if(end[0] >= gridSize[0] || end[1] >= gridSize[1] || end[0] < 0 || end[1] < 0) { //end is outside of grid
                        accEnd = accumulation_raster(endI, endJ); 
                        end=(currPoint + end)/2; //set end on border of grid
                        break;
                    }
                    
                    endI = round(end[0]); 
                    endJ = round(end[1]);
                    endDir = int(direction_raster(endI, endJ));
                    accEnd = accumulation_raster(endI, endJ);

                    if(std::abs(accStart-accEnd)>maxAccDiff){
                        if(currPoint==start) break;
                        end = currPoint;
                        endI = round(currPoint[0]); 
                        endJ = round(currPoint[1]);
                        skip[endJ*gridSize[0]+endI]=0; //end should not be skipped
                        accEnd = accumulation_raster(endI, endJ); 
                        break;
                    }

                }

                skip[j*gridSize[0]+i] = 0; //start should not be skipped (skipped in first loop iteration)

                double volume = (accStart+accEnd)/2;
                double depht = volume/scaleDephtFactor;
                depht = std::min(15.0, depht); //max depht 15m
                depht = std::max(0.2, depht); //min depht 0.2m

                depht = 90;

                double width = volume/(scaleWidthFactor); //width of flow in cells

                //std::cout << "accStart " << accStart << " accEnd " << accEnd << std::endl;

                //std::cout <<"breite: " << width << std::endl;

                if(fixedWidth){   
                    
                    if (volume < 500) width = 9.6;//1/8*90.0; // 9.4;
                    else if (volume < 2500) width = 46.355; //46.2;
                    else if (volume < 5000) width = 65.5;
                    else width = 90;
                    //else continue; //river big enough to not be refined
                    //std::cout << "fixed width " << width << std::endl;
                }
                else{
                    width = volume/(scaleWidthFactor);
                    width = std::min((realCellSize[0]+realCellSize[1])/2, width); //max width realCellSize
                    width = std::max(minWidth, width); //min width
                }
    

                //std::cout << start << " - " << end << std::endl;
                //std::cout << "width: " << width << " --> "   << width*realCellSize[0] <<" vs: "<< volume/scaleWidthFactor << std::endl;
                //std::cout << "depht: " <<depht << " --> "    << depht*realCellSize[0] << " vs: "<< volume/scaleDephtFactor << std::endl << std::endl;

                //std::cout << "final: " << start << " - " << end << std::endl;

                std::array<bool, 2> boxStart = {start[0] < xBisector, start[1] < yBisector};
                std::array<bool, 2> boxEnd = {end[0] < xBisector, end[1] < yBisector};

                Dune::FieldVector<double, 2>shift = {0.5, 0.5}; //fragments should start in middle of cell
                start = convertToGlobalCoordinates(start+shift, cellSize, realCellSize);
                end = convertToGlobalCoordinates(end+shift, cellSize, realCellSize);
                width = width/cellSize[0] * realCellSize[0];
                //std::cout << width << std::endl;
                
                flowFragment f = {start, end, width}; 
                f.depht = depht;

                if((boxStart[0] xor boxEnd[0]) && (boxStart[1] xor boxEnd[1])){ 
                    rivers00.push_back(f);
                    rivers01.push_back(f);
                    rivers10.push_back(f);
                    rivers11.push_back(f);
                    continue;
                }

                if((boxStart[0] && boxStart[1])|| (boxEnd[0] && boxEnd[1])) rivers00.push_back(f);
                if((boxStart[0] && !boxStart[1])|| (boxEnd[0] && !boxEnd[1])) rivers01.push_back(f);
                if((!boxStart[0] && boxStart[1])|| (!boxEnd[0] && boxEnd[1])) rivers10.push_back(f);
                if((!boxStart[0] && !boxStart[1])|| (!boxEnd[0] && !boxEnd[1])) rivers11.push_back(f);                
            }
        }
        
    }

    for (auto rivers : {rivers00, rivers01, rivers10, rivers11}){ //delete all fragments that can be skipped but were not skipped due to order	
        for (int i = rivers.size()-1; i>=0; i--){ //delete all fragments that can be skipped but were not skipped due to order
            auto start = rivers[i].start;
            int start_index_x = round(start[0]/(realCellSize[0] * ((gridSize[0]-1.0)/gridSize[0])) -1);
            int start_index_y = round(start[1]/(realCellSize[1] * ((gridSize[1]-1.0)/gridSize[1])) -1);
            if (skip[start_index_y*gridSize[0]+start_index_x]) rivers.erase(rivers.begin()+i);
        }
    }

    std::cout << "Find fragments done. " << std::endl;

    return {rivers00, rivers01, rivers10, rivers11};
}

RasterDataSet<float> removeUpwardsRivers(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                             std::array<int, 2> gridSize, double minAcc){


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

    std::cout << "Remove upward rivers done. " << std::endl;

    return elevation_raster;
}


std::vector<double> addRiversToMap(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, 
                                    const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView,
                                    std::array<double, 2> cellSize, std::array<int, 2> gridSize,
                                    RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                                    double minSizeFactor, double minAcc, double maxAccDiff, double scaleDephtFactor, double scaleWidthFactor, int maxIterations, 
                                    double minWidth, double minMinSize){



    elevation_raster = removeUpwardsRivers(accumulation_raster, direction_raster, elevation_raster, gridSize, minAcc);

    std::vector<std::vector<flowFragment>> fragments = detectFragments(accumulation_raster, direction_raster, cellSize, gridSize, minAcc, maxAccDiff, true, scaleDephtFactor, scaleWidthFactor, minWidth);
    refineGridwithFragments(grid, fragments, gridSize, minSizeFactor, cellSize, maxIterations);
    
    std::vector<double> height = overallHeight(gridView, elevation_raster, cellSize, gridSize);

    height= applyFlowHeightFragments(gridView, fragments, height, elevation_raster,cellSize, gridSize);
    return height;
}

std::array<double, 2> calcRealCellSize(std::array<double, 2> cellSize, std::array<int, 2> gridSize){
    return {cellSize[0] * (gridSize[0]-1)/gridSize[0], cellSize[1] * (gridSize[1]-1)/gridSize[1]};
}