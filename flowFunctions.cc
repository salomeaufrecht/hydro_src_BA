#include "flowFunctions.hh"

/**
 * @brief Calculates the real cell size based on the grid size and cell size.
 * 
 * @param cellSize The input size of each grid cell.
 * @param gridSize The dimensions of the grid.
 * 
 * @return std::array<double, 2> The actual size of each grid cell.
 */
std::array<double, 2> calcRealCellSize(std::array<double, 2> cellSize, std::array<int, 2> gridSize){
    return {cellSize[0] * (gridSize[0]-1)/gridSize[0], cellSize[1] * (gridSize[1]-1)/gridSize[1]};
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



//detect fragments

/**
 * @brief Determines the category of a given direction.
 * 
 * @param dir The direction value from direction_raster.
 * @return int The category of the direction: 
 *         - 3 for diagonal (2, 8, 32, 128)
 *         - 1 for horizontal (1, 16)
 *         - 2 for vertical (4, 64)
 *         - Otherwise, returns 0.
 */
int getDirectionCategory(int dir){
    if(dir == 2 || dir == 8 || dir == 32 || dir == 128) return 3; //diagonal
    else if(dir == 1 || dir == 16) return 1; //horizontal
    else if(dir == 4 || dir == 64) return 2; //vertical
    return 0;
}

/**
 * @brief Computes the step vector for a given direction.
 * 
 * @param dir The direction encoded as an integer.
 * @return Dune::FieldVector<int, 2> A 2D vector representing the movement step for the given direction.
 */
Dune::FieldVector<int, 2> getStepForDir(int dir){
    if(dir ==0 || dir >128) std::cerr << "Invalid direction" << std::endl;
    Dune::FieldVector<int, 2> step;
    switch (dir)
        {
        case 1:     step = Dune::FieldVector<int, 2>{1, 0};
            break;
        case 2:     step = Dune::FieldVector<int, 2>{1,-1};
            break;
        case 4:     step = Dune::FieldVector<int, 2>{0, -1};
            break;
        case 8:     step = Dune::FieldVector<int, 2>{-1, -1};
            break;
        case 16:    step = Dune::FieldVector<int, 2>{-1, 0};
            break;
        case 32:    step = Dune::FieldVector<int, 2>{-1, 1};
            break;
        case 64:    step = Dune::FieldVector<int, 2>{0, 1};
            break;
        case 128:   step = Dune::FieldVector<int, 2>{1, 1};
            break;
        default:
            break;
        }
    return step;
}

/**
 * @brief Calculates the depth of a flow based on accumulation and a scaling factor.
 * 
 * @param avgAccumulation The average accumulation value.
 * @param scaleDephtFactor The factor used to scale the depth.
 * @return double The computed depth, constrained between 0.2 and 15 meters.
 */
double calcDepht(double avgAccumulation, double scaleDephtFactor){
    double depht = avgAccumulation/scaleDephtFactor;
    depht = std::min(15.0, depht); //max depht 15m
    depht = std::max(0.2, depht); //min depht 0.2m
    return depht;
}

/**
 * @brief Computes the width of the flow and the maximum number of iterations for refinement.
 * 
 * @param avgAccumulation The average accumulation value.
 * @param cellSize The size of a cell in the grid.
 * @param realCellSize The actual size of a cell.
 * @param fixedWidth Whether a fixed width should be applied.
 * @param directionCategory The category of the flow direction.
 * @param scaleWidthFactor The factor used to scale the width.
 * @param minWidth The minimum allowed width.
 * @param maxWidth The maximum allowed width if no value is provided it is set to cellSize..
 * @return std::pair<double, int> A pair containing the computed width and the maximum number of iterations.
 */
std::tuple<double, int, double> calcWidth(double avgAccumulation, std::array<double, 2> cellSize, std::array<double, 2> realCellSize, bool fixedWidth, int directionCategory,
                                double scaleWidthFactor, double minWidth,  double maxWidth){

    if (maxWidth < 0) maxWidth = (cellSize[0]+cellSize[1])/2;
    double width = avgAccumulation/(scaleWidthFactor); //width of flow in cells

    int maxIterationsFragment = 10000;
    double minSize = -1;


    if(fixedWidth){   
        if (avgAccumulation < 500){
            width = 9.6;
            maxIterationsFragment = 6 + int(directionCategory==3);
            minSize = 2;
        }
        else if (avgAccumulation < 2500) {
            width = 46.3;
            maxIterationsFragment = 2 + int(directionCategory==3)*3;
            minSize = 4;
            }
        else if (avgAccumulation < 5000){ 
            width = 65.5;
            maxIterationsFragment = 1 + int(directionCategory==1 || directionCategory==2)*3;
            minSize = 5;
        }
        else {
            width = 90;
            maxIterationsFragment = 4;
            minSize = 5;
        }
    }
    else{
        width = avgAccumulation/(scaleWidthFactor); 
        width = std::min(maxWidth, width); //max width realCellSize
        width = std::max(minWidth, width); //min width
    }
    
    width = width/cellSize[0] * realCellSize[0];
    return {width, maxIterationsFragment, minSize};
}

/**
 * @brief Removes unnecessary flow fragments from a river network.
 * 
 * @param rivers A vector of vector of flow fragments representing the river network (divided in bouniding boxes).
 * @param skip A vector indicating which cells should be skipped.
 * @param cellSize The size of a grid cell.
 * @param realCellSize The real-world size of a cell.
 * @param gridSize The size of the grid.
 * @return std::vector<std::vector<flowFragment>> The cleaned river network with skippable fragments removed.
 */
std::vector<std::vector<flowFragment>> deleteSkippableFragments (std::vector<std::vector<flowFragment>> rivers, std::vector<int> skip, 
                                                                std::array<double, 2> cellSize, std::array<double, 2> realCellSize, std::array<int, 2> gridSize){
    for(int i = 0; i<rivers.size(); i++){
        for (int j = rivers[i].size()-1; j>=0; j--){ //delete all fragments that can be skipped but were not skipped due to order
            auto start = rivers[i][j].start;
            Dune::FieldVector<double, 2> rasterStart = convertToRasterCoordinates(start, cellSize, realCellSize); 
            rasterStart = rasterStart - Dune::FieldVector<double, 2>{0.5, 0.5};
            if (skip[round(rasterStart[1])*gridSize[0]+round(rasterStart[0])]){ 
                rivers[i].erase(rivers[i].begin()+j);}
        }
    }
    return rivers;
}

/**
 * @brief Determines the bounding box in which a given point or set of points is located.
 * 
 * @param corners A vector of corner points.
 * @param xBisector The x-coordinate of the bisecting line.
 * @param yBisector The y-coordinate of the bisecting line.
 * @return int An integer representing the box index (from bottom to top, left to right).
 */
int selectBox( std::vector<Dune::FieldVector<double, 2>> corners, int xBisector, int yBisector){
    int box = 0;
    if(corners.size() ==1){
        if(corners[0][0] >= xBisector){
            box += 2;
        }
        if(corners[0][1] >= yBisector){
            box += 1;
        }
        return box;
    }
    for (auto c : corners){
        if(c[0] == xBisector) continue;
        if(c[0] > xBisector){
            box += 2;
        }
        break;
    }
    for (auto c : corners){
        if(c[1] == yBisector) continue;
        if(c[1] > yBisector){
            box += 1;
        }
        break;
    }
    return box;
}

std::vector<std::vector<flowFragment>> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, 
                            std::array<double, 2> cellSize, std::array<int, 2> gridSize, double minAcc, double maxAccDiff, bool fixedWidth, double scaleDephtFactor, 
                            double scaleWidthFactor, double minWidth, double maxWidth){
    int size = int(gridSize[0] * gridSize[1]);
    std::vector<int> skip(size, 0);

    int xBisector = std::floor(gridSize[0]/2);
    int yBisector = std::floor(gridSize[1]/2);

    std::vector<std::vector<flowFragment>> rivers(4, std::vector<flowFragment>{});

    const std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize);

    std::cout << "Start finding flow fragments" << std::endl;
    for (int i = 0; i < gridSize[0]; i++){ 
        for (int j = 0; j < gridSize[1]; j++){
            if(skip[j*gridSize[0]+i]) continue;
            if(accumulation_raster(i, j) > minAcc ){
                Dune::FieldVector<double, 2> startWorld = {-100, -100};
                Dune::FieldVector<double, 2> endWorld = {-100, -100};

                Dune::FieldVector<int, 2> start = {i, j};
                Dune::FieldVector<int, 2> end = start;

                int dir = int(direction_raster(i, j));
                if(dir>128 || dir==0) continue;
            
                Dune::FieldVector<int, 2> step = getStepForDir(dir);
                int directionCategory = getDirectionCategory(dir);
                int endDir = dir;

                double accStart = accumulation_raster(i, j);
                double accEnd = accStart;

                while (endDir == dir ){ //straight flow
                    skip[end[1]*gridSize[0]+end[0]]=1; //elements in middle of fragment don't need to be checked again
                    Dune::FieldVector<int, 2> currPoint = end;
                    end += step;
                    
                    if(end[0] >= gridSize[0] || end[1] >= gridSize[1] || end[0] < 0 || end[1] < 0) { //end is outside of grid
                        accEnd = accumulation_raster(currPoint[0], currPoint[1]); 
                        endWorld = Dune::FieldVector<double, 2>{currPoint+end}/2; //set end on border of grid
                        break;
                    }

                    endDir = int(direction_raster(end[0], end[1]));
                    accEnd = accumulation_raster(end[0], end[1]);

                    if(std::abs(accStart-accEnd)>maxAccDiff){
                        if(currPoint==start) break;
                        end = currPoint;
                        accEnd = accumulation_raster(end[0], end[1]); 
                        break;
                    }

                }

                skip[j*gridSize[0]+i] = 0; //start should not be skipped (skipped in first loop iteration)
                if(endWorld == Dune::FieldVector<double, 2>{-100, -100}) skip[end[1]*gridSize[0]+end[0]]=0; //end should not be skipped

                double avgAccumulation = (accStart+accEnd)/2;
                double depht = calcDepht(avgAccumulation, scaleDephtFactor);
                
                std::tuple<double, int, double> widthInfo = calcWidth(avgAccumulation, cellSize, realCellSize, fixedWidth, directionCategory, scaleWidthFactor, minWidth, maxWidth);
                double width =std::get<0>(widthInfo);
                int maxIterationsFragment = std::get<1>(widthInfo);
                double minSize = std::get<2>(widthInfo);

                Dune::FieldVector<double, 2> shift = {0.5, 0.5}; //fragments should start in middle of cell
                startWorld = Dune::FieldVector<double, 2>(start)+shift;
                if(endWorld == Dune::FieldVector<double, 2>{-100, -100}) endWorld = Dune::FieldVector<double, 2>(end)+shift;
                else endWorld += shift;
            
                startWorld = convertToGlobalCoordinates(startWorld, cellSize, realCellSize);
                endWorld = convertToGlobalCoordinates(endWorld, cellSize, realCellSize);
                //std::cout << width << std::endl;
                
                flowFragment f = {startWorld, endWorld, width, width, depht, maxIterationsFragment, directionCategory, minSize}; 

                int boxStart = selectBox({start}, xBisector, yBisector);
                int boxEnd = selectBox({end}, xBisector, yBisector);

                if(boxStart+boxEnd == 3){
                    for (int i = 0; i< rivers.size(); i++){
                        rivers[i].push_back(f);
                    }
                }
                else{
                    rivers[boxStart].push_back(f);
                    if(boxStart != boxEnd) rivers[boxEnd].push_back(f);
                }
                              
            }
        }
        
    }

    int count = 0;

    std::cout << "Number of fragments before skip: " << rivers[0].size() + rivers[1].size() + rivers[2].size() + rivers[3].size() << std::endl; //TODO delete
    rivers = deleteSkippableFragments(rivers, skip, cellSize, realCellSize, gridSize);
    std::cout << "Number of fragments: " << rivers[0].size() + rivers[1].size() + rivers[2].size() + rivers[3].size() << std::endl; //TODO delete

    std::vector<flowFragment> allFragments; //TODO delete
        for (auto& f : rivers){
            for (auto& ff : f){
              bool dublicate = false;
              for(auto& fff : allFragments){
                if(ff == fff){
                  dublicate = true;
                  break;
                }
              }
              if(!dublicate) allFragments.push_back(ff);
            }
        }

    std::cout << "Number of unique fragments: " << allFragments.size() << std::endl; //TODO delete
    
    for(auto& s : skip){
        if(s) count +=1;
    }

    std::cout << "Find fragments done. " << count+allFragments.size() << " cells with river." << std::endl; //TODO delete

    return rivers;
}


//refine

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
            f.end[0] -=0.01;
            f.end[1] -=0.01;
        }

        Dune::FieldVector<double, 2> flowVector = f.end - f.start;
        Dune::FieldVector<double, 2> normalFlowVector = calculateNormal(flowVector);

        if (f.widthEnd < 0) f.widthEnd = f.widthStart;

        //boundaries
        Dune::FieldVector<double, 2> start1 = f.start + (f.widthStart/2) * normalFlowVector;
        Dune::FieldVector<double, 2> start2 = f.start - (f.widthStart/2) * normalFlowVector;
        Dune::FieldVector<double, 2> end1 = f.end + (f.widthEnd/2) * normalFlowVector;
        Dune::FieldVector<double, 2> end2 = f.end - (f.widthEnd/2) * normalFlowVector;
        
        //boundary vectors
        Dune::FieldVector<double, 2> b1 = end1 - start1;
        Dune::FieldVector<double, 2> b2 = end2 - start2;

        //axis aligned bounding box for each fragment
        double minX = std::min({start1[0], start2[0], end1[0], end2[0]}); 
        double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
        double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
        double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});


        double minSize =  minSizeFactor*std::min(f.widthStart, f.widthEnd);
        if(f.minSize > 0) minSize = f.minSize;

        fragmentBoundaries fB = fragmentBoundaries{start1, start2, b1, b2, minX, maxX, minY, maxY, normalFlowVector, f.direction, minSize, f.depht, f.maxIterationsFragment};
        fragmentsBoundaries.push_back(fB);
    }
    return fragmentsBoundaries;
}

/**
 * @brief Computes the distance along a diagonal direction vector starting at a start point to a given point.
 * 
 * @param start The starting point of the vector.
 * @param vec The direction vector (needs to be diagonal).
 * @param point The target point to measure the distance to.
 * @return double The computed diagonal distance.
 */
double diagDistance(const Dune::FieldVector<double, 2> start, const Dune::FieldVector<double, 2> vec, const Dune::FieldVector<double, 2> point){
    if (abs(abs(vec[0]) - abs(vec[1])) > 1e-8) std::cerr << "Not a diagonal vector" << std::endl;
    const Dune::FieldVector<double, 2> startToPoint = point - start;
    double d = std::copysign(1, vec[0]) * startToPoint[0]- std::copysign(1, vec[1]) * startToPoint[1];
    double dist = d / sqrt(2);
    return abs(dist);
}

/**
 * @brief Calculates the hypotenuse vector and its center from a set of corner points and distances.
 * 
 * @param corners A vector containing three corner points.
 * @param dist01 The distance between corner 0 and corner 1.
 * @param dist02 The distance between corner 0 and corner 2.
 * @return std::array<Dune::FieldVector<double, 2>, 2> An array containing the hypotenuse vector and its center.
 */
std::array<Dune::FieldVector<double, 2>, 2> calcHypo(std::vector<Dune::FieldVector<double, 2>> corners, double dist01, double dist02){
    Dune::FieldVector<double, 2> hypo;
    Dune::FieldVector<double, 2> hypoCenter;
    
    if (std::abs(dist01-dist02)< 1e-6){ //90deg at 0
        hypo = corners[2] - corners[1];
        hypoCenter = corners[1] + hypo/2;
    }
    else if ( dist01 > dist02){  //90deg at 2
        hypo = corners[0] - corners[1];
        hypoCenter = corners[1] + hypo/2;
    }
    else { //90deg at 1
        hypo = corners[0] - corners[2];
        hypoCenter = corners[2] + hypo/2;
    }
    return {hypo, hypoCenter};
}

void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid, 
        std::vector<std::vector<flowFragment>> fragments, std::array<int, 2> gridSize, std::array<double, 2> cellSize, int maxIterations, double minSizeFactor,
        bool exactCalc){
    
    std::cout << "Refining";
    typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize);

    Dune::FieldVector<double, 2> bisectors = convertToGlobalCoordinates({std::floor(gridSize[0]/2), std::floor(gridSize[0]/2)}, cellSize, realCellSize);
    double xBisector = bisectors[0];
    double yBisector = bisectors[1];

    std::vector<std::vector<fragmentBoundaries>> fragmentsBoundaries;
    for(auto fragment : fragments){
        fragmentsBoundaries.push_back(calcFragmentBoundaries(fragment, minSizeFactor));
    }
    
    int c=0;
    bool change = true;
    while( c<maxIterations && change){ 
        std::cout << ".";
        c++; 
        change = false;
         
        for (const auto& element : elements(gridView)){ //check for each element if its part of the border of a flow --> refine

            std::vector<Dune::FieldVector<double, 2>> corners;
            for(int i = 0; i < 3; i++){
                corners.push_back(element.geometry().corner(i));
            }

            
            int box = selectBox(corners, xBisector, yBisector);
            
            bool marked = false;

            for (const auto& fB : fragmentsBoundaries[box]){
                if(marked) break;
                if(c>fB.maxIterationsFragment) continue;

                bool between = false;
                bool onSide1 = false;
                bool onSide2 = false;

                bool continueFragment = false;
                int xOut = 0;
                int yOut = 0;

                
                for (int i = 0; i < 3; i++){
                    auto cornerI =corners[i];
                    if (cornerI[0] < fB.minX-realCellSize[0]-0.1 || cornerI[0] > fB.maxX+realCellSize[0]+0.1 || cornerI[1] < fB.minY-realCellSize[1]-0.1 ||
                            cornerI[1] > fB.maxY+realCellSize[1]+0.1){  //element too far away
                        continueFragment = true; 
                        break;
                    }
                    
                    if(cornerI[0] < fB.minX) xOut--; //test if (corner of) triangle is outside of desired range
                    else if (cornerI[0] > fB.maxX)xOut++;
                    if(cornerI[1] < fB.minY) yOut--;
                    else if (cornerI[1] > fB.maxY) yOut ++;
                 
                    if(fB.direction == 1){ //horizontal
                        if(cornerI[1] > fB.maxY) onSide1 = true;
                        else if (cornerI[1] < fB.minY) onSide2 = true;
                        else between = true;
                        continue;
                    }

                    if(fB.direction == 2){ //vertical
                        if(cornerI[0] > fB.maxX) onSide2 = true;
                        else if (cornerI[0] < fB.minX) onSide1 = true;
                        else between = true;
                        continue;
                    }
        
                    Dune::FieldVector<double, 2> vec1 = cornerI - fB.start1;
                    Dune::FieldVector<double, 2> vec2 = cornerI - fB.start2;
                    double cross1 = fB.b1[0] * vec1[1] - fB.b1[1] * vec1[0];
                    double cross2 = fB.b2[0] * vec2[1] - fB.b2[1] * vec2[0];
                    
                    if(cross1 > 0) onSide1 = true;
                    else if(cross2 < 0) onSide2 = true;
                    else between = true;
                }


                if(continueFragment || std::abs(xOut) == 3 || std::abs(yOut) == 3 || between + onSide1 + onSide2 <= 1){continue; } // all 3 corners outside x/y boundary or element too far away or flow not trough element

                if(onSide1+onSide2==2) { //whole flow inside triangle
                    grid->mark(1, element); 
                    marked = true;
                    change=true; 
                    continue;}

                double dist01 = squaredDistance(corners[0], corners[1]);
                double dist02 = squaredDistance(corners[0], corners[2]);

                if (std::min({dist01, dist02}) < fB.minSize*fB.minSize){ continue;} //small enough

                std::array<Dune::FieldVector<double, 2>, 2> hypoAndCenter = calcHypo(corners, dist01, dist02);
                Dune::FieldVector<double, 2> hypo = hypoAndCenter[0];
                Dune::FieldVector<double, 2> hypoCenter = hypoAndCenter[1];



                if(fB.direction == 1 || fB.direction == 2){ //flow (and border) horizontal or vertical
        
                    if (exactCalc && (std::abs(hypo[0]) < 1e-8 || std::abs(hypo[1]) < 1e-8)){ //hypotenuse vertical or horizontal
                        grid->mark(1, element);
                        marked = true;
                        change = true;
                        continue;
                    }
                    if(fB.direction == 1 && (std::abs(fB.start1[1] - hypoCenter[1]) < fB.minSize/2 || std::abs(fB.start2[1] - hypoCenter[1]) < fB.minSize/2)){continue;} //horizontal: border close to middle of triangle (in right angle)
                    else if(fB.direction == 2 && (std::abs(fB.start1[0] - hypoCenter[0]) < fB.minSize/2 || std::abs(fB.start2[0] - hypoCenter[0]) < fB.minSize/2)){continue;} //vertical: border close to middle of triangle (in right angle)
                }
                else if(fB.direction == 3){ //flow diagonal
                    if (exactCalc && (std::abs(hypo[0] + hypo[1]) < 1e-8 || std::abs(hypo[0] - hypo[1]) < 1e-8)){ //hypothenuse diagonal
                        grid->mark(1, element);
                        marked = true;
                        change = true;
                        continue;
                    }
                                        
                    if(abs(diagDistance(fB.start1, fB.b1, hypoCenter))<fB.minSize/2 || abs(diagDistance(fB.start2, fB.b1, hypoCenter))<fB.minSize/2){continue;}
                }
            
                if(squaredDistance(fB.start1, hypoCenter) < 625 || squaredDistance(fB.start1, hypoCenter) < 625) continue; //start and and don't need to be refined perfectly
                grid->mark(1,element);
                marked = true;
                change = true;
            }
        }
        
        grid->preAdapt();
        grid->adapt();  
        grid->postAdapt();            
         
    }
    std::cout << std::endl;
}


//overall height

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


//height rivers

/**
 * @brief Calculates the height of a vertex in a flow fragment for diagonal flows.
 * 
 * If the point lies outside of the diagonal cells where the fragment was detected, the height value is from elevation_raster is too high.
 * The value is therefore interpolated from the neighboring cells.
 * If the point lies on a diagonal cell, false is returned and the normal calculation is used.
 * 
 * @param corner The vertex to calculate the height for.
 * @param cellSize The size of the grid cells.
 * @param realCellSize The actual size of the grid cells.
 * @param gridSize The dimensions of the grid.
 * @param elevation_raster The raster dataset containing elevation values.
 * @param fB The boundaries of the flow fragment.
 * @param originalHeightVertex The original height of the vertex.
 * @return std::pair<bool, double> Returns if height value is calculated and the new height value.
 */
std::pair<bool, double>  diagSpecialCase(Dune::FieldVector<double, 2> corner, const std::array<double, 2>& cellSize, const std::array<double, 2>& realCellSize, 
                                        std::array<int, 2> gridSize, const RasterDataSet<float>& elevation_raster, 
                                        fragmentBoundaries fB, double originalHeightVertex){
    Dune::FieldVector<double, 2> rasterCorner = convertToRasterCoordinates(corner, cellSize, realCellSize);

    std::vector<int> neighbors= calcNeighbors(rasterCorner, gridSize);
    int neighborX1 = neighbors[0];
    int neighborX2 = neighbors[1];
    int neighborY1 = neighbors[2];
    int neighborY2 = neighbors[3];


    if(neighborX2==-2 && neighborY2 == -2){ // vertex in cell (of raster data)
        
        if(std::abs(fB.b1[0]-fB.b1[1])<1e-8){ //flow from bottom left to top right or vice versa
            Dune::FieldVector<double, 2> bottomLeftCellCorner = convertToGlobalCoordinates({neighborX1, neighborY1}, cellSize, realCellSize);
            
            if(!isPointInFlow(bottomLeftCellCorner, fB)){ //flow through corner of cell (bottom left corner not in flow)
                double x_inCell = rasterCorner[0] - neighborX1;
                double y_inCell = rasterCorner[1] - neighborY1;

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
                return std::pair<bool, double>{true, std::min(originalHeightVertex, double((elevation_raster(x1, y1)+elevation_raster(x2, y2))/2)) - fB.depht};
            }
        }
        else{//flow from bottom right to top left or vice versa
            Dune::FieldVector<double, 2> bottomRightCellCorner = convertToGlobalCoordinates({neighborX1+1, neighborY1}, cellSize, realCellSize);

            if(!isPointInFlow(bottomRightCellCorner, fB)){ //flow through corner of cell (bottom right corner not in flow)
                double x_inCell = rasterCorner[0]- neighborX1;
                double y_inCell = rasterCorner[1]- neighborY1;

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
                return std::pair<bool, double>{true, std::min(originalHeightVertex, double((elevation_raster(x1, y1)+elevation_raster(x2, y2))/2)) - fB.depht};
            }
        }
    } else if(neighborX2==-2 && neighborY2 != -2){ //vertex between two cells in y direction
        return std::pair<bool, double>{true, std::min(elevation_raster(neighborX1, neighborY1), elevation_raster(neighborX1, neighborY2)) - fB.depht};
    } else if(neighborX2 != -2 && neighborY2 == -2){ //vertex between two cells in x direction
        return std::pair<bool, double>{true, std::min(elevation_raster(neighborX1, neighborY1), elevation_raster(neighborX2, neighborY1)) - fB.depht};
    } else{ //vertex between four cells
        return std::pair<bool, double>{true, std::min((elevation_raster(neighborX1, neighborY1) + elevation_raster(neighborX2, neighborY2))/2, 
                                                                (elevation_raster(neighborX2, neighborY1) + elevation_raster(neighborX1, neighborY2))/2) - fB.depht};
    }
    return std::pair<bool, double>{false, 0.0};
}


std::vector<double> applyFlowHeightFragments(const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView& gridView, 
       std::vector<std::vector<flowFragment>> fragments, std::vector<double> height, RasterDataSet<float> elevation_raster, std::array<double, 2> cellSize, 
        std::array<int, 2> gridSize){

    std::vector<double> originalHeight = height;

    std::array<double, 2> realCellSize = calcRealCellSize(cellSize, gridSize);

    Dune::FieldVector<double, 2> bisectors = convertToGlobalCoordinates({std::floor(gridSize[0]/2), std::floor(gridSize[0]/2)}, cellSize, realCellSize);
    double xBisector = bisectors[0];
    double yBisector = bisectors[1];

    std::vector<std::vector<fragmentBoundaries>> fragmentsBoundaries;
    for (int i = 0; i < fragments.size(); i++){
        fragmentsBoundaries.push_back(calcFragmentBoundaries(fragments[i]));
    }

    for(auto& v : vertices(gridView)){ //check for each vertex if its part of a fragment
        auto corner = v.geometry().corner(0); //avoid recomputation

        int selectFB = selectBox({corner}, xBisector, yBisector);
 
        for (const auto& fB : fragmentsBoundaries[selectFB]){

            if (!isPointInFlow(corner, fB)) continue;
            int vertexIndex = gridView.indexSet().index(v);
            if(fB.direction == 3){ //diagonal special case
                std::pair<bool, double> dsc = diagSpecialCase(corner, cellSize, realCellSize, gridSize, elevation_raster, fB, originalHeight[vertexIndex]);
                if(dsc.first){
                    height[vertexIndex] = dsc.second;
                    continue;
                }
            }         
  
            height[vertexIndex] = std::min(originalHeight[vertexIndex] - fB.depht, height[vertexIndex]);
        }
    }
    return height;
}


//all

std::vector<double> addRiversToMap(std::shared_ptr<Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>> grid,
                                    std::array<double, 2> cellSize, std::array<int, 2> gridSize,
                                    RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, RasterDataSet<float> elevation_raster,
                                    double minAcc, double maxAccDiff, bool fixedWidth, double scaleDephtFactor, double scaleWidthFactor, double minSizeFactor, 
                                    int maxIterations, double minWidth, double maxWidth, bool exactCalc){

    const Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming>::LeafGridView gridView = grid->leafGridView();
    std::vector<std::vector<flowFragment>> fragments = detectFragments(accumulation_raster, direction_raster, cellSize, gridSize, minAcc, maxAccDiff, fixedWidth, 
                                                                        scaleDephtFactor, scaleWidthFactor, minWidth, maxWidth);
    refineGridwithFragments(grid, fragments, gridSize, cellSize, maxIterations, minSizeFactor, exactCalc);
    
    std::vector<double> height = overallHeight(gridView, elevation_raster, cellSize, gridSize);

    height = applyFlowHeightFragments(gridView, fragments, height, elevation_raster,cellSize, gridSize);
    return height;
}




