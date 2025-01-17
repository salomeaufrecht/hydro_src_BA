#include "flowFunctions.hh"


/**
 * @brief calculated the squared distance of two points
 * 
 * @param p1 
 * @param p2 
 * @return double 
 */
double squaredDistance(Dune::FieldVector<double, dim> p1, Dune::FieldVector<double, dim> p2){
    Dune::FieldVector<double, dim> diff = (p1 - p2);
    return diff[0]*diff[0] + diff[1]*diff[1];
}


std::vector<fragmentBoundaries> calcFragmentBoundaries(std::vector<flowFragment> fragments, double minSizeFactor=0.4){
    
    std::vector<fragmentBoundaries> fragmentsBoundaries;

        for (auto f : fragments){ //store relevant information for each fragment so that it needs to be computed only once
        if(f.end==f.start){
            std::cout << "flow without end" << std::endl; //debug message
            f.end[0] -=0.001;
            f.end[1] -=0.001;
        }

        Dune::FieldVector<double, dim> flowVector = f.end - f.start;
        Dune::FieldVector<double, dim> normal_flowVector = { -flowVector[1]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1]), 
                                                        flowVector[0]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1])} ;

        if (f.widthEnd < 0) f.widthEnd = f.widthStart;

        Dune::FieldVector<double, dim> start1 = f.start + (f.widthStart/2) * normal_flowVector;
        Dune::FieldVector<double, dim> start2 = f.start - (f.widthStart/2) * normal_flowVector;
        Dune::FieldVector<double, dim> end1 = f.end + (f.widthEnd/2) * normal_flowVector;
        Dune::FieldVector<double, dim> end2 = f.end - (f.widthEnd/2) * normal_flowVector;

        Dune::FieldVector<double, dim> b1 = end1 - start1;
        Dune::FieldVector<double, dim> b2 = end2 - start2;

        int direction = 0;
        direction = (std::abs(flowVector[1]) < 1e-8 && f.widthStart==f.widthEnd) ? 1 : 0;
        direction = (std::abs(flowVector[0]) < 1e-8 && f.widthStart==f.widthEnd) ? 2 : 0;
        direction = ((std::abs(flowVector[0] + flowVector[1]) < 1e-8 || std::abs(flowVector[0] - flowVector[1]) < 1e-8) && f.widthStart==f.widthEnd) ? 3 : 0;


        double minX = std::min({start1[0], start2[0], end1[0], end2[0]}); //axis aligned bounding box for each fragment
        double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
        double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
        double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});
        
        double minSize =  minSizeFactor*std::min(f.widthStart, f.widthEnd);

        fragmentBoundaries fB = fragmentBoundaries{start1, start2, b1, b2, minX, maxX, minY, maxY, direction, minSize};
        fB.depht = f.depht;
        fragmentsBoundaries.push_back(fB);
    }

    return fragmentsBoundaries;
}

std::vector<double> overallHeight(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                std::vector<double> height, RasterDataSet<float> map, std::array<double, 2> cellSize){

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    for(auto& v : vertices(gridView)){
        auto cornerI = v.geometry().corner(0);
        double map_val =map(int(cornerI[0]/cellSize[0]), int(cornerI[1]/cellSize[1])); 
        if(map_val < -5000 || map_val > 10000) map_val = 0; //wrong values

        height[gridView.indexSet().index(v)] += map_val; //add height value bc rivers are stored in height as -depht --> results in height-depht as value
    }
    return height;
}



std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, std::vector<double> height){

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    std::vector<fragmentBoundaries> fragmentsBoundaries = calcFragmentBoundaries(fragments);

    
    for(auto& v : vertices(gridView)){ //check for each vertex if its part of a fragment
        auto cornerI = v.geometry().corner(0); //avoid recomputation
 
        for (const auto& fB : fragmentsBoundaries){

            if(cornerI[0] < fB.minX) continue; //test if vertex is outside bounding box
            else if (cornerI[0] > fB.maxX)continue;
            if(cornerI[1] < fB.minY) continue;
            else if (cornerI[1] > fB.maxY) continue;

            // if flow is horizontal/vertical the bounding box describes it correctly --> following can be skipped
            if (fB.direction == 0 || fB.direction == 3){  //computes if vertex is inside a flow
                Dune::FieldVector<double, dim> vec1 = cornerI - fB.start1;
                Dune::FieldVector<double, dim> vec2 = cornerI - fB.start2;
                double cross1 = fB.b1[0] * vec1[1] - fB.b1[1] * vec1[0];
                double cross2 = fB.b2[0] * vec2[1] - fB.b2[1] * vec2[0];
                
                if(cross1 > 0) continue;
                else if(cross2 < 0) continue;
            }
            height[gridView.indexSet().index(v)] = -fB.depht;
        }
    }
    return height;
}

void refineGridwithFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor, std::array<double, 2> cellSize){
    
    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    std::vector<fragmentBoundaries> fragmentsBoundaries = calcFragmentBoundaries(fragments, minSizeFactor);

    double epsilonAngle = 1e-5; //TODO: till which angle?
   
    
    int c=0;
    bool change = true;
    while( c<100 && change){ 
        std::cout << c << std::endl;
        change = false;
        c++;  
        for (const auto& element : elements(gridView)){ //check for each element if its part of the border of a flow --> refine
            
            auto corner0 = element.geometry().corner(0);
            auto corner1 = element.geometry().corner(1);
            auto corner2 = element.geometry().corner(2); 
            std::vector<Dune::FieldVector<double, dim>> corners = {corner0, corner1, corner2};
    
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
                        if (abs(cornerI[1]-fB.maxY)< 1e-5 || abs(cornerI[1]-fB.minY)< 1e-5){ //corner on border TODO: coeresind doesnet work --> refine?
                            grid->mark(-1, element); 
                            continueFragment = true;
                            break;
                        }
                        continue;
                    }

                    if(fB.direction == 2){ //vertical
                        if(cornerI[0] > fB.maxX) out2 = true;
                        else if (cornerI[0] < fB.minX) out1 = true;
                        else between = true;
                        if (abs(cornerI[0]-fB.maxX)< 1e-5 || abs(cornerI[0]-fB.minX)< 1e-5){
                            grid->mark(-1, element); 
                            continueFragment = true;
                            break;
                        }
                        continue;
                    }
        
                    Dune::FieldVector<double, dim> vec1 = cornerI - fB.start1;
                    Dune::FieldVector<double, dim> vec2 = cornerI - fB.start2;
                    double cross1 = fB.b1[0] * vec1[1] - fB.b1[1] * vec1[0];
                    double cross2 = fB.b2[0] * vec2[1] - fB.b2[1] * vec2[0];
                
                    if (std::abs(cross1) < 1e-5 ||std::abs( cross2) < 1e-5){ //border on edge of triangle
                        grid->mark(-1, element); 
                        continueFragment = true;
                        break;
                    }
                    
                    if(cross1 > 0) out1 = true;
                    else if(cross2 < 0) out2 = true;
                    else between = true;
                }

                
                if(continueFragment || std::abs(xOut) == 3 || std::abs(yOut) == 3 || between + out1 + out2 <= 1){continue;} // all 3 corners outside x/y boundary/border on corner/border not within triangle
                if(out1+out2==2) { //whole flow inside triangle
                    grid->mark(1, element); 
                    change=true; 
                    continue;}

                double dist01 = squaredDistance(element.geometry().corner(0), element.geometry().corner(1));
                double dist02 = squaredDistance(element.geometry().corner(0), element.geometry().corner(2));

                
                if (std::min({dist01, dist02}) < fB.minSize*fB.minSize){ continue;} //small enough

                Dune::FieldVector<double, 2> hypo;
                double yCenterHeight;
                Dune::FieldVector<double, dim> hypoCenter;

                //calc hypothenuse v; y value middle of height of triangle; point on center of hypothenuse
                if (std::abs(dist01-dist02)< 1e-6){ //90deg at 0
                    hypo = element.geometry().corner(2) - element.geometry().corner(1);
                    hypoCenter = element.geometry().corner(1) + hypo/2;

                    if (element.geometry().corner(0)[1] == element.geometry().corner(1)[1]){
                        yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(1)[1]) / 2;
                    }
                    else yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(2)[1]) / 2;
                    }
                else if ( dist01 > dist02){  //90deg at 2
                    hypo = element.geometry().corner(0) - element.geometry().corner(1);
                    hypoCenter = element.geometry().corner(1) + hypo/2;
                    if(element.geometry().corner(2)[1] == element.geometry().corner(0)[1]) {
                        yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(0)[1]) / 2;}
                    else yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(1)[1]) / 2;
                    }
                else { //90deg at 1
                    hypo = element.geometry().corner(0) - element.geometry().corner(2);
                    hypoCenter = element.geometry().corner(2) + hypo/2;
                    if (element.geometry().corner(1)[1] == element.geometry().corner(2)[1]){
                        yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(2)[1]) / 2;}
                    else yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(0)[1]) / 2;
                    }
            

                if(fB.direction == 1 || fB.direction == 2){ //flow (and border) horizontal or vertical (dot product with {1,0}/{0,1})
                    if (std::abs(hypo[0]) < epsilonAngle || std::abs(hypo[1]) < epsilonAngle){ //hypothenuse vertical or horizontal
                        grid->mark(1, element);//hypothenuse vertical or horizontal
                        change = true;
                        continue;
                    }                
                    if(fB.direction == 1 && (std::abs(fB.start1[1] - hypoCenter[1]) < fB.minSize/2 || std::abs(fB.start2[1] - hypoCenter[1]) < fB.minSize/2)){ continue;} // border close to middle of triangle (in right angle)
                    else if(fB.direction == 2 && (std::abs(fB.start1[0] - hypoCenter[0]) < fB.minSize/2 || std::abs(fB.start2[0] - hypoCenter[0]) < fB.minSize/2)){ continue;} // border close to middle of triangle (in right angle)
                }
                else if(fB.direction == 3){ //flow diagonal
                    if (std::abs(hypo[0] + hypo[1]) < 1e-8 || std::abs(hypo[0] - hypo[1]) < 1e-8){ //hypothenuse diagonal
                        grid->mark(1, element);
                        change = true;
                        continue;
                    }                
                    Dune::FieldVector<double, dim> vec1_U = hypoCenter + Dune::FieldVector<double, dim> {-fB.minSize/2, fB.minSize/2} - fB.start1;
                    Dune::FieldVector<double, dim> vec1_D = hypoCenter - Dune::FieldVector<double, dim> {-fB.minSize/2, fB.minSize/2} - fB.start1;
                    Dune::FieldVector<double, dim> vec2_U = hypoCenter + Dune::FieldVector<double, dim> {-fB.minSize/2, fB.minSize/2} - fB.start2;
                    Dune::FieldVector<double, dim> vec2_D = hypoCenter - Dune::FieldVector<double, dim> {-fB.minSize/2, fB.minSize/2} - fB.start2;

                    double cross1_U = fB.b1[0] * vec1_U[1] - fB.b1[1] * vec1_U[0];
                    double cross1_D = fB.b1[0] * vec1_D[1] - fB.b1[1] * vec1_D[0];
                    double cross2_U = fB.b2[0] * vec2_U[1] - fB.b2[1] * vec2_U[0];
                    double cross2_D = fB.b2[0] * vec2_D[1] - fB.b2[1] * vec2_D[0];
    
                    if((cross1_U > 0 && cross1_D < 0) ||(cross2_U > 0 && cross2_D < 0)){ continue;} // points on different sides of border
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


std::vector<flowFragment> detectFragments(RasterDataSet<float> accumulation_raster, RasterDataSet<unsigned char> direction_raster, std::array<double, 2> cellSize, std::array<int, dim> gridSize){
    std::vector<flowFragment> rivers;
    int size = int(gridSize[0] * gridSize[1]);
    std::vector<int> skip(size, 0);

    std::cout << "Start finding flow fragments" << std::endl;
    for (int i = 0; i< gridSize[0]; i++){ 
        for (int j=0; j<gridSize[1]; j++){
            if(skip[j*gridSize[0]+i]) continue;
            if(accumulation_raster(i, j)>50){
                std::cout << std::endl;
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
                    //std::cout << "#";
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

                    endI = round(end[0]); 
                    endJ = round(end[1]);
                    endDir = int(direction_raster(endI, endJ));
                    accEnd = accumulation_raster(endI, endJ);
                    std::cout << end <<  " " << endDir << " , " << dir<< std::endl;

                    if(end[0] > gridSize[0] || end[1]>gridSize[1] || end[0] <0 || end[1]<0 || std::abs(accStart-accEnd)>200) { //end is outside of grid or accumulation changed too much --> cut fragment
                        if(currPoint==start &&  std::abs(accStart-accEnd)>200) break; //flow over two cells --> fragments cannot be splitted
                        end=currPoint; //undo changes
                        endI = round(end[0]); 
                        endJ = round(end[1]);
                        skip[endJ*gridSize[0]+endI]=0; //end should not be skipped
                        accEnd = accumulation_raster(endI, endJ); 
                        break;
                    }

                }

                skip[j*gridSize[0]+i] = 0; //start should not be skipped (skipped in first loop iteration)
                double volume = (accStart+accEnd)/2;
                double width = std::sqrt(volume/3000) * 3; //width:depht 3:1
                double depht = std::sqrt(volume/3000);
                width = std::min(0.6, width);
                depht = std::min(0.4, depht);

                start[0] = (start[0] + 1) * ((gridSize[0]-1.0)/gridSize[0]) * cellSize[0]; //cell centered date vs corner data (n cells on lenght n-1)
                start[1] = (start[1] + 1) * ((gridSize[1]-1.0)/gridSize[1]) * cellSize[1];
                end[0]   = (end[0]   + 1) * ((gridSize[0]-1.0)/gridSize[0]) * cellSize[0];
                end[1]   = (end[1]   + 1) * ((gridSize[1]-1.0)/gridSize[1]) * cellSize[1];

                flowFragment f = {start, end, width * cellSize[0]}; //TODO *cellSize[0]? or cellSize[1]
                f.depht = depht * cellSize[0];
                if(f.depht < 0.2) std::cout <<"shallow flow: " << f.depht << std::endl;

                rivers.push_back(f);
            }
        }
    }

    std::cout << "\n Find fragments done. " << std::endl;

    for (int i = rivers.size()-1; i>=0; i--){ //delete all fragments that can be skipped but were not skipped due to order
        auto start = rivers[i].start;
        int start_index_x = round(start[0]/(cellSize[0] * ((gridSize[0]-1.0)/gridSize[0])) -1);
        int start_index_y = round(start[1]/(cellSize[1] * ((gridSize[1]-1.0)/gridSize[1])) -1);
        if (skip[start_index_y*gridSize[0]+start_index_x]) rivers.erase(rivers.begin()+i);
    }
    return rivers;
}


std::vector<double> adjustFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f,
                                             std::vector<double> height){
    double cellSize = 90;
    if(f.widthStart>=1) return adjustFlowHeight(grid, f.start, f.end, f.widthStart/cellSize, f.widthEnd/cellSize, f.depht, height);
    else return adjustFlowHeight(grid, f.start, f.end, f.widthStart, f.widthEnd, f.depht, height);
}

std::vector<double> adjustFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                    Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, 
                                    double depht, std::vector<double> height){ //TODO: make more efficient, better
    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    Dune::FieldVector<double, dim> flowVector = end - start;
    Dune::FieldVector<double, dim> normal_flowVector = { -flowVector[1]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1]), 
                                                    flowVector[0]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1])} ;

    if (widthEnd < 0) widthEnd = widthStart;

    Dune::FieldVector<double, dim> start1 = start + (widthStart/2) * normal_flowVector;
    Dune::FieldVector<double, dim> start2 = start - (widthStart/2) * normal_flowVector;
    Dune::FieldVector<double, dim> end1 = end + (widthEnd/2) * normal_flowVector;
    Dune::FieldVector<double, dim> end2 = end - (widthEnd/2) * normal_flowVector;

    Dune::FieldVector<double, dim> b1 = end1 - start1;
    Dune::FieldVector<double, dim> b2 = end2 - start2;

    double minX = std::min({start1[0], start2[0], end1[0], end2[0]});
    double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
    double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
    double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});
   
   std::vector<double> heights;
   double count = 0;
    
    for(auto& v : vertices(gridView)){

        auto cornerI = v.geometry().corner(0);
        
        if(cornerI[0] < minX) continue; //test if vertex is outside of desired range
        else if (cornerI[0] > maxX)continue;
        if(cornerI[1] < minY) continue;
        else if (cornerI[1] > maxY) continue;

        Dune::FieldVector<double, dim> vec1 = cornerI - start1;
        Dune::FieldVector<double, dim> vec2 = cornerI - start2;
        double cross1 = b1[0] * vec1[1] - b1[1] * vec1[0];
        double cross2 = b2[0] * vec2[1] - b2[1] * vec2[0];
    
        
        if(cross1 > 0) continue;
        else if(cross2 < 0) continue;
        heights.push_back(height[gridView.indexSet().index(v)]);
        count ++;
    }

    double mean = accumulate(heights.begin(), heights.end(), 0)/count;
    double standardDeviation = 0;
    for(int i = 0; i < heights.size(); ++i) {
        standardDeviation += pow(heights[i] - mean, 2);
    }
    standardDeviation = standardDeviation/heights.size();
    //std::cout << "mean: " << mean << ", standardder: " << standardDeviation << std::endl;

    for(auto& v : vertices(gridView)){

        auto cornerI = v.geometry().corner(0);
        
        if(cornerI[0] < minX) continue; //test if vertex is outside of desired range
        else if (cornerI[0] > maxX)continue;
        if(cornerI[1] < minY) continue;
        else if (cornerI[1] > maxY) continue;

        Dune::FieldVector<double, dim> vec1 = cornerI - start1;
        Dune::FieldVector<double, dim> vec2 = cornerI - start2;
        double cross1 = b1[0] * vec1[1] - b1[1] * vec1[0];
        double cross2 = b2[0] * vec2[1] - b2[1] * vec2[0];
    
        
        if(cross1 > 0) continue;
        else if(cross2 < 0) continue;
        if (std::abs(height[gridView.indexSet().index(v)] - mean)>standardDeviation) height[gridView.indexSet().index(v)] = mean;
        //height[gridView.indexSet().index(v)] = -5;
    }
    
    return height;
}
