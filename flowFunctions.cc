#include "flowFunctions.hh"

const double pixelSize = 90;


double squaredDistance(Dune::FieldVector<double, dim> p1, Dune::FieldVector<double, dim> p2){
    Dune::FieldVector<double, dim> diff = (p1 - p2);
    return diff[0]*diff[0] + diff[1]*diff[1];
}



void flow(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd , 
        double minSizeFactor){



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

    bool horFlow = (std::abs(flowVector[1]) < 1e-8 && widthStart==widthEnd) ? true : false;
    bool vertFlow = (std::abs(flowVector[0]) < 1e-8 && widthStart==widthEnd) ? true : false;
    bool diagFlow = (std::abs(flowVector[0] + flowVector[1]) < 1e-8 || std::abs(flowVector[0] - flowVector[1]) < 1e-8);


    double minX = std::min({start1[0], start2[0], end1[0], end2[0]});
    double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
    double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
    double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});

    //std::cout <<"start: " << start1[0] << "|" << start1[1]<< std::endl;
    //std::cout <<"end: " << end1[0] << "|" << end1[1]<< std::endl;



    double minSize =  minSizeFactor*std::min(widthStart, widthEnd);
    double epsilon3 = minSize/2;
    //std::cout << "minsize " << minSize << std::endl;

    double epsilonAngle = 1e-5; //TODO: till which angle?
   
    
    int c=0;
    bool change = true;
    while( c<100 && change){
         change = false;
        c++;  
        for (const auto& element : elements(gridView)){

            bool between = false;
            bool out1 = false;
            bool out2 = false;

            bool continueElem = false;
            int xOut = 0;
            int yOut = 0;
            //std::cout << "\n" ;
            for (int i = 0; i < 3; i++){
                //std::cout << i << std::endl;
                auto cornerI = element.geometry().corner(i);
                //std::cout << cornerI << std::endl;
                if (cornerI[0] < minX-1 || cornerI[0] > maxX+1 || cornerI[1] < minY-1 || cornerI[1] > maxY+1){  continueElem = true; break;} //element too far away
                if(cornerI[0] < minX) xOut--; //test if (corner of) trangle is outside of desired range
                else if (cornerI[0] > maxX)xOut++;
                if(cornerI[1] < minY) yOut--;
                else if (cornerI[1] > maxY) yOut ++;

                if(horFlow){
                    if(cornerI[1] > maxY) out1 = true;
                    else if (cornerI[1] < minY) out2 = true;
                    else between = true;
                     if (abs(cornerI[1]-maxY)< 1e-5 || abs(cornerI[1]-minY)< 1e-5){
                        grid->mark(-1, element); 
                        continueElem = true;
                        //std::cout << "b-!" << std::endl;
                        break;
                    }
                    continue;
                }

                if(vertFlow){
                    if(cornerI[0] > maxX) out2 = true;
                    else if (cornerI[0] < minX) out1 = true;
                    else between = true;
                    if (abs(cornerI[0]-maxX)< 1e-5 || abs(cornerI[0]-minX)< 1e-5){
                        grid->mark(-1, element); 
                        continueElem = true;
                        //std::cout << "b!!" << c << std::endl;
                        break;
                    }
                    continue;
                }
    
                Dune::FieldVector<double, dim> vec1 = cornerI - start1;
                Dune::FieldVector<double, dim> vec2 = cornerI - start2;
                double cross1 = b1[0] * vec1[1] - b1[1] * vec1[0];
                double cross2 = b2[0] * vec2[1] - b2[1] * vec2[0];
            
                if (std::abs(cross1) < 1e-5 ||std::abs( cross2) < 1e-5){ //border on edge of triangle
                    grid->mark(-1, element); 
                    continueElem = true;
                    //change = true;
                    //std::cout << "on border"   << cornerI<< std::endl;
                    break;
                }

                //std::cout << out1 << between << out2 << std::endl;

                
                if(cross1 > 0) out1 = true;
                else if(cross2 < 0) out2 = true;
                else between = true;
                //std::cout << out1 << between << out2 << "cor" << std::endl;
            }
            //std::cout << "end elem \n " << std::endl;

            
            if(continueElem || std::abs(xOut) == 3 || std::abs(yOut) == 3 || between + out1 + out2 <= 1){continue;} // all 3 corners outside x/y boundry/border on corner/border not within triangle
            if(out1+out2==2) { //whole flow inside triangle
                grid->mark(1, element); 
                change=true; 
                continue;}

            double dist01 = squaredDistance(element.geometry().corner(0), element.geometry().corner(1));
            double dist02 = squaredDistance(element.geometry().corner(0), element.geometry().corner(2));

            
            if (std::min({dist01, dist02}) < minSize*minSize){ continue;} //small enough

            Dune::FieldVector<double, 2> hypo;
            double yCenterHeight;
            Dune::FieldVector<double, dim> hypoCenter;
            //calc hypothenuse v; y value middle of height of triangle; point on center of hypothenuse
            if (std::abs(dist01-dist02)< 1e-6){ //90deg at 0
                //std::cout << "1" << std::endl;
                hypo = element.geometry().corner(2) - element.geometry().corner(1);
                hypoCenter = element.geometry().corner(1) + hypo/2;

                if (element.geometry().corner(0)[1] == element.geometry().corner(1)[1]){
                    yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(1)[1]) / 2;
                }
                else yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(2)[1]) / 2;
                }
            else if ( dist01 > dist02){  //90deg at 2
                //std::cout << "2" << std::endl;
                hypo = element.geometry().corner(0) - element.geometry().corner(1);
                hypoCenter = element.geometry().corner(1) + hypo/2;
                if(element.geometry().corner(2)[1] == element.geometry().corner(0)[1]) {
                    yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(0)[1]) / 2;}
                else yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(1)[1]) / 2;
                }
            else { //90deg at 1
                //std::cout << "3" << std::endl;
                hypo = element.geometry().corner(0) - element.geometry().corner(2);
                hypoCenter = element.geometry().corner(2) + hypo/2;
                if (element.geometry().corner(1)[1] == element.geometry().corner(2)[1]){
                    yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(2)[1]) / 2;}
                else yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(0)[1]) / 2;
                }
           

            //std::cout << hypo << std::endl;
            //TODO: till which angle? 
            if(vertFlow || horFlow){ //flow (and border) horizontal or vertical (dot product with {1,0}/{0,1})
                if (std::abs(hypo[0]) < epsilonAngle || std::abs(hypo[1]) < epsilonAngle){ //hypothenuse vertical or horizontal
                    grid->mark(1, element);//hypothenuse vertical or horizontal
                    change = true;
                    continue;
                }                
                if(horFlow && (std::abs(start1[1] - hypoCenter[1]) < epsilon3 || std::abs(start2[1] - hypoCenter[1]) < epsilon3)){ continue;} // border close to middle of triangle (in right angle)
                else if(vertFlow && (std::abs(start1[0] - hypoCenter[0]) < epsilon3 || std::abs(start2[0] - hypoCenter[0]) < epsilon3)){ continue;} // border close to middle of triangle (in right angle)
            }
            else if(diagFlow){ //flow diagonal
                if (std::abs(hypo[0] + hypo[1]) < 1e-8 || std::abs(hypo[0] - hypo[1]) < 1e-8){ //hypothenuse diagonal
                    grid->mark(1, element);
                    change = true;
                    continue;
                }                
                Dune::FieldVector<double, dim> vec1_U = hypoCenter + Dune::FieldVector<double, dim> {-epsilon3, epsilon3} - start1;
                Dune::FieldVector<double, dim> vec1_D = hypoCenter - Dune::FieldVector<double, dim> {-epsilon3, epsilon3} - start1;
                Dune::FieldVector<double, dim> vec2_U = hypoCenter + Dune::FieldVector<double, dim> {-epsilon3, epsilon3} - start2;
                Dune::FieldVector<double, dim> vec2_D = hypoCenter - Dune::FieldVector<double, dim> {-epsilon3, epsilon3} - start2;

                double cross1_U = b1[0] * vec1_U[1] - b1[1] * vec1_U[0];
                double cross1_D = b1[0] * vec1_D[1] - b1[1] * vec1_D[0];
                double cross2_U = b2[0] * vec2_U[1] - b2[1] * vec2_U[0];
                double cross2_D = b2[0] * vec2_D[1] - b2[1] * vec2_D[0];
 
                if((cross1_U > 0 && cross1_D < 0) ||(cross2_U > 0 && cross2_D < 0)){ continue;} //epsilon points on different sides of border
            }
        
            change = true;
            grid->mark(1,element);
        }
       
        grid->preAdapt();
        grid->adapt();  
        grid->postAdapt();            
    }  
}

void flowWithFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        flowFragment f, double minSizeFactor, double pixelSize){
    if(f.widthStart>=1) flow(grid, f.start, f.end, f.widthStart/pixelSize, f.widthEnd/pixelSize, minSizeFactor);
    else flow(grid, f.start, f.end, f.widthStart, f.widthEnd, minSizeFactor);
}


std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f,
                                             std::vector<double> height){
    if(f.widthStart>=1) return applyFlowHeight(grid, f.start, f.end, f.widthStart/pixelSize, f.widthEnd/pixelSize, f.depht, height); //TODO take this out?? /adapt depht and min val in code
    else return applyFlowHeight(grid, f.start, f.end, f.widthStart, f.widthEnd, f.depht, height);
}

std::vector<double> applyFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                    Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, 
                                    double depht, std::vector<double> height){


    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    if(end==start){
        std::cout << "flow without end" << std::endl;
        end[0] -=0.001;
        end[1] -=0.001;
    }

    Dune::FieldVector<double, dim> flowVector = end - start;
    Dune::FieldVector<double, dim> normal_flowVector = { -flowVector[1]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1]), 
                                                    flowVector[0]/ std::sqrt(flowVector[0]*flowVector[0] + flowVector[1]*flowVector[1])} ;

    if(widthStart==0){
        std::cout << "width = 0" << std::endl;
        widthStart = 0.01;
    }
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
        if(cornerI[0]<2 && cornerI[1]<2) std::cout  << start << "-" << end << ", " << widthStart << "-||- " << start1 << " -- " << end1 << " -||- " << start2 << " -- " <<end2 << std::endl;
        height[gridView.indexSet().index(v)] = -depht;
    }
    
    return height;
}

std::vector<double> overallHeigth(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                std::vector<double> height, RasterDataSet<float> map){
    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();

    for(auto& v : vertices(gridView)){

        auto cornerI = v.geometry().corner(0);
        double map_val =map(int(cornerI[0]), int(cornerI[1])); 
        if(map_val < -5000 || map_val > 1000) map_val = 0;
        //map_val = std::min(map_val, 10000.0);
        //map_val = std::max(map_val, -3000.0);
        height[gridView.indexSet().index(v)] += map_val;
        
    }
    return height;
}

std::vector<double> adjustFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f,
                                             std::vector<double> height){
    if(f.widthStart>=1) return adjustFlowHeight(grid, f.start, f.end, f.widthStart/pixelSize, f.widthEnd/pixelSize, f.depht, height);
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



std::vector<double> applyFlowHeightFragments2(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, std::vector<double> height){

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    std::vector<fragmentBoundaries> fragmentsBoundaries;

    for (auto f : fragments){
        if(f.end==f.start){
            std::cout << "flow without end" << std::endl;
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


        double minX = std::min({start1[0], start2[0], end1[0], end2[0]});
        double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
        double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
        double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});

        

        fragmentBoundaries fB = fragmentBoundaries{start1, start2, b1, b2, minX, maxX, minY, maxY, direction};
        fB.depht = f.depht;

        fragmentsBoundaries.push_back(fB);
 
    }

   
    
    for(auto& v : vertices(gridView)){
        auto cornerI = v.geometry().corner(0);

 
        for (const auto& fB : fragmentsBoundaries){

    
            if(cornerI[0] < fB.minX) continue; //test if vertex is outside of desired range
            else if (cornerI[0] > fB.maxX)continue;
            if(cornerI[1] < fB.minY) continue;
            else if (cornerI[1] > fB.maxY) continue;

            if (fB.direction == 0 || fB.direction == 3){

                //TODO hor/vert flow
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

void flowWithFragments2(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor, double pixelSize){
    

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid->leafGridView();


    std::vector<fragmentBoundaries> fragmentsBoundaries;

    for (auto f : fragments){
        
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


        double minX = std::min({start1[0], start2[0], end1[0], end2[0]});
        double maxX = std::max({start1[0], start2[0], end1[0], end2[0]});
        double minY = std::min({start1[1], start2[1], end1[1], end2[1]});
        double maxY = std::max({start1[1], start2[1], end1[1], end2[1]});

        double minSize =  minSizeFactor*std::min(f.widthStart, f.widthEnd);

        fragmentBoundaries fB = fragmentBoundaries{start1, start2, b1, b2, minX, maxX, minY, maxY, direction, minSize};

        fragmentsBoundaries.push_back(fB);
 
    }


    double epsilonAngle = 1e-5; //TODO: till which angle?
   
    

    int c=0;
    bool change = true;
    while( c<100 && change){
        std::cout << c << std::endl;
        change = false;
        c++;  
        for (const auto& element : elements(gridView)){
            
            //auto corner0 = element.geometry().corner(0);
            //auto corner1 = element.geometry().corner(1);
            //auto corner2 = element.geometry().corner(2); 
            //std::vector<Dune::FieldVector<double, dim> corners = {corner0, corner1, corner2};
            for (const auto& fB : fragmentsBoundaries){

                bool between = false;
                bool out1 = false;
                bool out2 = false;

                bool continueFragment = false;
                int xOut = 0;
                int yOut = 0;
                //std::cout << "\n" ;
                for (int i = 0; i < 3; i++){
                    //std::cout << i << std::endl;
                    auto cornerI = element.geometry().corner(i); //TODO replace
                    //std::cout << cornerI << std::endl;
                    if (cornerI[0] < fB.minX-1 || cornerI[0] > fB.maxX+1 || cornerI[1] < fB.minY-1 || cornerI[1] > fB.maxY+1){  continueFragment = true; break;} //element too far away
                    if(cornerI[0] < fB.minX) xOut--; //test if (corner of) trangle is outside of desired range
                    else if (cornerI[0] > fB.maxX)xOut++;
                    if(cornerI[1] < fB.minY) yOut--;
                    else if (cornerI[1] > fB.maxY) yOut ++;

                    if(fB.direction == 1){ //horizontal
                        if(cornerI[1] > fB.maxY) out1 = true;
                        else if (cornerI[1] < fB.minY) out2 = true;
                        else between = true;
                        if (abs(cornerI[1]-fB.maxY)< 1e-5 || abs(cornerI[1]-fB.minY)< 1e-5){
                            grid->mark(-1, element); 
                            continueFragment = true;
                            break;
                        }
                        continue;
                    }

                    if(fB.direction == 2){ //vertikal
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
                        //change = true;
                        //std::cout << "on border"   << cornerI<< std::endl;
                        break;
                    }

                    //std::cout << out1 << between << out2 << std::endl;

                    
                    if(cross1 > 0) out1 = true;
                    else if(cross2 < 0) out2 = true;
                    else between = true;
                    //std::cout << out1 << between << out2 << "cor" << std::endl;
                }
                //std::cout << "end elem \n " << std::endl;

                
                if(continueFragment || std::abs(xOut) == 3 || std::abs(yOut) == 3 || between + out1 + out2 <= 1){continue;} // all 3 corners outside x/y boundry/border on corner/border not within triangle
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
                    //std::cout << "1" << std::endl;
                    hypo = element.geometry().corner(2) - element.geometry().corner(1);
                    hypoCenter = element.geometry().corner(1) + hypo/2;

                    if (element.geometry().corner(0)[1] == element.geometry().corner(1)[1]){
                        yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(1)[1]) / 2;
                    }
                    else yCenterHeight = element.geometry().corner(0)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(2)[1]) / 2;
                    }
                else if ( dist01 > dist02){  //90deg at 2
                    //std::cout << "2" << std::endl;
                    hypo = element.geometry().corner(0) - element.geometry().corner(1);
                    hypoCenter = element.geometry().corner(1) + hypo/2;
                    if(element.geometry().corner(2)[1] == element.geometry().corner(0)[1]) {
                        yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(1)[1] - element.geometry().corner(0)[1]) / 2;}
                    else yCenterHeight = element.geometry().corner(2)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(1)[1]) / 2;
                    }
                else { //90deg at 1
                    //std::cout << "3" << std::endl;
                    hypo = element.geometry().corner(0) - element.geometry().corner(2);
                    hypoCenter = element.geometry().corner(2) + hypo/2;
                    if (element.geometry().corner(1)[1] == element.geometry().corner(2)[1]){
                        yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(0)[1] - element.geometry().corner(2)[1]) / 2;}
                    else yCenterHeight = element.geometry().corner(1)[1] + (element.geometry().corner(2)[1] - element.geometry().corner(0)[1]) / 2;
                    }
            

                //std::cout << hypo << std::endl;
                //TODO: till which angle? 
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
    
                    if((cross1_U > 0 && cross1_D < 0) ||(cross2_U > 0 && cross2_D < 0)){ continue;} //epsilon points on different sides of border
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

