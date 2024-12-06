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

    //std::cout <<"start: " << start1[0] << "|" << start1[1]<< std::endl;
    //std::cout <<"end: " << end1[0] << "|" << end1[1]<< std::endl;



    double minSize =  minSizeFactor*std::min(widthStart, widthEnd);
    double epsilon3 = minSize/2;
    //std::cout << "minsize " << minSize << std::endl;

    double epsilonAngle = 1e-5; //TODO: till which angle?
    bool vOrHFlow = (std::abs(b1[0]) < epsilonAngle || std::abs(b1[1]) < epsilonAngle) && ((std::abs(b2[0]) < epsilonAngle || std::abs(b2[1]) < epsilonAngle));
   
    
    int c=0;
    bool change = true;
    while( c<100 && change){
        change = false;
        c++;  
        for (const auto& element : elements(gridView)){

            bool between = false;
            bool out1 = false;
            bool out2 = false;

            bool continueFor = false;
            int xOut = 0;
            int yOut = 0;
            //std::cout << "\n" ;
            for (int i = 0; i < 3; i++){
                auto cornerI = element.geometry().corner(i);
                //std::cout << cornerI[0] << "|" << cornerI[1] << std::endl;
                if(cornerI[0] < std::min({start1[0], start2[0], end1[0], end2[0]})) xOut--; //test if trangle is outside of desired range
                else if (cornerI[0] > std::max({start1[0], start2[0], end1[0], end2[0]}))xOut++;
                if(cornerI[1] < std::min({start1[1], start2[1], end1[1], end2[1]})) yOut--;
                else if (cornerI[1] > std::max({start1[1], start2[1], end1[1], end2[1]})) yOut ++;
    
                Dune::FieldVector<double, dim> vec1 = cornerI - start1;
                Dune::FieldVector<double, dim> vec2 = cornerI - start2;
                double cross1 = b1[0] * vec1[1] - b1[1] * vec1[0];
                double cross2 = b2[0] * vec2[1] - b2[1] * vec2[0];
            
                if (std::abs(cross1) < 1e-5 ||std::abs( cross2) < 1e-5){ //border on edge of triangle
                    grid->mark(-1, element); 
                    continueFor = true;
                    change = true;
                    break;
                }
                
                if(cross1 > 0) out1 = true;
                else if(cross2 < 0) out2 = true;
                else between = true;
            }

            
            if(std::abs(xOut) == 3 || std::abs(yOut) == 3 || continueFor || between + out1 + out2 <= 1) continue; // all 3 corners outside x/y boundry/border on corner/border not within triangle
            if(out1+out2==2) { //whole flow inside triangle
                grid->mark(1, element); 
                change=true; 
                continue;}

            double dist01 = squaredDistance(element.geometry().corner(0), element.geometry().corner(1));
            double dist02 = squaredDistance(element.geometry().corner(0), element.geometry().corner(2));

            
            if (std::min({dist01, dist02}) < minSize*minSize){std::cout << "small" <<std::endl; continue;} //small enough

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
           


            //TODO: till which angle? 
            if(vOrHFlow){ //flow (and border) horizontal or vertical (dot product with {1,0}/{0,1})
                if (std::abs(hypo[0]) < epsilonAngle || std::abs(hypo[1]) < epsilonAngle){ //hypothenuse vertical or horizontal
                    grid->mark(1, element);//hypothenuse vertical or horizontal
                    change = true;
                    continue;
                }                
                if(std::abs(start1[1] - yCenterHeight) < epsilon3 || std::abs(start2[1] - yCenterHeight) < epsilon3){std::cout << "-" << std::endl; continue;} // border close to middle of triangle (in right angle)
            }
            else if(std::abs(flowVector[0] + flowVector[1]) < 1e-8 || std::abs(flowVector[0] - flowVector[1]) < 1e-8){ //flow diagonal
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
 
                if((cross1_U > 0 && cross1_D < 0) ||(cross2_U > 0 && cross2_D < 0)){std::cout << "/" << std::endl; continue;} //epsilon points on different sides of border
            }
        
            change = true;
            grid->mark(1,element);
            //std::cout << "*";
        }
       
        grid->preAdapt();
        grid->adapt();  
        grid->postAdapt();            
    }  
}

void flowWithFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        flowFragment f, double minSizeFactor){
    if(f.widthStart>=1) flow(grid, f.start, f.end, f.widthStart/pixelSize, f.widthEnd/pixelSize, minSizeFactor);
    else flow(grid, f.start, f.end, f.widthStart, f.widthEnd, minSizeFactor);
}


std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f, std::vector<double> height){
    if(f.widthStart>=1) return applyFlowHeight(grid, f.start, f.end, f.widthStart/pixelSize, f.widthEnd/pixelSize, height);
    else return applyFlowHeight(grid, f.start, f.end, f.widthStart, f.widthEnd, height);
}

std::vector<double> applyFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, std::vector<double> height){
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


   
    for (const auto& element : elements(gridView)){

        bool between = false;

        for(auto& v : vertices(gridView)){
            

            auto cornerI = v.geometry().corner(0);
            if(cornerI[0] < std::min({start1[0], start2[0], end1[0], end2[0]})) continue; //test if trangle is outside of desired range
            else if (cornerI[0] > std::max({start1[0], start2[0], end1[0], end2[0]}))continue;
            if(cornerI[1] < std::min({start1[1], start2[1], end1[1], end2[1]})) continue;
            else if (cornerI[1] > std::max({start1[1], start2[1], end1[1], end2[1]})) continue;

            Dune::FieldVector<double, dim> vec1 = cornerI - start1;
            Dune::FieldVector<double, dim> vec2 = cornerI - start2;
            double cross1 = b1[0] * vec1[1] - b1[1] * vec1[0];
            double cross2 = b2[0] * vec2[1] - b2[1] * vec2[0];
        
            
            if(cross1 > 0) continue;
            else if(cross2 < 0) continue;
            height[gridView.indexSet().index(v)] = -1;
        }
    }
    return height;
}
