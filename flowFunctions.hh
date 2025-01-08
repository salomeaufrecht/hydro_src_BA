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

// Struktur zur Repräsentation eines Flussfragments
struct flowFragment
{
    Dune::FieldVector<double, dim> start;
    Dune::FieldVector<double, dim> end;
    double widthStart;
    double widthEnd = -1;
    double depht = -2;
};

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

// Funktionen zur Berechnung von Distanzen
double squaredDistance(Dune::FieldVector<double, dim> p1, Dune::FieldVector<double, dim> p2);

// Hauptfunktionen für die Flussbearbeitung
void flow(std::shared_ptr<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>> grid,
          Dune::FieldVector<double, dim> start, 
          Dune::FieldVector<double, dim> end,
          double widthStart, 
          double widthEnd = -1,
          double minSizeFactor = 0.5);

void flowWithFragments(std::shared_ptr<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>> grid,
                       flowFragment f,
                       double minSizeFactor = 0.4,  double pixelSize = 90);

std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f, 
                                            std::vector<double> height);

std::vector<double> applyFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, double depht, std::vector<double> height);

std::vector<double> applyFlowHeightFragments2(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, std::vector<double> height);

std::vector<double> overallHeigth(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                std::vector<double> height, RasterDataSet<float> map, std::array<double, 2> H);


std::vector<double> adjustFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f,
                                             std::vector<double> height);

std::vector<double> adjustFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
                                    Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, 
                                    double depht, std::vector<double> height);

void flowWithFragments2(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        std::vector<flowFragment> fragments, double minSizeFactor, std::array<double, 2> pixelSize);

#endif // FLOWFUNCTIONS_HH
