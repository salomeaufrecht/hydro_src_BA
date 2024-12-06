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
#include <array>


constexpr int dim = 2; // Dimension des Grids

// Struktur zur Repräsentation eines Flussfragments
struct flowFragment
{
    Dune::FieldVector<double, dim> start;
    Dune::FieldVector<double, dim> end;
    double widthStart;
    double widthEnd = -1;
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
                       double minSizeFactor = 0.4);

std::vector<double> applyFlowHeightFragments(std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, flowFragment f, std::vector<double> height);

std::vector<double> applyFlowHeight (std::shared_ptr<Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming>> grid, 
        Dune::FieldVector<double, dim> start, Dune::FieldVector<double, dim> end, double widthStart, double widthEnd, std::vector<double> height);

#endif // FLOWFUNCTIONS_HH
