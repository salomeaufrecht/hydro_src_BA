#include "config.h"
#include <iostream>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#include <iostream>
#include <array>
#include "flowFunctions.hh"

/*
for every triangle checks if corners are above, below or inbetween the edges --> if corners are on different sides border within trinagle

*/


// custom made grid function
template<typename GV>
class MyVTKFunction
    : public Dune::VTKFunction<GV>
{
  typedef typename GV::Grid::ctype DF;
  enum
  {
    n = GV::dimension
  };
  using Entity = typename GV::Grid::template Codim<0>::Entity;

public:
  MyVTKFunction(const GV& gv_, std::vector<double>& data_, std::string s_)
      : gv(gv_), data(data_), s(s_)
  {
  }

  virtual int ncomps() const override
  {
    return 1;
  }

  virtual double evaluate(int comp, const Entity &e, const Dune::FieldVector<DF, n> &local) const override
  {
    auto x = e.geometry().global(local);
    double minnorm = 1e10;
    int mincorner;
    for (int i=0; i<e.geometry().corners(); i++)
    {
      auto xc = e.geometry().corner(i);
      xc -= x;
      auto norm = xc.two_norm();
      if (norm<minnorm) { minnorm=norm; mincorner=i; }
    }
    return data[gv.indexSet().subIndex(e,mincorner,2)];
  }

  virtual std::string name() const override
  {
    return s;
  }

private:
  const GV& gv;
  std::vector<double>& data;
  std::string s;
};





using namespace Dune;



int main(int argc, char *argv[])
{
    // Set up MPI if available
    MPIHelper::instance(argc, argv);

    typedef Dune::ALUGrid< dim, dim, Dune::simplex, Dune::conforming > Grid;
    using GridView = Grid::LeafGridView;

    

    // Start with a structured grid
    const std::array<unsigned, dim> n = {6, 6};
    const FieldVector<double, dim> lower = {0, 0};
    const FieldVector<double, dim> upper = {6, 6};

    std::shared_ptr<Grid> grid = StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, n);

    const GridView gridView = grid->leafGridView();



    

    flowFragment f1 = {{3.5, 3.5}, {6, 3.5}, 9};
    flowFragment f2 = {{0, 0}, {3.5, 3.5}, 0.1, 0.1};
    flowFragment f3 = {{1.5, 1.5}, {1.5, 2}, 0.1};
    flowFragment f4 = {{1.5, 2}, {1.5, 3}, 0.2};
    flowFragment f5 = {{1.5, 3}, {1.5, 4}, 0.3};
    flowFragment f6 = {{1.5, 4}, {1.5, 5}, 0.4};
    flowFragment f7 = {{1.5, 5}, {1.5, 6}, 0.5};
    flowFragment f8 = {{1.5, 5}, {3.5, 5}, 0.1};

    
    std::vector<flowFragment> fragments1 = {f1, f2, f3, f4, f5, f6, f7,  f8};

    std::vector<flowFragment> fragments2 = {f8, f7, f6, f5, f4, f3, f2, f1};

    for(auto f : fragments2){
        flowWithFragments(grid, f);
    }

    std::vector<double> height(gridView.indexSet().size(2));
    for(auto f : fragments2){
        height = applyFlowHeightFragments(grid, f, height);
    }

    
    int subsampling = 2;
    Dune::SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementIntervals(subsampling));
     auto f = std::make_shared<MyVTKFunction<GridView>>(gridView,height,"test");
    // Write grid to file
    //VTKWriter<GridView> vtkWriter(gridView);
    //vtkWriter.addVertexData(height, "height");
    vtkWriter.addCellData(f);
    vtkWriter.write("refined_grid_with_height");
    
}
