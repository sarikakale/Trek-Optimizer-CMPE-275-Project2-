#include <vector>
#include <cmath>
#include <random>
#include <cassert>
#include <iostream>
#include <boost/chrono.hpp>
#include "omp.h"
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/utility.hpp>

/*
  Class that holds the implementation of the trekkers in the problem.
  They are used by the TrekOptimizer class to generate the
  bestCost and the bestPath
*/
class Trekker {
  int initSpot;
  double* footprints;
  double* dist;
  int dimension;
  std::vector<unsigned char> visited{};

public:
  /*
  Settings array that holds the initial settings of the problem.
  */
  std::vector<double> settings{1.0, 5.0, 0.5, 100};
  double tourLength;
  std::vector<int> builtPath{};

  int getInitSpot() const
  {
    return initSpot;
  }

  void setInitSpot(int initSpot)
  {
    Trekker::initSpot = initSpot;
    std::fill(std::begin(visited), std::end(visited), 0);
    visited[initSpot] = 1;
  }

  Trekker(int spot, double* footprints, double* distgraph, int dimension)
    : initSpot(spot), footprints(footprints), dist(distgraph), dimension(dimension)
  {
    visited.resize(dimension);
    visited[spot] = 1;
  }

  /*
    Function that generates a probability based on a formula using some randomized values.
  */
  double prob(int from, int to)
  {
    double nom = (std::pow(1 / dist[(from * dimension) + to], settings[1]) * std::pow(footprints[(from * dimension) + to], settings[0]));
    double denom = 0.0;

    //OpenMP reduction used
    #pragma omp parallel for reduction(+:denom)
    for (int i = 0; i < visited.size(); ++i) {
      if (visited[i] == 0) {
        denom += (std::pow(1 / dist[(from * dimension) + i], settings[1]) * std::pow(footprints[(from * dimension) + i], settings[0]));
      }
    }
    if (denom > 0.00001) {
      return nom / denom;
    }
    else {
      return 1.0;
    }
  }

  /*
  This method holds the core functionality of generating the tourLen of this node while travelling to all other nodes in the network.
  */
  double run()

  {
    int visitNum = 0;
    int size = dimension;
    int curSpot = initSpot;
    std::random_device rand;
    double tourLen = 0.0;
    visited.clear();
    visited.resize(size);
    visited[curSpot] = 1;
    builtPath.clear();
    builtPath.resize(size);
    builtPath[visitNum] = curSpot;
    while (visitNum < (size - 1)) {

      double probNum;
      double visitProb;
      int index;

      index = curSpot;
      do {
        ++index;
        if (index >= size) {
          index = 0;
        }
        if (visited[index]) {
          visitProb = 0;
        }
        else {
          visitProb = prob(curSpot, index);
        }
        probNum = static_cast<double>(rand()) / rand.max();
      }
      while (probNum > visitProb);
      tourLen += dist[(curSpot * dimension) + index];
      curSpot = index;
      visited[index] = 1;
      visitNum++;
      builtPath[visitNum] = curSpot;
    }
    tourLen += dist[(curSpot * dimension) + initSpot];
    tourLength = tourLen;
    return tourLen;
  }

};

/*
  Class that generates the best path and cost for the network
*/
class TrekOptimizer {
private:
  double* dist;
  double* footprints;
  int dimension;
  int threadnum;
  std::vector<int> bestPath{};
public:
  std::vector<double> settings{1.0, 5.0, 0.5, 100};
  int numTrekkers = -1;
  int maxIter;
  double bestCost;

  const std::vector<int>& getBestPath() const
  {
    return bestPath;
  }

  TrekOptimizer(int threadnum, double* distgraph, double* footprints, int iters = 10, int dimension = 4, int numTrekker = -1)
    : threadnum(threadnum), dist(distgraph), footprints(footprints), numTrekkers(numTrekker), dimension(dimension), maxIter(iters) { if (numTrekker == -1) numTrekkers = dimension;}

    /*
      Core function that performs calculation on the trekkers and updates the footprints[] array accordingly.
      The footprints[] array is weakened using 
    */
  double run()
  {
    using namespace boost::chrono;
    double size = dimension;
    assert(size > 0);
    double init_footprints = 1 / size;

    double init_trail = omp_get_wtime();

    //OpenMP for used
    #pragma omp parallel for num_threads(threadnum)
    for (int i = 0; i < dimension * dimension; ++i) {

      if (i % (dimension + 1) == 0)
        footprints[i] = 0.0;
      else
        footprints[i] = init_footprints;
    }
    double trail_intialize_time = omp_get_wtime() - init_trail;

    bestCost = std::numeric_limits<double>::max();



    double run_track_time = 0;
    double update_footprints_time = 0;

    for (int i = 0; i < maxIter; ++i) {
      std::vector<Trekker> trekkers(numTrekkers, {0, footprints, dist, dimension});
      int index = 0;

      // Run Trekkers
      double start_time = omp_get_wtime();

      //OpenMP for used
      #pragma omp parallel for num_threads(threadnum)
      for (int j = 0; j < trekkers.size(); j++) {
        Trekker trekker = trekkers[j];

        trekker.setInitSpot(index);
        trekker.settings = settings;
        double curCost = trekker.run();
        if (curCost <= bestCost) {
          bestCost = curCost;
          bestPath = trekker.builtPath;
        }
        if (++index >= size) index = 0;
        trekkers[j] = trekker;
      }
      double run_Trekkers = omp_get_wtime() - start_time;
      run_track_time = run_track_time + run_Trekkers;

      double start_time_foot = omp_get_wtime();

      //OpenMP for used
      #pragma omp parallel for num_threads(threadnum)
      for (int k = 0; k < trekkers.size(); k++) {

        Trekker trekker = trekkers[k];
        for (int l = 0; l < size; ++l) {
          int from = l;
          int to = l + 1;
          if ((l + 1) == size) {
            to = 0;
          }
          int trekkerFrom = trekker.builtPath[from];
          int trekkerTo = trekker.builtPath[to];
          //Footprints are updated/strengthened here
          footprints[(trekkerFrom * dimension) + trekkerTo] += (settings[2] * (settings[3] / trekker.tourLength));
          footprints[(trekkerTo * dimension) + trekkerFrom] = footprints[(trekkerFrom * dimension) + trekkerTo];


        }

      }
      double update_foot = omp_get_wtime() - start_time_foot;
      update_footprints_time = update_footprints_time + update_foot;

      // Footprints are obscured by the elements(wind, rain, snow) and weakened here
      #pragma omp parallel for num_threads(threadnum)
      for (int m = 0; m < (int)size ; m++) {
        for (int n = 0; n < (int)size; n++) {
          footprints[(m * dimension) + n] = footprints[(m * dimension) + n] * settings[2];
        }
      }
    }
    printf("TOTAL TIME IN TREK INITIALIZATION : %lf\n", trail_intialize_time);

    printf("TOTAL TIME IN RUN TREKKERS: %lf\n", run_track_time);
    printf("TOTAL TIME IN UPDATING FOOTPRINTS%lf\n", update_footprints_time );
    return bestCost;

  }


};

int trekPlanner() {

  using namespace boost::chrono;
  auto dt_s = high_resolution_clock::now();

  //TODO
  int dimension;
  printf("Enter the no. of spots to be covered : \n");
  scanf( "%d", &dimension);

  int threadnum;
  printf("Enter the no. of threads : \n");
  scanf( "%d", &threadnum);

  int    const NUM_ROWS = dimension;
  int    const NUM_COLS = dimension;
  size_t const N_BYTES = NUM_ROWS * NUM_COLS * sizeof(double);

  //Using 1-D arrays instead 
  double* distance   = (double*)malloc(N_BYTES);
  double* footprints   = (double*)malloc(N_BYTES);
  srand(time(NULL));
  //Populate 2d array for distance
  //OpenMP for used
  #pragma omp parallel for
  for (int i = 0; i < dimension; ++i)
  {
    for (int j = 0; j < dimension; ++j)
    {
      int index = i * dimension + j;
      if (i == j)
        distance[index] = 0;
      else
        distance[index] = index;

      std::cout <<  distance[index] << " ";
    }
    printf("\n");
  }


  TrekOptimizer trekOp(threadnum, distance, footprints, 100, dimension = dimension);

  double c = trekOp.run();
  printf("BestCost = %f\n", c);

  std::vector<int> p = trekOp.getBestPath();
  printf("Optimal path : \n");
  for (int i = 0; i < p.size(); i++) {
    printf("%d\t", p[i]);
  }


  //Freeing allocated memory
  free(distance);
  free(footprints);
  auto dt = duration_cast<nanoseconds> (high_resolution_clock::now() - dt_s);

  std::cout << "\n main() clocked at : " << "dt = " << dt.count() << " ns" << "\n";
  return 0;
}

/*
Boost Python module to drive the above C++ script.
*/
BOOST_PYTHON_MODULE(TrekOptimizer) {
  using namespace boost::python;

  // Expose the function hello().
  def("trekPlanner", trekPlanner);

  //Expose class TrekOptimizer
  class_<TrekOptimizer>("TrekOptimizer",
                        init<int, double* , double*, int, int , int >()).
  def("run", &TrekOptimizer::run);
}