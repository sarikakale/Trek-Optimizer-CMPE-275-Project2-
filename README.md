# Trek-Optimizer-CMPE-275-Project2-
Sometimes, it becomes necessary to commute to some desired places in an efficient and optimized way to save time, energy and expenses. A good approach to this problem is navigating through all the possible paths and wisely arrive to a solution. Shortest path algorithms are savior in this case.Many real time problems associated to shortest path algorithms that we face in day to day life are resolved. Algorithms like travelling salesman problem, ant colony optimization focus on solving these problems in an efficient and optimized way. We have come across such a real time problem where our solution provides a hiker with an optimized and best possible path to reach the destination.
#Algorithm 
Sometimes, it becomes necessary to commute to some desired places in an efficient and optimized way to save time, energy and expenses. A good approach to this problem is navigating through all the possible paths and wisely arrive to a solution. Shortest path algorithms are savior in this case.Many real time problems associated to shortest path algorithms that we face in day to day life are resolved. Algorithms like travelling salesman problem, ant colony optimization focus on solving these problems in an efficient and optimized way. We have come across such a real time problem where our solution provides a hiker with an optimized and best possible path to reach the destination.
In our scenario, ants are analogous to trekkers, pheromone trails are analogous to footprints and the visited cities are analogous to a hiking spot.
The ACO algorithm in our application works as follows : 

1.	Receive the number of hiking spots to be trekked from the user - (dimension).

2.	Receive the number of threads that will be used to process the application (threadnum).

3.	Generate a distance matrix of dimensions (dimension * dimension). We have initialized a 1-D array of size (dimension * dimension) so that there is no slowdown due to accessing vectors or 2-D arrays at non contiguous blocks. Matrix data generation is done by iterating within an OpenMP parallel for block.
4.	The TrekOptimizer class initializes the settings for the algorithm calculation and takes the generated distance graph and a footprint array which stores the strength of the footprint trail from 1 node to another.

5.	The TrekOptimizer initializes and generates the footprints 1-D contiguous matrix in run() as follows : 
For i == j values, the footprint is initialized to 0.0
For all other values, the initial footprint is set to the init_footprint setting (the default value is 0.5)

6.	We then iterate over the maxiters value, which by default it set to 100 in the settings array. 

7.	For each iteration, a trekkers[] array is created and the best path for each of the trekkers is calculated using the probability function prob().

8.	The prob() function calculates the probability of a trekker for travelling from 1 spot to another based on the footprint trail to that spot and also based on some randomization. The OpenMP reduction directive 

9.	The tourLen attribute that computes the tour length of travelling to all the hiking spots in the network, is calculated.

10.	Based on the best cost that is received from all the trekkers, bestCost and bestPath for the current iteration is calculated.
