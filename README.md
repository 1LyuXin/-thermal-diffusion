# -thermal-diffusion
There are three versions of this program developed in C/C++ for the solution of two-dimensional heat transfer problems
1. Serial version
2. MPICH blocking communication concurrent version
3. MPICH non-blocking communication concurrent version
Two-dimensional heat transfer is characterised by the fact that the temperature at a given coordinate at the current moment depends on the temperature at the four coordinates around the previous coordinate.
Therefore to calculate the data at the current location you need to know the temperature data around it at the previous time.
When writing with blocking communication, the calculation is done upon receipt of the data
When using non-blocking communication, it is more difficult than using blocking communication to take into account the dependencies between the data and to avoid confusion between the data.
By analysing the results of comparing the operation of the programs we can conclude that
1. for the concurrent version of MPI with the two-dimensional heat transfer problem, a better speed-up can be obtained with a reasonable number of processes.
2. For the concurrent version of MPICH, we can learn that the time taken to run the program decreases and then increases gradually as the number of processes increases, because the processes themselves take up system resources.   As the number of processes increases, the processes themselves take up more and more resources and each process is allocated fewer and fewer tasks, resulting in lower process utilization. It is therefore necessary to debug to find the number of processes for a given problem size.
3. A comparison of non-blocking and blocking communication shows that non-blocking communication takes less time on average than blocking communication for the same problem size and number of processes, because some of the time spent in blocking communication is waiting and not fully utilised.

