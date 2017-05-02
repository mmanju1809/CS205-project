# CS 205 Final Project: Mapping out Trajectories of Charged Defects

Click [here](https://rkuade.github.io/) to visit the website.  

## Introduction
The promise of the nitrogen-vacancy (NV) center in diamond as a system for implementing memory storage for a quantum computer has spawned interest in defect centers in related materials such as SiC. SiC is particularly interesting as it is a polymorphic material, exhibiting about 250 known polytypes, which imbues it with a degree of freedom unavailable in diamond. The three most common polytypes, 4H- and 6H-SiC and 3C-SiC, all have spin relaxation times ranging from 8 to 24 ms at 20 K (with 4H-SiC being the highest) and coherence persists up to room temperature [1]. In addition to the long spin coherence, a key feature is the ability to optically address (write in and read out) the spin states. However, much of the luminescence or emission of the defects is diverted into transitions involving scattering processes (and is not purely from the desired spin transitions) at ambient temperatures. Indeed, only about 4% of luminescence is from the desired transitions [2]. A potential solution is to place the defects near cavities on resonance with the desired transitions [3]. Positioning defects, however, is a non-trivial endeavor. The defects would be created at roughly the desired location using the process of focused ion beam implantation, but this process creates a lot of damage. In order to heal the damage the sample is annealed, causing the defects to diffuse and some to be lost through conversion to other species. The purpose of this study is then to assess the probability that the negatively charged silicon vacancy defect in 4H-SiC would be optimally positioned and exist given a certain initial position and a certain number of time steps.

## Background
Work has already been done on obtaining barriers to diffusion and the vibrational frequencies that provide the directionality and time scale for the motion [4,5], but these studies lack a comprehensive understanding of the potential diffusion pathways and as a result do not perform the kinetic Monte Carlo simulations necessary to truly establish diffusion pathways and probabilities with the Coulomb interaction.

## Method
There are two main stages to the calculation of the trajectory probability maps. In the first stage, density functional theory calculations will be carried out to obtain the barriers to diffusion for the various pathways and in the second stage kinetic Monte Carlo simulations will be carried out to map out trajectories. We have existing code for the second part in Matlab, which will be converted to C code. The kinetic Monte Carlo simulations will be carried out with Coulomb interaction and for 1D and 2D random walk test cases (to be compared with theoretical calculations). Without the Coulomb interaction the probability map will be calculated using a breadth first search down the tree of possible transitions (which is O(log(N)), where N is the number of possible pathways).

## Parallelization of BFS
We initially used a recursive crawling method for performing the breadth-first search (BFS) required to map probabilities in the lattice. This method involved traversing to a particular position in the lattice and directly running our BFS function from that position (a recursive call). We found that applying OpenMP directives to the recursive portion was difficult, especially since the probability that a defect appeared in a particular position on the lattice likely depended on other probabilities which may not have been updated.

Instead, we aimed to change the model of how we traverse through the 12-ary tree. We were considering switching entirely to the matrix-vector multiplication method of BFS with the boolean semiring model, but instead we opted for a serial traversal of BFS that still adopted the semiring model. We use two different semirings - one for updating the probability and one for updating the position and time steps. This semiring model is described in more detail in the advanced features section. This made it easier for us to directly write parallelizing pragmas because our serial implementation resulted in several for loops which are simple to parallelize.

## SPMD shared memory with OpenMP
We used the SPMD paradigm of parallelism in order to achieve initial speedup. We used the standard library for shared memory multiprocessing, OpenMP. When parallelizing for loops, we used the standard compiler directive pragma omp parallel, designating shared matrices and iterators appropriately with the shared and private keywords.

Initially, we tried applying OpenMP pragmas to every easily parallelizable section of the code. We realized this created too much parallel overhead, since a lot of the operations we aimed to parallelize involved 3x3 matrices (matrix-matrix multiplication, vector-matrix multiplication, etc.), and the cost of starting up threads comparing to the speedup from distributing the work simply was not worth it.

With the serial implementation of BFS in place, we were simply able to apply an OpenMP pragma to a loop in order to quickly iterate through the particles in the lattice. Additionally, we were able to quickly iterate through printing the probabilities of transitioning to positions by parallelizing that portion (this was especially simple since we don't care about the order in which we print things after they are calculated!).

## Hybrid SPMD with OpenMP MPI
We used the common hybrid model of combining message passing with MPI with OpenMP. By splitting the work that needs to be done across several nodes, we can achieve further speedup.

In our particular usage of MPI, we use MPI_Isend calls with the swtcher variable to communicate to a node that calculations corresponding to a particular transition (as determined by the swtcher variable) needs to be done.

Once all of the calculations are done, if we are in the root process, we print out the probabilities of transitioning to a particular location.

## Advanced features
### Semiring model
One of the advanced features we used was to convert our initially recursive implementation of BFS to a model similar to the semiring model. The initial, immediate advantage to this was that we could convert a recursive implementation into a serial one, which is an advantage in the ease of parallelization.

A disclaimer is that we cannot use the pure boolean semiring with 0, 1 elements and addition as binary OR and multiplication as binary AND. The reason for this is that the boolean semiring is used for determining whether or not we are able to traverse to a particular node in BFS. We care about more than that - we care about the probability of transition to a certain point, the coordinates associated with those points, and the time steps required to reach those points. Because of this, using a pure semiring is not feasible.

However, there are several important takeaways from the boolean semiring model that we can use to construct our own implementation of BFS. One of them is that we can use a semiring - the semiring of the non-negative reals, with addition and multiplication as usual which will generally cover the cases we require for our lattice traversal. Another very important idea is the construction of a graph as a matrix. Initially, we were considering tree traversal in a context conducive to human understanding - traversing down a branch of the tree. But by converting our tree to a matrix representation, we were able to gain insight on how to construct a serial implementation of BFS and therefore parallelize this operation.

### Kinetic Monte Carlo
The naive BFS traversal does not take into account the particle-particle Coulomb interactions that may be occurring as a defect moves through the lattice. In order to take into account this process, we use a simulation called kinetic Monte Carlo, which simulates some process with known transition rates as evolving through time. We need to map out all transitions and obtain their rates. Then, we construct a cumulative vector of transition rates. We choose a random number between 0 and 1 and multiply by the sum of the transition rates, and choose the transition corresponding to the element of the cumulative vector mapped to by this multiplication. Then we update time and restart the process at the new position.

1. A. L. Falk, B. B. Buckley, G. Calusine, W. F. Koehl, V. V. Dobrovitski, A. Politi, C. A. Zorman, P. X. L. Feng, and D. D. Awschalom, “Polytype control of spin qubits in silicon carbide,” Nature Communications, vol. 4, p. 1819, 05 2013.

2. I. Aharonovich, S. Castelletto, D. A. Simpson, C.-H. Su, A. D. Greentree, and S. Prawer, “Diamond-based single-photon emitters,” Reports on Progress in Physics, vol. 74, no. 7, p. 076501, 2011.

3. D. O. Bracher, X. Zhang, and E. L. Hu, “Selective purcell enhancement of two closely linked zero-phonon transitions of a silicon carbide color center,” arXiv:1609.03918, 2016.

4. E. Rauls, T. Frauenheim, A. Gali, and P. Deák, “Theoretical study of vacancy diffusion and vacancy-assisted clustering of antisites in sic,” Phys. Rev. B, vol. 68, p. 155208, Oct 2003.

5. X. Wang, M. Zhao, H. Bu, H. Zhang, X. He, and A. Wang, “Formation and annealing behaviors of qubit centers in 4h-sic from first principles,” Journal of Applied Physics, vol. 114, no. 19, p. 194305, 2013.
