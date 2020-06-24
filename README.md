# QuantumManyBody
 A Mathematica package to perform numerical computations on quantum many body physics problems.

 In recent days, an increasing number of physicists with a high-energy background (or, more specifically, with a background in string theory) are getting interested in condensed matter problems and quantum computing. In the string theory community, Wolfram Mathematica is the standard programming language that people use. This package hopes to be a useful tool for string theorists who are moving to condensed matter and quantum computing problems, by providing functionalities to perform numerical computations using their favorite language. In the near future, I hope to add more support to parallel computations, in order to make this package useful also for people who want to use Mathematica to perform high performance computations.
 Moreover, I hope to add more functionalities performing quantum algorithms.

 An auxiliary Mathematica notebook and a Jupyter notebook as well, containing some examples of usage, are included. These notebooks can be seen as tutorials, explaining the main functionalities provided by the package. Look at the [Jupyter notebook](https://github.com/Dario-Rosa85/QuantumManyBody/blob/master/QuantumManyBody_tutorial.ipynb) for an online preview of the current functionalities.

 Of course, the package is still in its first version. For any bugs or suggestions, you may contact me at [dario.rosa85@gmail.com](mailto:dario.rosa85@gmail.com).

 ## Main functionalities

 The main functions included so far in the package are

  - `SpinOperators`: this function defines the spin operators on a spin chain of given length.

  - `SpinChainHamiltonian`: this function receives a graph (or an hyper-graph), a corresponding set of coupling constants and provides a spin-chain Hamiltonian defined on the graph with the couplings assigned. It also gives a Suzuki-Trotter decomposition of the same Hamiltonian.

  - `FindGroundState`: this function, provided an Hamiltonian, computes the corresponding ground state.

  - `EnergyStored`: this function, provided a ket and an observable, computes the mean value of the observable in the given state.

  - `FindBandWidth`: this function, given an observable, computes its bandwidth, i.e. the difference between the maximum and the minimum eigenvalues.

  - `QITE`: this function implements a quantum version of the imaginary time evolution algorithm to find the ground state of a given Hamiltonian.

  - `EvolvedStateList`: this function computes the time evolution of a given initial state, as dictated by a given Hamiltonian.

  - `PartialTrace`: this function, given a certain pure state and a bi-partition of the spin chain on which the state is defined, computes the mixed state resulting by tracing out the degrees of freedom belonging to one of the subsystems.

  - `EntanglementEntropy`: this function computes the entanglement entropy associated to a given density matrix.

  - `GammaMatrices`: this function produces a set of Majorana operators.

  - `SYKHamiltonian`: this function generates the SYK Hamiltonian for a system of Majorana fermions.

 ## Revision history

 #### Version 0.1

  - Initial release.

 ## References

 The quantum imaginary time evolution (QITE) function, makes use of the results of

 #### [Motta, M. et al. (2019), "Determining eigenstates and thermal states on a quantum computer using quantum imaginary time evolution" , Nature Physics, 16(2), 205–210.](https://www.nature.com/articles/s41567-019-0704-4)
