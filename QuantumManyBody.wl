(* ::Package:: *)

(* :Context: QuantumManyBody` *)

(* :Title: QuantumManyBody *)

(* :Author: Dario Rosa *)

(* :Version: Mathematica 12.0 *)

(* :Package Version: 1.0 *)

(* :Copyright: 

   Copyright 2020 QuantumManyBody (https://github.com/Dario-Rosa85/QuantumManyBody) 

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*)

(* :Keywords:
    QMB , Hamiltonians , PartialTraces, QITE
*)

(* :Summary:
This package implements several functionalities for numerical computations of quantum
many-body systems.
*)


BeginPackage["QuantumManyBody`"]

SpinOperators::usage = "SpinOperators[L] produces a set of spin operators, \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(a\)]\), for a spin chain of length L."

GammaMatrices::usage = "GammaMatrices[n] produces a set of 2 * n gamma matrices, \!\(\*SuperscriptBox[\(\[Gamma]\), \(i\)]\), satisfying the algebra {\!\(\*SuperscriptBox[\(\[Gamma]\), \(\(i\)\(\\\ \)\)]\), \!\(\*SuperscriptBox[\(\[Gamma]\), \(j\)]\)} = 2 \!\(\*SuperscriptBox[\(\[Delta]\), \(ij\)]\)."

EnergyStored::usage = "EnergyStored[x , hamiltonian] computes numerically the mean value of hamiltonian in the state x."

EvolvedStateList::usage = "EvolvedStateList[stateIn , hamiltonian , tMin , tMax , nPoints] computes the time evolution of stateIn, from tMin to tMax, with the hamiltonian and with nPoints intermediate 
steps, equally spaced. 
EvolvedStateList[stateIn , hamiltonian , tMin , tMax , nPoints , Log] computes the time evolution of stateIn, from tMin to tMax, with the hamiltonian and with nPoints intermediate 
steps, equally spaced in \!\(\*SubscriptBox[\(Log\), \(10\)]\) scale." 

FindGroundState::usage = "FindGroundState[hamiltonian] computes numerically the ground state of a given hamiltonian."

FindBandwidth::usage = "FindBandwidth[hamiltonian] computes numerically the bandwidth of a given hamiltonian."

PartialTrace::usage = "PartialTrace[stateIn , tracedSites_List , nSites] computes the partial trace of the state stateIn, defined over a spin chain of length nSites, over the lattice sites defined in the list tracedSites."

EntanglementEntropy::usage = "EntanglementEntropy[densityMatrix] computes the von Neumann entanglement entropy of the given densityMatrix."

XCouplingRegular::usage = "XCouplingRegular[nGamma , kK , qQ] creates an \!\(\*SubscriptBox[\(x\), \(\(\*SubscriptBox[\(i\), \(1\)] ... \) \*SubscriptBox[\(i\), \(qQ\)]\)]\)coupling, for a sparse SYK model defined on a (qQ * kK)-regular, qQ-hypergraph."

SYKHamiltonian::usage = "SYKHamiltonian[qQ , nGamma , gamma , xCoupling : {1} , seedRandom : -1] computes the \!\(\*SubscriptBox[\(SYK\), \(qQ\)]\) Hamiltonian for a system of nGamma Majorana fermions, represented by the gamma matrices gamma.
SYKHamiltonian[qQ , nGamma , gamma , xCoupling , seedRandom : -1] computes the \!\(\*SubscriptBox[\(SYK\), \(qQ\)]\) Hamiltonian for a system of nGamma Majorana fermions, represented by the gamma matrices gamma and with the non-vanishing couplings xCoupling turned on only."

SYKHamiltonianBoson::usage = "SYKHamiltonianBoson[qQ , nGamma , gamma] computes the \!\(\*SubscriptBox[\(SYK\), \(qQ\)]\) Hamiltonian for a system of nGamma real hard-core bosons, represented by the matrices gamma." 

SpinChainHamiltonian::usage = "SpinChainHamiltonian[listCouplings , listCoefficients] creates the spin chain Hamiltonian, coupling the spins according to listCouplings, with anisotropic couplings given by listCoefficients.
The output is given as {hamiltonianSpinChainUnitaries , hamiltonianSpinChain}. hamiltonianSpinChainUnitaries contains each term of the decomposition in operators proportional to unitaries {{\!\(\*SuperscriptBox[SubscriptBox[\(U\), \(1\)], \(x\)]\) , \!\(\*SuperscriptBox[SubscriptBox[\(U\), \(1\)], \(y\)]\) , \!\(\*SuperscriptBox[SubscriptBox[\(U\), \(1\)], \(z\)]\)}, ... , {\!\(\*SuperscriptBox[SubscriptBox[\(U\), \(L\)], \(x\)]\) , \!\(\*SuperscriptBox[SubscriptBox[\(U\), \(L\)], \(y\)]\) , \!\(\*SuperscriptBox[SubscriptBox[\(U\), \(L\)], \(z\)]\)}}, with L = Length @ listCouplings.
hamiltonianSpinChain is the total Hamiltonian: hamiltonianSpinChain = Total @ Flatten[hamiltonianSpinChainUnitaries , 1]."

QITE::usage = "QITE[latticeSize , hamiltonianDecomposed , listCouplings , \[Beta] , \[CapitalDelta]t , initialState , dD , toleranceNumeric :10^-6 , deltaDiagonal : 0.1] performs the quantum version of the imaginary time evolution of initialState, up to imaginary time \[Beta] (with a time step \[CapitalDelta]t) with a Trotterized hamiltonian given by hamiltonianDecomposed. 
Each term of the Trotterized Hamiltonian takes the form of hamiltonianDecomposed[[pos1 , pos2]] where pos1 denotes the sites of the spin chain involved and pos2 keeps into account of how many terms have the same support.
The list listCouplings takes into account all the connections among different sites of the underlying spin chain (we have Length @ listCouplings equal to Length @ hamiltonianDecomposed).
dD sets how larger the domain of the unitary approximation of non-unitary operator Exp[- \[CapitalDelta]t * hamiltonianDecomposed[[pos1 , pos2]]] should be."

Begin["`Private`"]

s1 = Symbol["s1"];
s2 = Symbol["s2"];
s3 = Symbol["s3"];
id = Symbol["id"];

toMat[list_] := KroneckerProduct @@ (list /. {s1 -> SparseArray @ PauliMatrix @ 1 , s2 -> SparseArray @ PauliMatrix @ 2 , s3 -> SparseArray @ PauliMatrix @ 3 , id -> SparseArray @ IdentityMatrix @ 2});

evolutionFunction[stateIn_ , hamiltonian_ , \[CapitalDelta]_] := Return[MatrixExp[- I * hamiltonian * \[CapitalDelta] , stateIn]]

stepEvolutionQITE[initialKet_ , hTrotterStep_ , sitesExtended_ , sS_ , identityOperator_ , toleranceNumeric_ , deltaDiagonal_ , \[CapitalDelta]t_]:= Block[{unitaryOperatorsExtended , \[CapitalDelta] , sSMatHalf , bB ,sSMat , aAOperator} , 
unitaryOperatorsExtended = Prepend[Flatten[(Dot @@@ Tuples @ sS[[#]] &)/@ sitesExtended , 1] , identityOperator];
\[CapitalDelta] = (- hTrotterStep + (Conjugate @ initialKet . hTrotterStep . initialKet) * identityOperator) . initialKet; (* expansion at the first order in \[CapitalDelta]t *)
{sSMatHalf , bB} = (Chop[{# . initialKet , -2 Im[Conjugate @ initialKet . Conjugate @ # . \[CapitalDelta]]} , toleranceNumeric] &)/@ unitaryOperatorsExtended // Transpose;
sSMat =  Chop[(Level[sSMatHalf , 1] . Conjugate @ # &)/@ sSMatHalf  , toleranceNumeric];
aAOperator = SparseArray[Chop[Total[LeastSquares[(sSMat + Transpose @ sSMat) + deltaDiagonal *  IdentityMatrix[Length @ bB] , - bB ] * unitaryOperatorsExtended] ,toleranceNumeric] , {Length @ initialKet , Length @ initialKet} , 0];
Return[Chop[ MatrixExp[-I * \[CapitalDelta]t * aAOperator , initialKet] , toleranceNumeric]]
]

partialTrace[latticeSizeLocal_ , tracedLocal_ , \[Psi]Local_]:=Block[{reducedNumbers = 2^(latticeSizeLocal - tracedLocal) , tracedHilbert = 2^tracedLocal , vectorReduced},
vectorReduced[x_] := vectorReduced[x] = Take[\[Psi]Local , {tracedHilbert(x-1) + 1 , tracedHilbert * x }];  
Return[Array[Total[vectorReduced[#1] * Conjugate @ vectorReduced[#2]] & , {reducedNumbers , reducedNumbers}]]
]

permutationSpins[startingIndices_ , finalIndices_] := Block[{tracedSitesStarting = Range @ Length @ finalIndices , rulePermutation} ,
rulePermutation = (#1 -> #2 &) @@@ Transpose @ {Complement[tracedSitesStarting , finalIndices] , Complement[finalIndices , tracedSitesStarting]};
Return[Select[Transpose @ {startingIndices , tracedSitesStarting /. rulePermutation} , #[[1]] != #[[2]] &]]
]

xCouplingMoving[nGamma_ , kK_ , qQ_ , listIndices_] :=  Block[
{indicesRemoved , elementsRemoved , candidateList , indicesChosen = listIndices[[1]] , indicesToChoose = listIndices[[2]]}, 
indicesRemoved = If[Length @ indicesToChoose >= qQ , RandomSample[Range @ Length @ indicesToChoose , qQ] , Return @ listIndices];
elementsRemoved = indicesToChoose[[indicesRemoved]] // Sort;
candidateList = {PadRight[indicesChosen , Length @ indicesChosen + 1 , {elementsRemoved}] , Delete[indicesToChoose , Partition[indicesRemoved , 1]]};
If[(DeleteDuplicates @ elementsRemoved == elementsRemoved) && (DeleteDuplicates @ candidateList == candidateList) , Return[candidateList] , Return[listIndices]]
]

parityFunc[i_] := (1 + ((-1)^i) )/2;

SpinOperators[latticeSize_] := Return[Partition[(toMat @ PadRight[PadLeft[{#2} , #1 , id] ,latticeSize , id]  &) @@@ Tuples[{Range @ latticeSize , {s1 , s2 , s3}}] , 3]]

GammaMatrices[n_] := Return[( toMat @ PadRight[PadLeft[If[OddQ @ # , {s1} , {s2}] , IntegerPart[(# - 0.1)/2] + 1 , s3] , n , id] &) /@ Range[2 * n]] 

EnergyStored[x_List , hamiltonian_] := Return[Conjugate[x] . hamiltonian . x]

FindGroundState[hamiltonian_] := Return[- Eigenvectors[- hamiltonian // N , 1 , Method -> {"Arnoldi" , "Criteria" -> "RealPart"}][[1]]]

FindBandwidth[hamiltonian_] := Return[Eigenvalues[ hamiltonian // N  , 1 , Method -> {"Arnoldi" , "Criteria" -> "RealPart"}][[1]] + Eigenvalues[- hamiltonian // N , 1 , Method -> {"Arnoldi" , "Criteria" -> "RealPart"}][[1]]]

EvolvedStateList[stateIn_ , hamiltonian_ , tMin_ , tMax_ , nPoints_ , spacingType_ : "Linear"] := Block[{\[Epsilon] = 10^-9 , points , \[Delta]t , evolvedState} , 
If[spacingType == "Log" , 
	(points = (10^# &) /@ N @ Range[Log[10 , Max[tMin , \[Epsilon]]] , Log[10 , tMax] , (Log[10 , tMax] - Log[10 , tMin]) / (nPoints - 1)];
	\[Delta]t = Differences[points];
	PrependTo[\[Delta]t , N @ tMin];) , 
	(points=N @ Range[tMin , tMax , (tMax - tMin) / (nPoints - 1)];
	\[Delta]t = Differences[points];
	PrependTo[\[Delta]t , N @ tMin];)
];
evolvedState = Delete[FoldList[(evolutionFunction[#1 , hamiltonian , #2] &) , stateIn , \[Delta]t] , 1];
Return[Transpose @ {points , evolvedState}]
]

PartialTrace[stateIn_ , tracedSites_List , nSites_] := Block[{listInitial , listDigit , nTraced , densityMatrix} , 
listInitial = Range[0 , 2^nSites - 1];
listDigit = IntegerDigits[# , 2 , nSites] & /@ listInitial;
nTraced = Length @ tracedSites;
Return @ partialTrace[nSites , nTraced , Permute[stateIn , FindPermutation[listInitial , (FromDigits[Permute[#, Cycles[permutationSpins[Range @ nTraced , tracedSites]]] , 2] &) /@ listDigit]]]
]

EntanglementEntropy[densityMatrix_] := Return[- Total[(# * Log[2 , #] &) /@ Select[Re @ Eigenvalues[densityMatrix // N] , Chop @ # != 0. &]]]

XCouplingRegular[nGamma_ , kK_ , qQ_ ] := Block[{xCouplingCandidate , listIndices = {{} , (Range @ nGamma &) /@ Range[kK * qQ] // Flatten}} , 
xCouplingCandidate = Nest[xCouplingMoving[nGamma , kK , qQ , #] & , listIndices , Round[3 * kK * nGamma]];
If[xCouplingCandidate[[2]] == {} , Return[xCouplingCandidate[[1]]] , XCouplingRegular[nGamma , kK , qQ]]
]

SYKHamiltonian[qQ_ , nGamma_ , gamma_ , xCoupling_ : {1} , seedRandom_ : -1] := Block[{indices , iterators , hamiltonianSYK , length = 2^(nGamma / 2)} ,
If[seedRandom == -1 , SeedRandom[] , SeedRandom[seedRandom]];
hamiltonianSYK = If[xCoupling == {1} ,
(
indices = Table[Symbol["i" <> ToString[i]] , {i , 1 , qQ}];
iterators = Transpose @ {indices , Range[qQ , 1 , -1] , Flatten[{nGamma , Take[indices ,  qQ -  1] -1}]}; 
(I/2)^(qQ / 2) * Sum[RandomVariate @ NormalDistribution[0. , Sqrt[Factorial[qQ - 1]/(nGamma^(qQ - 1))]] * Dot @@ gamma[[##]] & @ Evaluate @ Reverse @ indices  ,  ##] &  @@ iterators 
), 
SparseArray[Fold[#1 +  (I/2)^(qQ / 2) * RandomVariate @ NormalDistribution[0. , Sqrt[Factorial[qQ - 1]/((Length @ xCoupling * nGamma^(qQ - 1))/Binomial[nGamma , qQ])]] * Dot @@ gamma[[#2]] & , SparseArray[{1 , 1} -> 0 , {length , length}] , xCoupling]  , {length , length} , 0 ]
];
Return @ hamiltonianSYK
]

SYKHamiltonianBoson[qQ_ , nGamma_ , gamma_] := Block[{indices , iterators , hamiltonianSYKBoson} , 
indices = Table[Symbol["i" <> ToString[i]] , {i , 1 , qQ}];
iterators = Transpose @ {indices , Range[qQ , 1 , -1] , Flatten[{nGamma , Take[indices ,  qQ -  1] -1}]}; 
hamiltonianSYKBoson = (1/2)^(qQ / 2) * Sum[I^Total[(parityFunc @ Last @ ## * KroneckerDelta[First @ ## + 1 , Last @ ##] &) /@ Evaluate @ Partition[Reverse @ indices , 2 , 1]] *RandomVariate @ NormalDistribution[0. , Sqrt[Factorial[qQ - 1]/(nGamma^(qQ - 1))]] * (Dot @@ gamma[[##]] &) @ Evaluate @ Reverse @ indices  ,  ##] &  @@ iterators;
Return @ hamiltonianSYKBoson 
]

SpinChainHamiltonian[listCouplings_ , listCoefficients_] := Block[{latticeSize , length , sS , identityOperator , hamiltonianSpinChainUnitaries , hamiltonianSpinChain , nullList , hamiltonianSpinChainUnitariesNaive} , 
latticeSize = Max @ DeleteDuplicates @ Flatten @ listCouplings;
length = 2^latticeSize;
identityOperator = SparseArray[( {# , #} -> 1 &) /@ Range @ length , {length , length} , 0];
nullList = Position[listCoefficients , 0];
sS = Partition[(toMat @ PadRight[PadLeft[{#2} , #1 , id] ,latticeSize , id]  &) @@@ Tuples[{Range @ latticeSize , {s1 , s2 , s3}}] , 3];
hamiltonianSpinChainUnitariesNaive = Apply[Dot , Function[{x} , Map[sS[[# , x]] & , listCouplings , {2}]] /@ Range @ 3 // Transpose , {2}] * listCoefficients;
hamiltonianSpinChainUnitaries = Delete[hamiltonianSpinChainUnitariesNaive , nullList];
hamiltonianSpinChain = SparseArray[Total @ Flatten[hamiltonianSpinChainUnitaries , 1] , {length , length} , 0];
Return @ {hamiltonianSpinChainUnitaries , hamiltonianSpinChain}
]

QITE[latticeSize_ , hamiltonianDecomposed_ , listCouplings_ , \[Beta]_ , \[CapitalDelta]t_ , initialState_ , dD_ , toleranceNumeric_ :10^-6 , deltaDiagonal_ : 0.1] := Block[{sS , identityOperator , listHamiltonianIndices , evolvedState  , length = 2^latticeSize},

sS = Partition[(toMat @ PadRight[PadLeft[{#2} , #1 , id] ,latticeSize , id]  &) @@@ Tuples[{Range @ latticeSize , {s1 , s2 , s3}}] , 3];
identityOperator =SparseArray[( {# , #} -> 1 &) /@ Range @ length , {length , length} , 0];

listHamiltonianIndices = Catenate[(Reverse @ Flatten[Function[{x , y} ,({x , #} &) /@  Range @ y ]@@@ Transpose @ {Range @ Length @ hamiltonianDecomposed , Length /@ hamiltonianDecomposed} , 1] &) /@ Range[\[Beta] /\[CapitalDelta]t]];

evolvedState = FoldList[Function[{x , y} , stepEvolutionQITE[x , hamiltonianDecomposed[[First @ y , Last @ y]] , DeleteCases[(Mod[Range[#- dD , # + dD] , latticeSize , 1] &) /@ listCouplings[[First @ y]]//Flatten //DeleteDuplicates // Subsets , {}], sS , identityOperator , toleranceNumeric , deltaDiagonal , \[CapitalDelta]t]] , Normal @ initialState , listHamiltonianIndices];

Return[Transpose @ {Range[0 , \[Beta] , \[CapitalDelta]t] , Extract[Normal @ evolvedState , Partition[Range[1 , Length @ evolvedState , Length @ Flatten[hamiltonianDecomposed , 1]] , 1]]}]
]

End[]  (* QuantumManyBody`Private`*)

Protect[GammaMatrices , SpinOperators , EnergyStored , EvolvedStateList , FindGroundState , FindBandwidth , PartialTrace , EntanglementEntropy , XCouplingRegular , SYKHamiltonian , SYKHamiltonianBoson , SpinChainHamiltonian , QITE] (* Protect names of new functions*)

EndPackage[] (* QuantumManyBody` *)

