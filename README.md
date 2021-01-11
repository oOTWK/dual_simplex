# dual_simplex
Dual simplex method for LP relaxation of Set Cover Problem.

## How to run
1. `cd dual_simplex`
1. `make`
1. `./build/bin/dual_simplex file_path`
1. `make clean`

## References
[1] Koberstein, A. (2005). The dual simplex method, techniques for a fast and stable implementation. *Unpublished doctoral thesis, Universität Paderborn, Paderborn, Germany*.  
[2] Suhl, U. H., & Suhl, L. M. (1990). Computing sparse LU factorizations for large-scale linear programming bases. *ORSA Journal on Computing*, 2(4), 325-335.  
[3] Suhl, L. M., & Suhl, U. H. (1993). A fast LU update for linear programming. *Annals of Operations Research*, 43(1), 33-47.  
- Test problem sets  
scpnre, scpnrf, scpnrg, scpnrh from [Beasley’s OR library](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html)

## Computational results
scpnre1
itr: 263,
LP soln: 21.379416,
CPU time 0.218s

scpnre2
itr: 299,
LP soln: 22.360045,
CPU time 0.211s

scpnre3
itr: 304,
LP soln: 20.486142,
CPU time 0.245s

scpnre4
itr: 291,
LP soln: 21.352715,
CPU time 0.206s

scpnre5
itr: 274,
LP soln: 21.321921,
CPU time 0.182s

------------------------
scpnrf1
itr: 193,
LP soln: 8.795264,
CPU time 0.203s

scpnrf2
itr: 197,
LP soln: 9.993615,
CPU time 0.209s

scpnrf3
itr: 188,
LP soln: 9.492377,
CPU time 0.198s

scpnrf4
itr: 195,
LP soln: 8.471190,
CPU time 0.207s

scpnrf5
itr: 219,
LP soln: 7.835527,
CPU time 0.237s

------------------------
scpnrg1
itr: 1088,
LP soln: 159.886241,
CPU time 1.090s

scpnrg2
itr: 1130,
LP soln: 142.073321,
CPU time 1.327s

scpnrg3
itr: 1056,
LP soln: 148.269135,
CPU time 1.070s

scpnrg4
itr: 1133,
LP soln: 148.946521,
CPU time 1.098s

scpnrg5
itr: 1083,
LP soln: 148.231466,
CPU time 1.082s

------------------------
scpnrh1
itr: 774,
LP soln: 48.124555,
CPU time 1.171s

scpnrh2
itr: 870,
LP soln: 48.637625,
CPU time 1.381s

scpnrh3
itr: 739,
LP soln: 45.197462,
CPU time 1.208s

scpnrh4
itr: 721,
LP soln: 44.042108,
CPU time 1.152s

scpnrh5
itr: 729,
LP soln: 42.370359,
CPU time 1.168s
