# Reversal Distance model

# length of permutation
param N > 0;

# permutation pi
param Pi{1 .. N};

# permutation sigma
param Sigma{1 .. N};

# UB
param UB;

# LB
param LB;

# var B_{ijk}
var B{1 .. N, 1 .. N, 0 .. N-1} binary;

# var R_{abk}
var R{1 .. N, 1 .. N, 1 .. N-1} binary;

# var RD_{k}
var RD{0 .. N-1} binary;

# Constraints

# Common constraints
# Initial and final permutation are correct
s.t. initial_correct{i in 1 .. N}: 
     B[i, Pi[i],0] = 1;	  
s.t. final_correct{i in 1 .. N}: 
     B[i, Sigma[i],N-1] = 1;

# any position of permutation has exactly one value associated to it
s.t. pos1val{i in 1 .. N, k in 0 .. N-1}: 
     sum{j in 1 .. N} B[i,j,k] = 1;

# every value is assigned to one position of each permutation
s.t. val1pos{j in 1 .. N, k in 0 .. N-1}: 
     sum{i in 1 .. N} B[i,j,k] = 1;

# Revesal constraint
# if kth permutation does not alter the permutation, 
# none of subsequent will do so
s.t. rd0:
     RD[0] = 1;
s.t. alter_permutation{k in 1 .. N-1}: 
     RD[k] <= RD[k-1];

# at most one reversal is done at each step
s.t. one_reversal_is_done{k in 1 .. N-1}:
     (sum{a in 1 .. N-1} (sum{b in a+1 .. N} R[a,b,k])) = RD[k];

# changes in permutation by reversal
# i < a or i > b
s.t. rev1{i in 1 .. N, j in 1 .. N, k in 1 .. N-1}: 
     (sum{a in i+1 .. N-1} sum{b in a+1 .. N} R[a,b,k]) + 
     (sum{a in 1 ..  N-1} sum{b in a+1 .. i-1} R[a,b,k]) +
     B[i,j,k-1] + (1 - RD[k]) - B[i,j,k] <= 1;

# a <= i <= b
s.t. rev2{a in 1 .. N-1, b in a+1 .. N, i in a .. b, j in 1 .. N, k in 1 .. N-1}:
     	R[a,b,k] + B[b+a-i,j,k-1] - B[i,j,k] <= 1;

# objective
minimize rev_dist: sum{k in 1 .. N-1} RD[k];

end;
