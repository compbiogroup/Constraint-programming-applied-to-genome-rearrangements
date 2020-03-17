# Transposition Distance model

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
var B{1 .. N+1, 1 .. N, 0 .. N-1} binary;

# var T_{abck}
var T{1 .. N+1, 1 .. N+1, 1 .. N+1, 1 .. N-1} binary;

# var TD_{k}
var TD{0 .. N-1} binary;

# var R_{abk}
var R{1 .. N, 1 .. N, 1 .. N-1} binary;

# var RD_{k}
var RD{0 .. N-1} binary;

# var ZD_{k}
var ZD{0 .. N-1} binary;


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

# Transposition constraint
# at most one transposition is done at each step
s.t. one_transpositionl_is_done{k in 1 .. N-1}: 
     sum{a in 1 .. N-1} 
     sum{b in a+1 .. N} 
     sum{c in b+1 .. N+1} T[a,b,c,k] = TD[k];

# changes in permutation by transposition
# i < a or i >= c
s.t. trans1{i in 1 .. N, j in 1 .. N, k in 1 .. N-1}: 
     (sum{a in i+1 .. N-1} sum{b in a+1 .. N} sum{c in b+1 .. N+1} T[a,b,c,k])
     + (sum{a in 1 .. N-1} sum{b in a+1 .. N} sum{c in b+1 .. i} T[a,b,c,k])
     + B[i,j,k-1] + (1 - ZD[k]) - B[i,j,k] <= 1;

# a <= i < a + c - b
s.t. trans2{a in 1 .. N+1, b in 1 .. N+1, c in 1 .. N+1,
     i in a .. a+c-b-1, j in 1 .. N, k in 1 .. N-1}: 
     if a < b and b < c then
     	T[a,b,c,k] + B[b-a+i,j,k-1] - B[i,j,k] <= 1;

# a + c - b <= i < c
s.t. trans3{a in 1 .. N+1, b in 1 .. N+1, c in 1 .. N+1,
     i in a+c-b .. c-1, j in 1 .. N, k in 1 .. N-1}: 
     if a < b and b < c then
     	T[a,b,c,k] + B[b-c+i,j,k-1] - B[i,j,k] <= 1;

# Revesal constraint
# at most one reversal is done at each step
s.t. one_reversal_is_done{k in 1 .. N-1}:
     (sum{a in 1 .. N-1} (sum{b in a+1 .. N} R[a,b,k])) = RD[k];

# changes in permutation by reversal
# i < a or i > b
s.t. rev1{i in 1 .. N, j in 1 .. N, k in 1 .. N-1}: 
     (sum{a in i+1 .. N-1} sum{b in a+1 .. N} R[a,b,k]) + 
     (sum{a in 1 ..  N-1} sum{b in a+1 .. i-1} R[a,b,k]) +
     B[i,j,k-1] + (1 - ZD[k]) - B[i,j,k] <= 1;

# a <= i <= b
s.t. rev2{a in 1 .. N-1, b in a+1 .. N, 
     i in a .. b, j in 1 .. N, k in 1 .. N-1}:
     	R[a,b,k] + B[b+a-i,j,k-1] - B[i,j,k] <= 1;

# Trans and rev
s.t. zd_0:
     ZD[0] = 1;
s.t. alter_permutation{k in 1 .. N-1}: 
     ZD[k] <= ZD[k-1];

s.t. trans_or_rev_is_done{k in 1 .. N-1}: 
     TD[k] + RD[k] = ZD[k];

# objective
minimize trans_dist: sum{k in 1 .. N-1} ZD[k];

end;
