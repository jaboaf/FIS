X = [
104 0 0 0 0;
4.5 104.5 0 0 0;
4.5 4.5 104.5 0 0;
3.5 3.5 3.5 103.5 0;
4.0 4.0 4.0 4.0 104
]

p = [
100;
99;
98;
92;
90
]

b = (X'X)^-1 * X' * p
r = [ (1/b[t])^(1/t) - 1 for t in 1:5 ]
β = b .^ -1
𝒻 = [ (β[t]/β[s]) for s in 1:5,t in 1:5 ]
f = [ (β[t]/β[s])^(1/(t-s)) for s in 1:5,t in 1:5 ] .- 1


