# Fixed Income Securities
# PSET 1, 2/20/23
# Jojo Aboaf (jww262)
# Directory: https://github.com/jaboaf/FIS
# I like writing in programming languages with comments, it'll be a pretty pdf by some point this semester, (when I write a parser for my .jl style)

using Plots # grateful

## Here are 5 bonds with annual maturities, prices, and coupon information (the coupon rate is paid on a face value of 100).
P = [
99.39369844;
99.96463778;
100.5740901;
99.84682315;
100.1703903]

C = [
0;
0.01;
0.015;
0.015;
0.0175]


# # Q1) "Bootstrap" the 5 year yield curve r1 through r5.

# Cash flow matrix
X = Float64[ j<=i ? 100*C[i] + 100*(i==j) : 0 for i in 1:5, j in 1:5]
# Columns of X are linearly independent, therefore LS estimate of linear regression coefficients may be computed as below
β = (X'X)^-1 * X' * P
# B is the spot rate curve, naked i.e. 100 $_{τ+T} / P(τ,T) $_τ
B = β .^ -1
# r is the spot annualized return curve, naked rates de-compounded to 1 period rate less 1
r = [ (B[t])^(1/t) - 1 for t in 1:5 ]

# # Q2) Perform a Nelson Siegel interpolation of the same curve.
# I haven't written a program to do this yet, but heres my working method:
# 0: We have points in time. We want to find a NS function closest to the spot return curve
# Step: Picking a notion of 'close'; lets go with euclidean distance
# Want: Estimate r with a linear combination of functions of b; the vector of coefficients α=[α_1;α_2;α_3]   
# Want: Arguments α,b minimizing |M(b)α-r|_2
M(b) = Float64[ f(t) for t in 1:5, f in [x-> x, x-> (1-exp(-x/b))/x, x-> (1-(x+1)*exp(-x/b))/x] ]

# Step 1: M(b) has linearly independent columns so LS estimate of α (as a function of b) is
a(b) = (M(b)'M(b))^-1 * M(b)' * r
# Step 2: minimize 2-norm of M(b)â(b)-r
n(b) = sqrt(sum( (M(b)*a(b) - r) .^2 ) )
dom = 0.1:0.00001:10
plt1 = plot(dom,n,xlabel="x",ylabel="n(b)",label="n")
savefig(plt1, "n.png")
# Analytically: we can find a real zero of n'(b) by contour integration (take the edge of a thin strip around the positive real line)
# I made a janky computational approximation. (Newtons method is preferable but n is a super long formula divided by another super long formula and I didn't want to write that all out just to differentiate a rational function and I couldn't figure out how to tell Wolfram Alpha to do what I wanted it to do so here we are.)
# N0TE: Here 0.00001 is the increment, the precision is controlled by Float64 arithmetic 
img = diff(n.(dom)) ./ 0.00001
min = minimum(abs.(img))
argmin_I = findall( <=(min), abs.(img))
if length(argmin_I)==1
	b̂ = dom[argmin_I[1]]
	α = a(b̂)
	println("b̂ = $b̂ ")
	println("α = $(a(b̂))")
else 
	println("error: non-unique approximate zero of n' on domain ")
end

# # Q3) Using your answer from 2, what do you expect r10 and r20 to be?
NS(t) = ([1 (1-exp(-t/b̂))/t (1-(t+1)*exp(-t/b̂))/t] * a(b̂))[1]
plt2 = plot(0.01:0.01:20,NS,legend=false, title="Estimated Nelson-Seigel (b̂=$(round(b̂,digits=2)), α≈$(round.(a(b̂),digits=4))",xlabel="t")
scatter!(collect(1:5),r)
savefig(plt2, "NS.png")
println("I expect r(10) to be NS(10)=$(NS(10)) and r(20) to be NS(20)=$(NS(20))")
# BOOOOO JOJO this is not good. Why not? Figure it out man. Shape of the graph looks about right, but its shifted in the plane is in order.
# SOOOOO looks like I'll be expanding out this massive rational function to figure out whats up (increment too big? Float64 too imprecise ? Mistep in estimation procedure)
# Prof Note: I'll sort this out. I learn best by implementing things myself, and that way once I figure it out I'll have a program for an arbitrary cashflow matrix. 

# To check your answers, look for info on the January 2016 yield curve.

# # Q4) Suppose you are a holder of the 2 year zero.
#=
What would your one year holding period return (from t=0 to t=1) be under different assumptions for r1,2 ?
Let P(τ,T) denote the price of a bill at τ maturing in T time units.
Supposing P(t+1,1) := 100/(1+r_{1,2}),
Buying the 2 year zero for P(τ,2) $_τ at τ and selling it for P(t+1,1) := 100/(1+r_{1,2}) $_{τ+1} at τ+1 yields a return of 100/(1+r_{1,2}) $_{τ+1} - P(τ,2) $_τ. 
=#

# # Q5) Suppose you are a holder of the 5 year zero.
#=
What would your two year holding period return (from t=0 to t=2) be under different assumptions for r2,5?
Supposing P(t+2,3) := 100/(1+r_{2,5}),
Buying the 5 year zero for P(τ,5) $_τ at τ and selling it for P(t+2,3) := 100/(1+r_{2,5}) $_{τ+2} at τ+2 yields a return of 100/(1+r_{2,5}) $_{τ+2} - P(τ,5) $_τ. 
=#

# # Q6) Solve for the forward rate curve that gives one year forward rates from t=0 to t=5 (i.e f0,1, f1,2, f2,3, f3,4. f4,5).
# Forward rates implied are f_{t,T} = B(T)/B(t) where t<=T. (Forward prices implied are Face*f_{T,t} where t<=T)
f = [ B[T]/B[t] for t in 1:5, T in 1:5]

for T in reverse(2:5)
println([ "f_{$t,$T} =$(round(B[T]/B[t];digits=15))  " for t in 1:(T-1) ]...)
end

# Q7) Under expectations hypothesis, what is the expected 1 year spot rate from period 1 to 2 (r1,2)?
B[2]/B[1]

# Q8) Under expectations hypothesis, what is the expected three year spot rate from period 2 to 5 (r2,5)?
B[5]/B[2]

# Q9) Suppose these expectations are correct. Compare returns in questions 4 and 5 to r1 and r2.  Can you explain the relationship in words?
#=
The expectations hypothesis with a rate curve at one spot treats prices as a byproduct of forward rate surface, which we can invert to get prices, so returns (on a per face basis) become differences in rates.
If we use rates instead of returns, we get
(B[1]/B[2] ($_1/$_0)/($_2/$_0)) / (1/B[2] 1/($_2/$_0) )
= (B[1]/B[2] $_1/$_2 / (1/B[2] $_0/$_2 )
= B[1] $_1 / 1 $_0 
= B[1] $_1/$_0

Or using naked rates (or dressed units), we have
((B[1]/B[2]) $_1) / ((1/B[2]) $_0 )
= ((B[1]/B[2]) $_1) (B[2]/1) / ($_0)
= ((B[1]/1) $_1 )/ ($_0)
= (B[1]/1) $_1/$_0

Or using dressed units, we have
B[1] ($_1/$_0)

i) So, expectations hypothesis with a rate curve at one spot, yields no-arbitrage.
ii) Similarly, naked rates with dressed units delivers rates.
iii) Likewise, a spot-depenent forward surface.
iii) Restricting t and T to terms (i.e. amounts of time, changes in time), P(τ+t,T-t) F_{τ,t,T} = some ratio of distinct dollars
ii) When t and T any oriented terms (i.e. differences between in time) P(τ+t,T-t) F_{τ,t,T} = P((τ+t),T-t) F_{τ,t,T} = 1 USD at τ or any maturity spot = some ratio of distinct non-τ and non-maturity-spot dollars
i) When t and T vary over the same set, p(t,T) f_{t,T} = 1

Taking τ as naked time and assuming time is ordered degenerates to i).

Breifly, expectation hypothesis and exactly one spot curve forces no-arbitrage.
=#

# Q10) Finally, calculate all the expectations of 1 year spot rates from the current period through E(r4,5).  
[ f[t,t+1] for t in 1:4 ]

# Q11) Taking the expectations above as given, now recalculate the yield curve after the inclusion of a liquidity risk premium (k) of .5%.
# 50bp premiums at each forward rate gives
[ B[T]*(1.005)^(T-1) for T in 1:5]
# and 50bp premiums only at forwards with a term of 1 year gives
B .* 1.005