# t::Vetcor{Rational{UInt64}}
# C::Vector{UInt64}
# P::Rational{Int64}
local epsilon = 0.01

function YTM(P,C,t)
	push!(C,-P)
	push!(t,0)
	minimize( Polynomial(C,t), epsilon )
end
