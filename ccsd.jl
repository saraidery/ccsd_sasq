using SpinAdaptedSecondQuantization
include("ccsd_helpers.jl")

H = get_H()

eq = cc_ket(H, T2, 2, 2)

Ω_1 =  hf_expectation_value(μ_1(1,2)*occupied(2)*virtual(1) * eq) |> simplify_heavy
Ω_2 =  hf_expectation_value(μ_2(1,2,3,4)*occupied(2,4)*virtual(1,3) * eq) |> simplify_heavy

trans = translate(VirtualOrbital => 1:2:10, OccupiedOrbital => 2:2:10)
Ω_1 = (Ω_1, trans)
Ω_2 = (Ω_2, trans)

@show Ω_1
@show Ω_2


