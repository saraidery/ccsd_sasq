using SpinAdaptedSecondQuantization
include("ccsd_helpers.jl")
include("../TensorOperation-eT-code/src/omeinsum_impl.jl")

H = get_H()

eq = cc_ket(H, T2, 2, 2)

Ω1 =  hf_expectation_value(μ_1(1,2)*occupied(2)*virtual(1) * eq) |> simplify_heavy
Ω2 =  hf_expectation_value(μ_2(1,2,3,4)*occupied(2,4)*virtual(1,3) * eq) |> simplify_heavy

Ω1 = look_for_tensor_replacements_smart(Ω1, make_exchange_transformer("g", "L"))
Ω1 = look_for_tensor_replacements_smart(Ω1, make_exchange_transformer("t", "u")) |> simplify_heavy

Ω2 = look_for_tensor_replacements_smart(Ω2, make_exchange_transformer("g", "L"))
Ω2 = look_for_tensor_replacements_smart(Ω2, make_exchange_transformer("t", "u")) |> simplify_heavy

Ω2_s, Ω2_ss, Ω2_ns = desymmetrize(Ω2, make_permutation_mappings([(1,2),(3, 4)]))

trans = translate(VirtualOrbital => 1:2:10, OccupiedOrbital => 2:2:10)

@show (Ω2_s, trans)
@show (Ω2_ss, trans)
@show (Ω2_ns, trans)

# trans = translate(VirtualOrbital => 1:2:10, OccupiedOrbital => 2:2:10)
# Ω1 = (Ω1, trans)
# Ω2 = (Ω2, trans)
# @show Ω1
# @show Ω2
