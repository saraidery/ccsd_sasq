using SpinAdaptedSecondQuantization
include("ccsd_helpers.jl")
include("code_generation_helpers.jl")

H = get_H()

# Write equations for Ω1 and Ω2
eq = cc_ket(H, T2, 2, 2)

Ω1 =  hf_expectation_value(μ_1(1,2)*occupied(2)*virtual(1) * eq) |> simplify_heavy
Ω2 =  hf_expectation_value(μ_2(1,2,3,4)*occupied(2,4)*virtual(1,3) * eq) |> simplify_heavy

# Replace g -> L and t -> u
Ω1 = look_for_tensor_replacements_smart(Ω1, make_exchange_transformer("g", "L"))
Ω1 = look_for_tensor_replacements_smart(Ω1, make_exchange_transformer("t", "u")) |> simplify_heavy

Ω2 = look_for_tensor_replacements_smart(Ω2, make_exchange_transformer("g", "L"))
Ω2 = look_for_tensor_replacements_smart(Ω2, make_exchange_transformer("t", "u")) |> simplify_heavy

# Find symmetries
Ω2_s, Ω2_ss, Ω2_ns = desymmetrize(Ω2, make_permutation_mappings([(1,2),(3, 4)]))

# Nice print
trans = translate(VirtualOrbital => 1:2:10, OccupiedOrbital => 2:2:10)
@show (Ω2_s, trans)
@show (Ω2_ss, trans)
@show (Ω2_ns, trans)


# Generate code
generate_eT_code("test", "omega_s", Ω1, "omega_ai", [1, 2], trans, "ccsd")

