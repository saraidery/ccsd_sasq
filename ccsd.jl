using SpinAdaptedSecondQuantization

function get_H()

        # Define Hamiltonian in terms of F and g
        Φ = 1//2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
        ∑(-2 * psym_tensor("g", 1,2,3,3) * E(1,2) * occupied(3), 1:3) +
        ∑( 1 * psym_tensor("g", 1,3,3,2) * E(1,2) * occupied(3), 1:3)

        F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)

        H = F + Φ

        return H

end

H = get_H()

deex_1(a, i) =
    (1 // 2 * E(i, a))

deex_2(a, i, b, j) =
    (1 // 3 * E(i, a) * E(j, b) + 1 // 6 * E(i, b) * E(j, a))

T2 = ∑(1 // 2 * virtual(1, 3) * occupied(2, 4) *
       psym_tensor("t", 1, 2, 3, 4) * E(1, 2) * E(3, 4), 1:4)

eq = bch(H, T2, 2) |> simplify_heavy

Ω_1 =  hf_expectation_value(deex_1(1,2)*occupied(2)*virtual(1) * eq) |> simplify_heavy
Ω_2 =  hf_expectation_value(deex_2(1,2,3,4)*occupied(2,4)*virtual(1,3) * eq) |> simplify_heavy

trans = translate(VirtualOrbital => 1:2:10, OccupiedOrbital => 2:2:10)
@show (Ω_1, trans)
@show (Ω_2, trans)


