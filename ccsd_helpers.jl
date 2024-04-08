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

function cc_ket(H, T, n, r)
    # Construct bch(H, T, n) | HF > and keep only terms of excitation rank r or lower
    HT = bch(H, T, n) |> act_on_ket |> SpinAdaptedSecondQuantization.simplify
    terms = [length(t.operators) <= r for t in HT.terms]
    HT_rank = SASQ.Expression(HT[terms])
    return HT_rank
end


μ_1(a, i) = (1 // 2 * E(i, a))

μ_2(a, i, b, j) =
    (1 // 3 * E(i, a) * E(j, b) + 1 // 6 * E(i, b) * E(j, a))

T2 = ∑(1 // 2 * virtual(1, 3) * occupied(2, 4) *
       psym_tensor("t", 1, 2, 3, 4) * E(1, 2) * E(3, 4), 1:4)