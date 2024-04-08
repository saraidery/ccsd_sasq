using SpinAdaptedSecondQuantization


function Tn(n)
	if n == 1
		return τ(1,2) * virtual(1) * occupied(2)
	elseif n > 1
		X(m) = prod(E(2i-1,2i) for i = 2:m) * occupied(4:2:2m...) * virtual(3:2:2m...)
		return τ(1,2) * virtual(1) * occupied(2) * X(n)
	else
		throw("n cannot be negative")
	end
end

function Tdn(n)
	if n == 1
		return τ(1,2) * occupied(1) * virtual(2)
	elseif n > 1
		X(m) = prod(E(2i-1,2i) for i = 2:m) * virtual(4:2:2m...) * occupied(3:2:2m...)
		return τ(1,2) * occupied(1) * virtual(2) * X(n)
	else
		throw("n cannot be negative")
	end
end

χ(n) = prod(E(2i-1,2i) for i = 1:n) * occupied(2:2:2n...) * virtual(1:2:2n...)
Xn(n) = 1 // factorial(n) * ∑(psym_tensor("t", 1:2n...) * χ(n), 1:2n)

Rn(n) = ∑(real_tensor("c", 1:2n...) * Tn(n), 1:2n)


function A_transpose_transform(H, X, Td, R, n)
    # <order| [bch(H, X, n),T] | HF >

    Ω = Td * commutator(bch(H, X, n), R)
    ex = simplify_heavy(act_on_ket(Ω,0))

    return ex
end

# Define Hamiltonian in terms of F and g
Φ = 1//2 * ∑(psym_tensor("g", 1,2,3,4) * e(1,2,3,4), 1:4) +
      ∑(-2 * psym_tensor("g", 1,2,3,3) * E(1,2) * occupied(3), 1:3) +
      ∑( 1 * psym_tensor("g", 1,3,3,2) * E(1,2) * occupied(3), 1:3)

F = ∑(real_tensor("F", 1, 2) * E(1, 2), 1:2)

H = F + Φ

ρ_ai_1 = A_transpose_transform(H, Xn(2), 1//2 * Tdn(1), Rn(1), 2)
ρ_ai_2 = A_transpose_transform(H, Xn(2), 1//2 * Tdn(1), Rn(2), 2)

trans = translate(VirtualOrbital => 2:2:10, OccupiedOrbital => 1:2:10)
println("\n ρ_ai (singles):\n", (ρ_ai_1, trans))
println("\n ρ_ai (non-symmetrized doubles):\n", (ρ_ai_2, trans))


