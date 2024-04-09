
function generate_eT_code(file_name, routine_name, eq, eq_name, symmetry, translation, wf)
    open(file_name*".jl", "w") do io
        println(io, "include("*"\"../TensorOperation-eT-code/src/omeinsum_impl.jl\""*")")
        println(io, "open(\""*file_name*".F90\", \"w\") do io")
        println(io, "s = "*print_eT_function_generator(routine_name, eq, eq_name, symmetry, translation, wf,
            Dict([
                "u_vovo" => "wf%u_aibj",
                "t_vovo" => "wf%t_aibj",
                "F_oo" => "wf%fock_ij",
                "F_vv" => "wf%fock_ab",
                "F_ov" => "wf%fock_ia",
                "F_vo" => "wf%fock_ai",
            ]), ["F_oo", "F_vv", "F_ov", "F_vo", "u_vovo", "t_vovo"]))
        println(io, "println(io, s)")
        println(io, "end")
    end
end
