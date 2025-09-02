# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test
using Plots
using LinearAlgebra

# test case 
scale = 100.0
r = 0.25*scale  # radius of the cylinder in mm
h = 0.5*scale  # height of the cylinder in mm
ne = 4
ndim = 3
FunctionClass_x = "Q2"
FunctionClass_u = "Q2"
FunctionClass_p = "Q1"
nDof_u = ndim  # number of degree of freedom per node
nDof_p = 1
β_list = [1e-5, 1e0, 1e2, 1e5, 1e10]
ν = 40.0
μu_tp = -3.0
μu_btm = 0
μu_side = 0

exp_ne = [2, 4, 8, 16]  # number of elements in radial and height directions,
# error_iga_list = zeros(length(exp_ne),length(β_list))
# error_fem_list = zeros(length(exp_ne),length(β_list))
error_dif_list = zeros(length(exp_ne),length(β_list))

# q_ref, model_ref = simulate_single_tstep_stokes(r, h, 16, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                        # μu_side, FunctionClass_x="Q2")

filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/caovergence_analysis/stokes_convergence"
for (i, ne) in enumerate(exp_ne)    
    for (j,β) in enumerate(β_list)    
        println("Running for ne = $ne and β = $β")
        #check if the nodes are the same
        q_fem, model_fem = simulate_single_tstep_stokes(r, h, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                                μu_side, FunctionClass_x="Q2")

        q_iga, model_iga = simulate_single_tstep_stokes(r, h, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                                μu_side, FunctionClass_x="S2")

        T = Matrix{Float64}(I, size(model_iga.mesh_u.NodeList,2), size(model_iga.mesh_u.NodeList,2))
        if model_iga.mesh_x.FunctionClass == "S2" && model_iga.mesh_u.FunctionClass != "S2"
            T = get_nurbs_2_lagrange_proj(model_iga.mesh_x.IEN, model_iga.mesh_u.IEN, model_iga.mesh_x.C_vol, model_iga.mesh_x.NodeList, model_iga.mesh_x.W)
        end
        T_ = T'*inv((T*T'))

        q_iga_proj = q_iga*T_*T
        NodeList_iga_proj = model_iga.mesh_x.NodeList*T
        IEN_iga_proj = model_iga.mesh_u.IEN

        @assert size(model_fem.mesh_x.NodeList) == size(NodeList_iga_proj) "NodeList sizes do not match"
        iter = 1:size(NodeList_iga_proj, 2)

        error_iga = 0.0
        error_fem = 0.0
        error_dif = 0.0
        for i in iter
            # if β == 1e-5
            #     r_iga = sqrt(NodeList_iga_proj[1,i]^2 + NodeList_iga_proj[2,i]^2)
            #     r_fem = sqrt(model_fem.mesh_x.NodeList[1,i]^2 + model_fem.mesh_x.NodeList[2,i]^2)

            #     error_iga = error_iga + abs(sqrt(q_iga_proj[1,i]^2 + q_iga_proj[2,i]^2) + 0.5*μu_tp*r_iga/h)
            #     error_fem = error_fem + abs(sqrt(q_fem[1,i]^2 + q_fem[2,i]^2) + 0.5*μu_tp*r_fem/h)
            # else
            #     error_iga = error_iga + abs(sqrt(q_iga_proj[1,i]^2 + q_iga_proj[2,i]^2) - sqrt(q_ref[1,i]^2 + q_ref[2,i]^2))
            #     error_fem = error_fem + abs(sqrt(q_fem[1,i]^2 + q_fem[2,i]^2) - sqrt(q_ref[1,i]^2 + q_ref[2,i]^2))
            # end

            # compare with the solution from the compatible Lagrange mesh
            error_dif = error_dif + abs(sqrt(q_iga_proj[1,i]^2 + q_iga_proj[2,i]^2) - sqrt(q_fem[1,i]^2 + q_fem[2,i]^2))
        end

        # error_iga_list[i] = error_iga/length(iter)
        # error_fem_list[i] = error_fem/length(iter)
        error_dif_list[i,j] = error_dif/length(iter)
    end
end

plotpath = string(filepath,"/plots")
# set_plot(22)
# set_file(plotpath) # create the directory to store the VTK files
# Plots.plot!(exp_ne, error_iga_list, label="IGA", marker_shape=:circle, yscale=:log10, xlabel="Number of elements (ne)", ylabel="Average error")
# Plots.plot!(exp_ne, error_fem_list, label="FEM", marker_shape=:diamond, yscale=:log10, xlabel="Number of elements (ne)", ylabel="Average error")
# Plots.savefig(string(plotpath,"/stokes_convergence.pdf"))

set_plot(22)
for (j,β) in enumerate(β_list)
    Plots.plot!(exp_ne, error_dif_list[:,j], label="β=$(β)", marker_shape=:circle, yscale=:log10, xlabel="Number of elements (ne)", ylabel="Average error")
end
Plots.savefig(string(plotpath,"/stokes_convergence_difference.pdf"))
