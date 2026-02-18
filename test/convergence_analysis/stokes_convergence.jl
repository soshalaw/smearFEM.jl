# test 3D FEM solution for 2x2x2 elements
using smearFEM

using Plots
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames

function main()
    # test case 
    scale = 50
    r = 0.5*scale  # radius of the cylinder in mm
    h = 1*scale  # height of the cylinder in mm
    ne = 4
    ndim = 3
    FunctionClass_x = "Q2"
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1

    # β_list = [1e-5 10 100 500 1000 1e5]
    β_list = [1e1, 1e2, 1e3, 1e5]
    exp_ne = [2, 4, 8]  # number of elements in radial and height directions,

    ν = 40.0
    μu_tp = -1.0
    μu_btm = 0
    μu_side = 0

    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm

    # error_iga_list = zeros(length(exp_ne),length(β_list))
    # error_fem_list = zeros(length(exp_ne),length(β_list))
    error_dif_list = zeros(length(exp_ne),length(β_list))

    contour_error_list = zeros(length(exp_ne),length(β_list))
    contour_error_list_iga = zeros(length(exp_ne),length(β_list))
    contour_error_list_fem = zeros(length(exp_ne),length(β_list))

    contour_pts_error_list = zeros(length(exp_ne),length(β_list))
    contour_pts_error_list_iga = zeros(length(exp_ne),length(β_list))
    contour_pts_error_list_fem = zeros(length(exp_ne),length(β_list))

    # q_ref, model_ref = simulate_single_tstep_stokes(r, h, 16, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                            # μu_side, FunctionClass_x="Q2")
    set_plot(22)
    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/single_step"
    counter::Int = 1
    for (j,β) in enumerate(β_list)  
        # read nutils data
        gt_filepath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen/slip_$counter/"
        @info "reading nutils files from $gt_filepath"
        
        sol_u_df = CSV.read(string(gt_filepath,"sol_u.csv"), DataFrame; header=false)
        sol_u = Matrix{Float64}(sol_u_df)

        node_list_df = CSV.read(string(gt_filepath,"node_list.csv"), DataFrame; header=false)
        node_list = Matrix{Float64}(node_list_df)

        new_node_list = node_list + sol_u

        BorderPts2D, SurfacePts2D_fem = extract_borders(new_node_list', camera_matrix, camera_pose)
        
        # display(Plots.plot!(BorderPts2D[1,:],BorderPts2D[2,:],label="gt"))
        pi, qi = fit_curve(border=BorderPts2D)
        # set_plot(22)
        display(Plots.plot!(pi,qi,label="gt"))
        for (i, ne) in enumerate(exp_ne)     
            
            println("Running for ne = $ne and β = $β")

            # check if the nodes are the same
            q_fem, model_fem = simulate_single_tstep_stokes(r, h, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                                    μu_side, FunctionClass_x="Q2")

            new_nodes_FEM = model_fem.mesh_x.NodeList + q_fem

            BorderPts2D_fem, SurfacePts2D_fem = extract_borders(new_nodes_FEM, camera_matrix, camera_pose, BorderNodesList = model_fem.mesh_u.side_nodes)
            pi_fem, qi_fem = fit_curve(border=BorderPts2D_fem)

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

            new_nodes_iga = (model_iga.mesh_x.NodeList + q_iga*T_)*T

            BorderPts2D_iga, SurfacePts2D_iga = extract_borders(new_nodes_iga, camera_matrix, camera_pose, BorderNodesList = model_iga.mesh_u.side_nodes)
            pi_iga, qi_iga = fit_curve(border=BorderPts2D_iga)

            @assert size(model_fem.mesh_x.NodeList) == size(NodeList_iga_proj) "NodeList sizes do not match"
            iter = 1:size(NodeList_iga_proj, 2)

            # display(Plots.plot!(pi_iga,qi_iga,label="IGA, ne = $ne"))
            d, pairs = closest_point([vcat(pi_iga', qi_iga')], [vcat(pi_fem', qi_fem')])
            d_iga, pairs_iga = closest_point([vcat(pi_iga', qi_iga')], [vcat(pi', qi')])
            d_fem, pairs_fem = closest_point([vcat(pi_fem', qi_fem')], [vcat(pi', qi')])
            
            pts_d, pts_pairs = closest_point([BorderPts2D_iga], [vcat(pi_fem', qi_fem')])
            pts_d_iga, pts_pairs_iga = closest_point([BorderPts2D_iga], [vcat(pi', qi')])
            pts_d_fem, pts_pairs_fem = closest_point([BorderPts2D_fem], [vcat(pi', qi')])

            norm_cost = pairs[2]
            norm_cost_iga = pairs_iga[2]
            norm_cost_fem = pairs_fem[2]

            norm_pts_cost = pts_pairs[2]
            norm_pts_cost_iga = pts_pairs_iga[2]
            norm_pts_cost_fem = pts_pairs_fem[2]
            
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
                # error_iga = error_iga + abs(sqrt(q_iga_proj[1,i]^2 + q_iga_proj[2,i]^2) - sqrt(q_ref[1,i]^2 + q_ref[2,i]^2))
                # error_fem = error_fem + abs(sqrt(q_fem[1,i]^2 + q_fem[2,i]^2) - sqrt(q_ref[1,i]^2 + q_ref[2,i]^2))
                # end

                # compare with the solution from the compatible Lagrange mesh
                error_dif = error_dif + abs(sqrt(q_iga_proj[1,i]^2 + q_iga_proj[2,i]^2) - sqrt(q_fem[1,i]^2 + q_fem[2,i]^2))
            end

            # error_iga_list[i] = error_iga/length(iter)
            # error_fem_list[i] = error_fem/length(iter)\
            error_dif_list[i,j] = error_dif/length(iter)

            contour_error_list[i,j] = sum(norm_cost)
            contour_error_list_iga[i,j] = sum(norm_cost_iga)
            contour_error_list_fem[i,j] = sum(norm_cost_fem)

            contour_pts_error_list[i,j] = sum(norm_pts_cost)
            contour_pts_error_list_iga[i,j] = sum(norm_pts_cost_iga)
            contour_pts_error_list_fem[i,j] = sum(norm_pts_cost_fem)

        end
        counter = counter + 1
    end
    set_file(filepath)
    write_csv(string(filepath,"/error_dif"),error_dif_list)
    write_csv(string(filepath,"/contour_pts_error"),contour_pts_error_list)
    write_csv(string(filepath,"/contour_pts_error_fem"),contour_pts_error_list_fem)
    write_csv(string(filepath,"/contour_pts_error_iga"),contour_pts_error_list_iga)
    write_csv(string(filepath,"/contour_error"),contour_error_list)
    write_csv(string(filepath,"/contour_error_fem"),contour_error_list_fem)
    write_csv(string(filepath,"/contour_error_iga"),contour_error_list_iga)

    plotpath = string(filepath,"/plots")
    set_file(plotpath) # create the directory to store the VTK files

    # set_plot(22)
    # Plots.plot!(exp_ne, error_iga_list, label="IGA", marker_shape=:circle, yscale=:log10, xlabel="Number of elements (ne)", ylabel="Average error")
    # Plots.plot!(exp_ne, error_fem_list, label="FEM", marker_shape=:diamond, yscale=:log10, xlabel="Number of elements (ne)", ylabel="Average error")
    # Plots.savefig(string(plotpath,"/stokes_convergence.pdf"))

    plt1 = set_plot(22)
    for (j,β) in enumerate(β_list)
        Plots.plot!(plt1, exp_ne, error_dif_list[:,j], label="β=$(β)", yscale=:log10, xlabel="Number of elements (ne)", ylabel="normalized error")
    end
    Plots.savefig(plt1, string(plotpath,"/stokes_pts_convergence_difference.pdf"))

    plt2 = set_plot(22)
    for (j,β) in enumerate(β_list)
        Plots.plot!(plt2, exp_ne, contour_pts_error_list[:,j], label="β=$(β)", yscale=:log10, xlabel="Number of elements (ne)", ylabel="Normalized error")
    end
    Plots.savefig(plt2, string(plotpath,"/stokes_pts_contour_convergence_difference.pdf"))

    plt3 = set_plot(22)
    for (j,β) in enumerate(β_list)
        Plots.plot!(plt3, exp_ne, contour_pts_error_list_iga[:,j], label="IGA, β=$(β)", yscale=:log10, legend_column=2)
        Plots.plot!(plt3, exp_ne, contour_pts_error_list_fem[:,j], label="FEM, β=$(β)", yscale=:log10, legend_column=2)
        Plots.xlabel!(plt3, "Number of elements (ne)")
        Plots.ylabel!(plt3, "Normalized error")
    end
    Plots.savefig(plt3, string(plotpath,"/stokes_pts_contour_convergence.pdf"))

    plt5 = set_plot(22)
    for (j,β) in enumerate(β_list)
        Plots.plot!(plt5, exp_ne, contour_error_list[:,j], label="β=$(β)", yscale=:log10, xlabel="Number of elements (ne)", ylabel="Normalized error")
    end
    Plots.savefig(plt5, string(plotpath,"/stokes_contour_convergence_difference.pdf"))

    plt4 = set_plot(22)
    for (j,β) in enumerate(β_list)
        Plots.plot!(plt4, exp_ne, contour_error_list_iga[:,j], label="IGA, β=$(β)", yscale=:log10, legend_column=2)
        Plots.plot!(plt4, exp_ne, contour_error_list_fem[:,j], label="FEM, β=$(β)", yscale=:log10, legend_column=2)
        Plots.xlabel!(plt4, "Number of elements (ne)")
        Plots.ylabel!(plt4, "Normalized error")
    end
    Plots.savefig(plt4, string(plotpath,"/stokes_contour_convergence.pdf"))

end

function mesh_convergence_analysis()

    height_list = AbstractArray[]
    height_error_list = AbstractArray[]
    border_error_list = Float16[]
    μ_list = AbstractArray[]
    time_list = []

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm

    β::Float64 = 100.0
    η::Float64 = 100.0

    viscosity_type::String = "constant" # "constant" or "bulk_viscosity"
    FunctionClass_x::String = "Q2" # Function space for the ground truth
    FunctionClass_p::String = "Q1"
    FunctionClass_u::String = "Q2"
    control::String = "force" # "force" or "velocity"
    nDof_u::Int = 3
    nDof_p::Int = 1
    ndim::Int = 3

    F_ext::Float64 = 0.2*9.812*1e3 # force applied to the cylinder in N
    sim_time::Float64 = 0.1 # simulation time in seconds
    step_size = 0.1
    steps = round(Int, sim_time/step_size)  

    obj_pose = zeros(Float64, 4, 4)
    obj_pose[1,1] = -1.0
    obj_pose[2,3] = -1.0
    obj_pose[3,2] = -1.0
    obj_pose[1:3,4] = [0.0, h/2, 150.0]
    camera_matrix::AbstractMatrix{Float64} = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$FunctionClass_x/convergence_analysis")

    conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=false, filepath=filepath, ANIMATE=false)
    F = -F_ext*ones(Float64, round(Int, steps)) # force applied to the cylinder in N

    model_ref, scene_ref = def_problem(r, h, 2, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, 
                                        sim_time, step_size)

    est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq, elapsed_time = stokes_single_step_force(model_ref, scene_ref, conditions)   

    mesh_sz_ = [2, 4, 6, 8, 10, 12, 14, 16]
    mesh_sz = reverse(mesh_sz_)
    iter_index = 1
    h_ref = 0.0
    border_ref = borderPts2DList
    for mesh in mesh_sz
        println("Running for mesh size = $mesh")
        # run the simulation for the given mesh size and store the results
        model, scene = def_problem(r, h, mesh, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, 
                                    sim_time, step_size)

        est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq, elapsed_time = stokes_single_step_force(model, scene, conditions)
        
        h_mesh = get_height(est_μ_list, h)
        if iter_index == 1
            h_ref = h_mesh
            border_ref = borderPts2DList
        end
        δh = abs.(h_ref - h_mesh)
        d, _ = closest_point(borderPts2DList, border_ref)

        println("Height error: ", δh)
        push!(height_list, h_mesh)
        push!(height_error_list, δh)
        push!(border_error_list, sum(d)/length(d))
        push!(μ_list, est_μ_list)
        push!(time_list, elapsed_time.value)
        iter_index += 1
    end

    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/mesh_convergence_analysis/height_list", height_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/mesh_convergence_analysis/height_error_list", height_error_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/mesh_convergence_analysis/border_error_list", border_error_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/mesh_convergence_analysis/μ_list", μ_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/mesh_convergence_analysis/mesh_sz", mesh_sz)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/convergence_analysis/stokes_convergence/mesh_convergence_analysis/time_list", time_list)
end
mesh_convergence_analysis()