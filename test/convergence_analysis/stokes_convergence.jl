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

        BorderPts2D_gt, SurfacePts2D_fem = extract_borders(new_node_list', camera_matrix, camera_pose)
        
        # display(Plots.plot!(BorderPts2D_gt[1,:],BorderPts2D_gt[2,:],label="gt"))
        pi_gt, qi_gt = fit_curve(border=BorderPts2D_gt)
        # set_plot(22)
        display(Plots.plot!(pi_gt,qi_gt,label="gt"))
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
            d_iga, pairs_iga = closest_point([vcat(pi_iga', qi_iga')], [vcat(pi_gt', qi_gt')])
            d_fem, pairs_fem = closest_point([vcat(pi_fem', qi_fem')], [vcat(pi_gt', qi_gt')])
            
            pts_d, pts_pairs = closest_point([BorderPts2D_iga], [vcat(pi_fem', qi_fem')])
            pts_d_iga, pts_pairs_iga = closest_point([BorderPts2D_iga], [vcat(pi_gt', qi_gt')])
            pts_d_fem, pts_pairs_fem = closest_point([BorderPts2D_fem], [vcat(pi_gt', qi_gt')])

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

main()