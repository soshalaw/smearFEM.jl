using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

using smearFEM
using StatsPlots
using Distributions
using Dates

function run_exp(exp_params::Dict)

    WRITE_GT = exp_params["WRITE_GT"] 
    r::Float64 = 1  # radius of the cylinder in mm
    h::Float64 = 1  # height of the cylinder in
    filepath::String = ""
    β_gt::Float64 = 1.0
    η_gt::Float64 = 1.0
    sim_time_gt::Float64 = 1.0# simulation time in seconds
    steps_gt::Float64 = 1.0 # number of time steps
    t_steps_gt::Float64 = 1.0
    F::Vector{Float64} = ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

    ndim::Int = 3
    FunctionClass_x::String = ""
    FunctionClass_u::String = exp_params["FunctionClass_u"]
    nDof_u::Int = ndim  # number of degree of freedom per node
    FunctionClass_p::String = exp_params["FunctionClass_p"]
    nDof_p::Int = 1  # number of degree of freedom per node

    noiseLevel::Float64 = 0
    SIDES::Bool = false
    scale = 50
    r = 0.5*scale  # radius of the cylinder in mm
    h = 1*scale  # height of the cylinder in mm

    FunctionClass_x = exp_params["FunctionClass_x_gt"]
    control = exp_params["control"]
    viscosity_type = "constant" # "constant" or "bulk_viscosity"
    filepath = exp_params["filepath"]

    sim_time_gt = 20.0 # simulation time in seconds
    steps_gt = 20.0 # number of time steps
    t_steps_gt = sim_time_gt/steps_gt

    time = collect(Float64, range(start=t_steps, stop=sim_time, step=t_steps))

    β_gt = exp_params["β_gt"]
    η_gt = exp_params["η_gt"]
    ne_gt::Int = exp_params["ne_gt"] # number of elements in the mesh for the ground truth

    F_ext::Float64 = 1500.0*β_gt # force applied to the cylinder in N
    F = -F_ext*ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.5 4]'   # camera position in mm

    # Write the ground truth
    printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
    model_gt, scene_gt = def_problem(r, h, ne_gt, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                    sim_time_gt, t_steps_gt)

    @info "Writing ground truth gt data to with $ne_gt elements to $filepath"
    write_sim_data(model_gt, scene_gt, camera_matrix, camera_pose, filepath)

    ObsDataList, splinexObs, splineyObs = read_csv(string(filepath,"/data/sim_data/contour_data"))  
    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
    #sim data
    model_iga, scene_iga = def_problem(r, h, ne_gt, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                sim_time_gt, t_steps_gt)
    conditions_iga = Conditions(camera_matrix=camera_matrix, camera_pose=camera_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)

    est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model_iga, scene_iga, conditions_iga)
    d, ∂d, ∂2d, pairs = closest_point(borderPts2DList, obsBorderPts, gradList)

    # η = ηtst*(1-dev)
    # β = βtst*(1-dev)
    # cSample = zeros(steps+1,nSamples)

    # for n = 1:nSamples
    #     nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    #     ObsData = [nScene, nSplinex, nSpliney]

    #     if n == 1
    #         plot(x->pdf(pd, x))
    #         savefig(string(filepathi,"/Results"))
    #     end

    #     # hcost, cpCost = compare_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, CameraMatrix, endTime, tSteps, Control, 
    #     #             ObsData, SIDES, PLOT_MATCHES, filepathi)

    #     μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
    #                                                                         nDof_p, β, CameraMatrix, endTime, tSteps, Control, SIDES=SIDES)

    #     # test the closest point function
    #     d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

    #     println("sum ", sum(d))
    #     cSample[:,n] = d 
    # end

    plot_path  = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/caovergence_analysis/stokes_convergence/plots")
    # if nSamples == 1
    # cost = cSample/nSamples
    Plots.plot(time, d, label="Cost") 
    # else
    #     errorline(cSample, errorstyle=:ribbon, label="Cost")
    # end
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(plot_path,"/cost_cp.png"))

    params = Dict("Paramter type" => "Lame", "η" => ηtst, "β" => β, "Control" => Control, "Noise Level" => noiseLevel, 
                    "Sample_no" => nSamples)

    # write_json(string(filepathi,"/Results/params"), params) 
end

function main()
    ne_gt::Int = 8 # number of elements in the mesh for the ground truth
    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    β_gt_list = [5, 10, 50, 100.0, 200.0, 500.0, 1000.0, 10000.0]
    η_gt_list = [40.0]
    FunctionClass_x_List = ["S2", "Q2"]
    refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    FunctionClass_x_gt_list = ["S2", "Q2"] # Function space for the ground truth

    exp_size = size(FunctionClass_x_List,1)*size(refine_list,1)
    η_mat = zeros(30, exp_size)
    β_mat = zeros(30, exp_size)

    WRITE_GT = true
    for FunctionClass_x_gt in FunctionClass_x_gt_list
        run_id = 1
        file_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",control,"/model_acc_test_",FunctionClass_x_gt,"_gt_S2_updated_2")
        @info "Running optimization with FunctionClass_x_gt = $FunctionClass_x_gt with $ne_gt elements for the ground truth ..."
        for β_gt in β_gt_list
            for η_gt in η_gt_list
                filepath = string(file_path,"/",run_id)
                for ref in refine_list
                    ne = ne_exp^ref
                    @info "Running optimization with ne = $ne"
                    for FunctionClass_x in FunctionClass_x_List
                        @info "Running optimization with FunctionClass_x = $FunctionClass_x"
                        exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "ne_exp" => ne, 
                                    "β_gt" => β_gt, "η_gt" => η_gt, "filepath" => filepath, "control" => control, "FunctionClass_x_gt" => FunctionClass_x_gt)
                        test_opt_const(exp_params)
                        WRITE_GT = false
                    end
                end
                WRITE_GT = true
                run_id += 1
                post_analysis(filepath)
            end
        end
    end
end

main()