using smearFEM
using LaTeXStrings
using Plots
using Plots.PlotMeasures

global fs = 40  # font size for plots

function compare_init(exp_params::Dict{String, Any})

    ndim::Int = 3
    nDof_p::Int = 1  # number of degree of freedom per node
    nDof_u::Int = ndim  # number of degree of freedom per node

    
    # simulation parameters for the ground truth
    filepath_gt::String = exp_params["filepath_gt"]
    # filepath_res::String = exp_params["filepath_res"]
    
    obj_pose::AbstractArray = zeros(Float64, 4,4) # initial pose of the object
    camera_matrix::AbstractArray = zeros(Float64, 4,4)
    
    # # sim_time_exp::Float64 = exp_params["sim_time_exp"]
    # if haskey(exp_params, "steps_exp")
    #     steps_exp = exp_params["steps_exp"]
    # else
    #     steps_exp = sim_time_exp*10.0 # assuming 30 fps for the experiments
    # end
    # t_steps_exp::Float64 = sim_time_exp/steps_exp
    
    # ne_exp::Int = exp_params["ne_exp"] # number of elements in the mesh for the experiment
    
    # SIDES::Bool = false
    
    # outlier_frames = Int[]
    
    params = read_json(joinpath(filepath_gt,"data","sim_params.json"))        
    r = params["r"]
    h = params["h"]
    
    sim_time_gt = params["simulation_time"]
    t_steps_gt = params["time_steps"]
    
    FunctionClass_x::String = "Q2"
    camera_matrix = reshape(Array(params["camera_matrix"]), 3, 3)
    obj_pose = reshape(Array(params["obj_pose"])*1.0,4,4)

    ObsDataList_sim, splinexObs_sim, splineyObs_sim = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
    ObsDataList_syn, splinexObs_syn, splineyObs_syn = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  

    obsBorderPts_sim, nSplinex, nSpliney, pd = add_noise(ObsDataList_sim, nFactor=0.0)
    symBorderPts_syn, splinex, spliney, pd = add_noise(ObsDataList_syn, nFactor=0.0)

    if length(obsBorderPts_sim) != length(symBorderPts_syn)
        n_time = min(length(obsBorderPts_sim), length(symBorderPts_syn))
        @warn "Mismatched number of time steps between observed border points ($(length(obsBorderPts_sim))) and simulated border points ($(length(symBorderPts_syn))). Truncating to $n_time time steps."
        obsBorderPts_sim = obsBorderPts_sim[1:n_time, :]
        symBorderPts_syn = symBorderPts_syn[1:n_time, :]
        nSplinex = nSplinex[1:n_time, :]
        nSpliney = nSpliney[1:n_time, :]

        splinex = splinex[1:n_time, :]
        spliney = spliney[1:n_time, :]
    end

    d, pairs = closest_point(symBorderPts_syn, obsBorderPts_sim)
    
    path = dirname(@__FILE__)
    exp_path = joinpath(path, "plots")
    set_file(exp_path)
    plt_cnt_error = set_plot(fs, sz=(2000, 1800))
    Plots.plot!(plt_cnt_error, d, label="Closest point distance error", dpi=400, lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -30mm)
    Plots.xlabel!(plt_cnt_error, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(plt_cnt_error, L"\mathrm{Closest\;Point\;Distance\;(px)}")
    Plots.savefig(plt_cnt_error, joinpath(exp_path,"closest_point_distance_error.pdf"))
    
    # write_json(joinpath(exp_path,"Results","data","experiment_parameters"), exp_params)

    # initialize_mesh(r, h, ne_exp, FunctionClass_x, camera_matrix, obj_pose, filepath, SIDES)
    
end

function main()
    exp_params = Dict{String, Any}()
    exp_params["filepath_gt"] = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/force/bulk_viscosity/Q2_16/1/"
    compare_init(exp_params) 
end

main()