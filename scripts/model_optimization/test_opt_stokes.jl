using LinearAlgebra
using ProgressMeter
using SparseArrays
using NearestNeighbors

using smearFEM
using StatsPlots
using Distributions
using Dates
using Plots
using LaTeXStrings
using DelimitedFiles

using ArgCheck
using Base.Threads
using Random
using Colors

using Dates
using Plots.PlotMeasures

include("../ParallelExecution.jl")
using .ParallelExecution

global def_orange = RGB(245/255,118/255,0)
global def_blue = RGB(5/255,79/255,185/255)
global def_red = RGB(196/255,70/255,1/255)
global def_green = RGB(2/255,147/255,86/255)
global end_obs_win = 20.1

# Plot configuration, mirroring the PLOT_CONFIG pattern in
# test/convergence_analysis/forward_convergence.jl. PLOT_CONFIG is the
# fallback/default geometry (used to seed the globals below at load time);
# PLOT_PRESETS holds the full per-(data_type[, viscosity_type]) override
# (geometry + y-limits) applied by set_plot_config() — that's the single
# place to edit sizes for a specific data type. Restored from the last
# committed values (git show 65269dc), where "synthetic"/"bulk_viscosity"
# used a narrower 1/3-linewidth plot than the other cases.
# for 1/2 linewidth
const PLOT_CONFIG = Dict(
    :font_size => 12,
    :plot_height => 320,
    :plot_width => 480,
    :left_margin => 1pt,
    :right_margin => 5pt,
    :top_margin => 1pt,
)

# for 1/3 linewidth
# const PLOT_CONFIG = Dict(:font_size => 12, :plot_height => 360, :plot_width => 330,
#                           :left_margin => -6pt, :right_margin => 10pt, :top_margin => 0pt)

# for 1/4 linewidth
# const PLOT_CONFIG = Dict(:font_size => 10, :plot_height => 300, :plot_width => 239,
#                           :left_margin => -6pt, :right_margin => 0pt, :top_margin => -1pt)

global fs::Int = PLOT_CONFIG[:font_size]
global plt_height::Int = PLOT_CONFIG[:plot_height]
global plt_width::Int = PLOT_CONFIG[:plot_width]
global plt_lft_margin = PLOT_CONFIG[:left_margin]
global plt_right_margin = PLOT_CONFIG[:right_margin]
global plt_top_margin = PLOT_CONFIG[:top_margin]

# Per-(data_type[, viscosity_type]) geometry + y-limit presets, applied by
# set_plot_config(). Keyed by data_type alone except for "synthetic", which
# is further split by viscosity_type since its geometry differs.
const PLOT_PRESETS = Dict(
    "physical" => Dict(
        :font_size => 12, :plot_height => 350, :plot_width => 477,
        :left_margin => 1pt, :right_margin => 5pt, :top_margin => 1pt,
        :y_lims_h_norm => (0.8, 1.05), :y_lims_rel_error => (-0.05, 20),
    ),
    ("synthetic", "constant") => Dict(
        :font_size => 12, :plot_height => 350, :plot_width => 477,
        :left_margin => 1pt, :right_margin => 5pt, :top_margin => 1pt,
        :y_lims_h_norm => (0.97, 1.02), :y_lims_rel_error => (-0.1, 3.0),
    ),
    ("synthetic", "bulk_viscosity") => Dict(
        :font_size => 12, :plot_height => 360, :plot_width => 330,
        :left_margin => -6pt, :right_margin => 10pt, :top_margin => 0pt,
        :y_lims_h_norm => (0.97, 1.02), :y_lims_rel_error => (-0.1, 3.0),
    ),
    "simulated" => Dict(
        :font_size => 12, :plot_height => 350, :plot_width => 477,
        :left_margin => 1pt, :right_margin => 5pt, :top_margin => 1pt,
        :y_lims_h_norm => (0.995, 1.005), :y_lims_rel_error => (-0.05, 0.1),
    ),
)

# for synthetic data (default preset until set_plot_config() is called)
global y_lims_h_norm = PLOT_PRESETS[("synthetic", "constant")][:y_lims_h_norm]
global y_lims_rel_error = PLOT_PRESETS[("synthetic", "constant")][:y_lims_rel_error]

# Extended colorblind-friendly palette (20 colors, optimized for maximum distinction)
# Colors are spaced across hue spectrum with varying brightness/saturation
global color_palette = [
    RGB(230/255, 159/255, 0/255),       # 1. orange (warm, bright)
    RGB(0/255, 114/255, 178/255),       # 2. medium blue (cool)
    RGB(204/255, 0/255, 0/255),         # 3. red (warm)
    RGB(0/255, 158/255, 115/255),       # 4. green (cool, mid-tone)
    RGB(204/255, 41/255, 130/255),      # 5. magenta/pink (warm)
    RGB(0/255, 0/255, 153/255),         # 6. dark blue (cool, dark)
    RGB(213/255, 94/255, 0/255),        # 7. vermillion (warm, dark)
    RGB(75/255, 0/255, 130/255),        # 8. indigo (cool, dark)
    RGB(255/255, 127/255, 80/255),      # 9. coral (warm, light)
    RGB(86/255, 180/255, 233/255),      # 10. sky blue (cool, light)
    RGB(178/255, 140/255, 51/255),      # 11. tan/brown (neutral, warm)
    RGB(0/255, 150/255, 200/255),       # 12. cyan (cool, bright)
    RGB(153/255, 0/255, 76/255),        # 13. deep magenta (warm, dark)
    RGB(51/255, 153/255, 51/255),       # 14. forest green (cool, dark)
    RGB(200/255, 50/255, 100/255),      # 15. deep rose (warm)
    RGB(0/255, 120/255, 100/255),       # 16. dark teal (cool, dark)
    RGB(102/255, 0/255, 204/255),       # 17. purple (cool, bright)
    RGB(150/255, 75/255, 0/255),        # 18. dark brown (neutral, dark)
    RGB(220/255, 80/255, 150/255),      # 19. light pink (warm, light)
    RGB(89/255, 89/255, 89/255)         # 20. dark gray (neutral)
]

"""
    optimize(exp_params)

Run the full parameter-estimation pipeline for one experiment configuration:
build the ground-truth (or physical) model and observation data described by
`exp_params`, fit viscosity `η` (and slip `β`) to the observed contours via
Gauss-Newton (`fit_model`), then re-simulate with the estimated parameters and
write all results (fitted parameters, cost/iteration history, 2D/3D fields,
heights) to `exp_params["filepath_res"]`.

Supports `exp_params["data_type"] ∈ ("simulated"|"synthetic", "physical")`,
`viscosity_type ∈ ("constant", "bulk_viscosity")` (the latter optionally
windowed in time via `exp_params["mode"]`), and multiple camera view angles
(`exp_params["z_angle_list"]` in the ground-truth data).

# Arguments
- `exp_params::Dict`: experiment configuration (data/viscosity type, mesh,
  camera, and output-path parameters).

# Returns
None. All outputs are side effects written to disk.
"""
function optimize(exp_params::Dict)

    η_gt::Vector{Float64} = [0.0]
    β_gt::Vector{Float64} = [0.0]
    F_ext::Float64 = 0.0
    sim_time_gt::Float64 = 0.0
    t_steps_gt::Float64 = 0.0
    steps_exp::Float64 = 0.0
    outlier_frames::Vector{Int} = Int[]

    ndim::Int = 3
    nDof_p::Int = 1
    nDof_u::Int = ndim

    element_shape_u = Symbol(exp_params["element_shape_u"])
    basis_order_u   = Int(exp_params["basis_order_u"])
    element_shape_p = Symbol(exp_params["element_shape_p"])
    basis_order_p   = Int(exp_params["basis_order_p"])
    element_shape_x = Symbol(exp_params["element_shape_x"])
    basis_order_x   = Int(exp_params["basis_order_x"])
    
    _mesh_path_kw(d) = haskey(d, "mesh_path") ? (mesh_path=d["mesh_path"],) : (;)
    _fem = (element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x)

    iterList::Vector{Float64} = Vector{Float64}()
    costList::Vector{Float64}  = Vector{Float64}()
    ηpList::Vector{Float64} = Vector{Float64}()
    βpList::Vector{Float64} = Vector{Float64}()
    
    # simulation parameters for the ground truth
    start_time = Dates.now()
    filepath_gt::String = exp_params["filepath_gt"]
    filepath_res::String = exp_params["filepath_res"]
    
    obj_pose::Vector{Float64} = zeros(Float64, 3)
    camera_matrix::AbstractArray = exp_params["camera_matrix"]
    
    sim_time_exp::Float64 = exp_params["sim_time_exp"]
    t_steps_exp::Float64 = if haskey(exp_params, "dt")
        steps_exp = sim_time_exp / Float64(exp_params["dt"])
        Float64(exp_params["dt"])
    elseif haskey(exp_params, "steps_exp")
        steps_exp = Float64(exp_params["steps_exp"])
        sim_time_exp / steps_exp
    else
        steps_exp = sim_time_exp / 0.1
    end
    
    ne_exp::Int = exp_params["ne_exp"]
    
    r::Float64 = 0.0
    h::Float64 = 0.0
    data_type ::String = exp_params["data_type"] # "simulated" or "physical" or "real"
    geometry::Union{Symbol, Nothing} = nothing
    z_angles::Vector{Float64} = [0.0]
    edge_radius::Union{Float64,Nothing} = nothing
    noiseLevel::Float64 = 0.0
    viscosity_model::String = "" # "constant" or "bulk_viscosity"

    obs_data_angles = AbstractArray[]
    splinex_angles = AbstractArray[]
    spliney_angles = AbstractArray[]

    if data_type == "simulated" || data_type == "synthetic"
        
        WRITE_GT = exp_params["WRITE_GT"] 
        noiseLevel = exp_params["noise_level"]
        outlier_frames = Int[1]

        sim_params = read_json(joinpath(filepath_gt,"data","sim_params"))
        
        r = sim_params["r"]
        h = sim_params["h"]
        geometry = get!(sim_params, "geometry", :cylinder)
        z_angles = get!(sim_params, "z_angle_list", [0.0])
        edge_radius = get!(sim_params, "edge_radius", nothing)
        gt_viscosity_type = sim_params["viscosity_type"]
        F = sim_params["cParam"]

        sim_time_gt = sim_params["simulation_time"]
        t_steps_gt = sim_params["time_steps"]
        
        camera_matrix = sim_params["camera_matrix"]
        _obj_pose = sim_params["obj_pose"]
        obj_pose = if size(_obj_pose, 1) == 4 && size(_obj_pose, 2) == 4
                        _obj_pose[[3,1,2],4]
                    else
                        _obj_pose
                    end
        
        control = sim_params["control_type"]

        if WRITE_GT == true
            @info "Writing ground truth gt data to with $ne_exp elements to $exp_path"
            write_gt_data(exp_params)
        end
        
        if gt_viscosity_type == "bulk_viscosity"
            if haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
                viscosity_model = sim_params["model_type"]
            else
                viscosity_model = "power_law"
                η_gt = sim_params["η"]
                β_gt = sim_params["β"]
            end
        else
            η_gt = sim_params["η"]
            β_gt = sim_params["β"]
        end

        _shared = (η_gt[1], element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x, β_gt[1], F, control, gt_viscosity_type, sim_time_gt, t_steps_gt)
        model_gt, scene_exp = if geometry === :cylinder
                                def_problem(Cylinder(r, h), ne_exp, _shared...; _mesh_path_kw(exp_params)...)
                            elseif geometry === :cube
                                lx = Float64(get(exp_params, "lx", r))
                                ly = Float64(get(exp_params, "ly", r))
                                lz = Float64(get(exp_params, "lz", h))
                                def_problem(Cuboid(lx, ly, lz, edge_radius), ne_exp, _shared...; _mesh_path_kw(exp_params)...)
                            else
                                error("Unsupported geometry type: $geometry")
                            end
        
        for i::Int in 1:length(z_angles)
            if data_type == "synthetic"
                ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","view_$i","contour_data"))
            elseif data_type == "simulated"
                ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","view_$i","contour_data"))  
            end

            push!(obs_data_angles, ObsDataList)
            push!(splinex_angles, splinexObs)
            push!(spliney_angles, splineyObs)
        end

    elseif data_type == "physical"
        r = exp_params["r"]  # radius of the cylinder in mm
        h = exp_params["h"]  # height of the cylinder in mm
        
        control::String = exp_params["control"]
        gt_viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"

        t_obs, t_top_plt, t_btm_plt = read_perception_data(joinpath(filepath_gt, "data", "sequence.hdf5"))

        _ObsDataList, _splinexObs, _splineyObs = read_csv(joinpath(filepath_gt, "data", "img_data", "contour_data"))
        meta_data = read_json(joinpath(filepath_gt, "data", "video_metadata"))

        frame_rate::Float64 = round(meta_data["frame_rate"], digits=1)
        frame_width::Int = meta_data["frame_width"]
        frame_height::Int = meta_data["frame_height"]
        compression_frames = meta_data["compressed_frames"]
        
        sim_time_gt = length(compression_frames)/frame_rate  # seconds
        steps_exp = sim_time_exp*frame_rate
        t_steps_exp = 1/frame_rate
        t_steps_gt = t_steps_exp
        
        F_ext = exp_params["F_ext"]
        println("Applied force: $F_ext N")
        
        F = -F_ext*ones(Float64, round(Int, sim_time_gt*frame_rate)) # force applied to the cylinder in N
        
        obj_pose = get_pose(t_obs)
        obj_pose_ = zeros(Float64, 4,4)
        obj_pose_[1,1] = -1.0
        obj_pose_[2,3] = -1.0
        obj_pose_[3,2] = -1.0
        obj_pose_[1:3,4] = obj_pose[1:3,4]
        
        push!(ObsDataList, _ObsDataList[compression_frames])
        push!(splinex_angles, _splinexObs[compression_frames])
        push!(spliney_angles, _splineyObs[compression_frames])

        function get_gt_height(t_top_plt::AbstractArray, t_btm_plt::AbstractArray, h_gt_init::Float64)
            n_frames = size(t_top_plt,2)
            Δh_gt = zeros(Float64, n_frames)
            h_init = abs(t_top_plt[2,1] - t_btm_plt[2,1])
            for i::Int in 1:n_frames
                Δh_gt[i] = h_init - abs(t_top_plt[2,i] - t_btm_plt[2,i])
            end

            h_gt = h_gt_init .- Δh_gt

            return h_gt, Δh_gt
        end

        h_gt , _ = get_gt_height(t_top_plt[1:3,4,compression_frames], t_btm_plt[1:3,4,compression_frames], 38.5)

        sim_data = Dict("camera_matrix" => camera_matrix,
                    "obj_pose" => obj_pose,
                    "time_steps" => t_steps_exp,
                    "simulation_time" => sim_time_gt,
                    "cParam" => F,
                    "r" => r,
                    "h" => h,
                    "control_type" => control,
                    "viscosity_type" => gt_viscosity_type,
                    "data_type" => data_type)

        write_json(joinpath(filepath_gt,"data","sim_params"), sim_data)
        write_csv(joinpath(filepath_gt,"data","h"), h_gt)

        valid_frames, _outlier_frames = detect_outlier_observations(_ObsDataList[compression_frames])
        outlier_frames = vcat(outlier_frames, _outlier_frames)
        geometry = :cylinder
    else
        error("data_type should be either simulated or physical")
    end

    geom_exp = if geometry === :cylinder
        Cylinder(r, h)
    elseif geometry === :cube
        lx = Float64(get(exp_params, "lx", r))
        ly = Float64(get(exp_params, "ly", r))
        lz = Float64(get(exp_params, "lz", h))
        Cuboid(lx, ly, lz, edge_radius)
    else
        error("Unsupported geometry type: $geometry")
    end

    if sim_time_gt < sim_time_exp
        @warn "Ground truth simulation time $sim_time_gt is less than experimental simulation time $sim_time_exp , switching to ground truth simulation time"
        sim_time_exp = sim_time_gt
    end

    if t_steps_gt > t_steps_exp
        @warn "time resolution of the ground truth $t_steps_gt is larger than the experimental $t_steps_exp, switching to ground truth resolution"
        t_steps_exp = t_steps_gt
    end
    
    for (i,z_angle) in enumerate(z_angles)
        
        printstyled("Processing view $i with z_angle = $z_angle degrees\n"; color = :blue)
        if data_type == "physical"  || viscosity_model == "carreau"
            η_start = exp_params["η_start"]
            β_start = exp_params["β_start"]
        else
            dev::Float64 = 0.3
            
            dev_η::Float64 = dev*η_gt[1]
            η_start::Float64 = abs(η_gt[1] - dev_η)
            
            if β_gt[1] >= 200.0 # setting the slip to partial slip for no slip cases
                β_start = 50.0
            elseif β_gt[1] <= 1 # setting the slip to partial slip for free slip cases
                β_start = 10.0
            else
                dev_β::Float64 = dev*β_gt[1]
                β_start::Float64 = abs(β_gt[1] - dev_β)
            end
        end
        
        θ::Vector{Float64} = [η_start, β_start]

        ObsDataList = obs_data_angles[i]
        splinexObs = splinex_angles[i]
        splineyObs = spliney_angles[i]
        
        exp_params["z_angle"] = z_angle
        exp_path = joinpath(dirname(filepath_res), "view_$i", basename(filepath_res))
        write_json(joinpath(exp_path, "data","experiment_parameters"), exp_params)

        if gt_viscosity_type == "constant"
            
            printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
            _range = collect(range(start = 1, stop = (round(Int,sim_time_exp/t_steps_gt)+1), step = round(Int,t_steps_exp/t_steps_gt)))
            println(round(Int,sim_time_exp/t_steps_gt)+1," ",round(Int,t_steps_exp/t_steps_gt)," ",length(_range))
            @info "Considering from frame $(first(_range)) to frame $(last(_range)) in the observations"
            ObsDataList = ObsDataList[_range] # align the observation points with the simulation time
            printstyled("Observation data length: $(length(ObsDataList))\n"; color = :blue)
            model, scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F, control, gt_viscosity_type, sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)
            est_model, est_scene = def_problem(geom_exp, ne_exp, η_start, _fem..., β_start, F, control, gt_viscosity_type, sim_time_gt, t_steps_exp; _mesh_path_kw(exp_params)...)
            conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, filepath=exp_path, ANIMATE=false, viewing_angles=[z_angle])
            
            gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
            gt_h = gt_h_[1:(round(Int,sim_time_gt/t_steps_exp)+1)]

            if noiseLevel == 0.0
                    
                obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
                stats = fit_model(model, scene, conditions, obs_border_pt_lst, θ, outliers=outlier_frames, method=:gn)

                iterList = stats["iterList"]
                costList = stats["cost_list"]
                ηpList = stats["ηList"]
                βpList = stats["βList"]

                η = stats["η"]
                β = stats["β"]

                printstyled("Estimated η : $(η), estimated β: $(β)\n"; color = :green)

                η_accuracy = (1-abs((η_gt[1]-η)/η_gt[1]))*100
                β_accuracy = (1-abs((β_gt[1]-β)/β_gt[1]))*100
                printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

                set_file(joinpath(exp_path,"plots"))
                
                write_json(joinpath(exp_path,"data","view_$i","stats"), stats)
                write_csv(joinpath(exp_path,"data","view_$i","η"), ηpList)
                write_csv(joinpath(exp_path,"data","view_$i","β"), βpList)
                write_csv(joinpath(exp_path,"data","view_$i","gt_h"), gt_h)
                write_csv(joinpath(exp_path,"data","view_$i","cost_iter"), costList)
                
                reset_model!(est_model)
                est_model.η = [η]
                est_scene.β = [β]
                
                conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, filepath=exp_path, ANIMATE=true, viewing_angles=[z_angle])
                est_μ_list, gradList, borderPts2DList, fields, surface_pts_3D_gt, pos2D, pos3D, _, _, _ = simulate(est_model, est_scene, conditions)
                est_h = get_height(est_μ_list, h)

                write_csv(joinpath(exp_path,"data","est_h"), est_h)
                write_data(joinpath(exp_path,"data","sim_data","3D_point_est"), pos3D)
                write_data(joinpath(exp_path,"data","sim_data","motion_fields_est "), fields)
                write_data(joinpath(exp_path,"data","sim_data","3D_surface_points_est"), surface_pts_3D_gt)

                write_2d_data(joinpath(exp_path,"data","sim_data","2D_surface_points_est"), pos2D)
                write_2d_data(joinpath(exp_path,"data","sim_data","2D_border_points_est"), borderPts2DList)

                if maximum(ηpList) > η+dev_η
                    ηStop = maximum(ηpList)*1.1
                else
                    ηStop = η+dev_η
                end

                if minimum(ηpList) < η-dev_η
                    if minimum(ηpList) <= 0.0
                        η_start = 0.1
                    else
                        η_start = minimum(ηpList)*0.9
                    end
                else
                    η_start = η-dev_η
                end

            else
                n_samples = 10
                η_pred = zeros(Float64, n_samples)
                β_pred = zeros(Float64, n_samples)
                costnList = Vector{Vector{Float64}}(undef, n_samples)
                iternList = Vector{AbstractVector}(undef, n_samples)
                est_h_list = Matrix{Float64}(undef, n_samples, round(Int,sim_time_gt/t_steps_exp)+1)
                ANIMATED = false

                # Pre-allocate containers for batch file I/O (collect phase)
                sim_params_list = Vector{Dict}(undef, n_samples)
                pos2D_list = Vector{Any}(undef, n_samples)
                pos3D_list = Vector{Any}(undef, n_samples)
                fields_list = Vector{Any}(undef, n_samples)
                borderPts2D_list = Vector{Any}(undef, n_samples)
                η_steps_list = Vector{Vector}(undef, n_samples)
                β_steps_list = Vector{Vector}(undef, n_samples)
                cost_steps_list = Vector{Vector}(undef, n_samples)
                iter_steps_list = Vector{Vector}(undef, n_samples)

                set_file(joinpath(exp_path,"plots"))
                for n::Int in 1:n_samples
                    obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

                    stats = fit_model(model, scene, conditions, obs_border_pt_lst, θ, outliers=outlier_frames)
                                        
                    iterList = stats["iterList"]
                    costList = stats["cost_list"]
                    ηpList = stats["ηList"]
                    βpList = stats["βList"]

                    η = stats["η"]
                    β = stats["β"]
                    
                    reset_model!(est_model)
                    est_model.η = [η]
                    est_scene.β = [β]
                    
                    est_μ_list, gradList, borderPts2DList, fields, surface_pts_3D_gt, pos2D, pos3D, _, _, _ = simulate(est_model, est_scene, conditions)

                    est_h = get_height(est_μ_list, h)
                    
                    η_accuracy = (1-abs((η_gt[1]-η)/η_gt[1]))*100
                    β_accuracy = (1-abs((β_gt[1]-β)/β_gt[1]))*100
                    printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                    printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)
                    
                    η_pred[n] = η
                    β_pred[n] = β
                    est_h_list[n,:] = est_h
                    costnList[n] = costList
                    iternList[n] = iterList

                    # Collect phase: Store all data in memory instead of writing immediately
                    sim_params = Dict("gt_η" => η_gt,
                    "gt_β" => β_gt,
                    "η" => η,
                    "β" => β,
                    "η_accuracy" => η_accuracy,
                    "β_accuracy" => β_accuracy)
                    
                    sim_params_list[n] = sim_params
                    cost_steps_list[n] = costList
                    η_steps_list[n] = ηpList
                    β_steps_list[n] = βpList
                    iter_steps_list[n] = iterList
                    pos2D_list[n] = pos2D
                    pos3D_list[n] = pos3D
                    fields_list[n] = fields
                    borderPts2D_list[n] = borderPts2DList

                end

                # Write phase: Batch write all collected data at once (after loop completes)
                @info "Beginning batch file write phase ($n_samples runs)..."
                for n::Int in 1:n_samples
                    write_json(joinpath(exp_path,"data","stats","run_$n"), sim_params_list[n])
                    write_csv(joinpath(exp_path,"data","opt_data","cost_steps","run_$n"), cost_steps_list[n])
                    write_csv(joinpath(exp_path,"data","opt_data","eta_steps","run_$n"), η_steps_list[n])
                    write_csv(joinpath(exp_path,"data","opt_data","beta_steps","run_$n"), β_steps_list[n])
                    write_csv(joinpath(exp_path,"data","opt_data","iter","run_$n"), iter_steps_list[n])
                    write_csv(joinpath(exp_path,"data","opt_data","est_height","run_$n"), est_h_list[n,:])
                    write_2d_data(joinpath(exp_path,"data","sim_data","2D_points","run_$n"), pos2D_list[n])
                    write_data(joinpath(exp_path,"data","sim_data","3D_points","run_$n"), pos3D_list[n])
                    write_data(joinpath(exp_path,"data","sim_data","motion_fields ","run_$n"), fields_list[n])
                    write_2d_data(joinpath(exp_path,"data","sim_data","2D_border_points","run_$n"), borderPts2D_list[n])
                end

                write_csv(joinpath(exp_path,"data","eta_est"), η_pred)
                write_csv(joinpath(exp_path,"data","beta_est"), β_pred)
                write_csv(joinpath(exp_path,"data","h_est"), est_h_list)

                @info "Data writing completed. Results saved to $exp_path"
            end

        elseif gt_viscosity_type == "bulk_viscosity"
            
            av_η::Float64 = 0.0
            obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
            mode::String = exp_params["mode"] # "single_window" or "multiple_window" or "full_time"
            sim_time_window::Float64 = 30.0 # time window size for optimization in seconds
            model_gt.η = η_gt

            viscosity_type = "constant"
            
            conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, viewing_angles=[z_angle])
            model, scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F, control, viscosity_type, sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)

            set_file(joinpath(exp_path,"plots"))
            
            time_windows, windows, data_ranges_, t_windows = set_time_window(1/t_steps_exp, obs_border_pt_lst, method="linear", window_size=sim_time_exp)
            _, splinexObs_win, _, _ = set_time_window(1/t_steps_exp, splinexObs, method="linear", window_size=sim_time_exp)
            _, splineyObs_win, _, _ = set_time_window(1/t_steps_exp, splineyObs, method="linear", window_size=sim_time_exp)
            println("Time windows: $(time_windows)")
            obs_time = sum(time_windows)

            if obs_time < sim_time_gt
                @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time_gt, switching to observation time frame"
                sim_time_gt = obs_time
            end

            if obs_time < sim_time_exp
                @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                sim_time_exp = obs_time
            end
            
            data_pt_len = round(Int,obs_time/t_steps_exp)
            est_ηpList = Vector{Float64}(undef,data_pt_len)
            avg_ηList = Vector{Float64}(undef,data_pt_len)
            est_βpList = Vector{Float64}(undef,data_pt_len)
            cost_list = AbstractArray[]
            pred_h_list = AbstractArray[]

            if mode == "single_window"
                @info "Optimizing over a single time window"
                ti = 1
                data_range_ = data_ranges_[ti]
                scene.sim_time = time_windows[ti]
                if data_type == "physical"
                    _F = -F_ext*ones(Float64, round(Int, scene.sim_time*frame_rate)) # force applied to the cylinder in N
                    scene.cParam = _F
                else
                    scene.cParam = F[data_range_]
                end
                @info "Time window $(time_windows[ti])"

                println("Data frame : $(data_range_)")
                println("Time frame : $(scene.sim_time)")

                printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)

                obs_border_pt_lst_t = windows[ti] # align the observation points with the simulation time
                splinexObs_t = splinexObs_win[ti]
                splineyObs_t = splineyObs_win[ti]

                println("observation Window size: $(size(obs_border_pt_lst_t,1)) seconds")

                println("Time frame : $data_range_")
                printstyled("Time window: $(ti)\n"; color = :green)

                if data_type != "physical"
                    η_gt_ = model_gt.η[data_range_]
                    av_η = mean(η_gt_)
                    avg_ηList[data_range_] .= av_η
                    printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
                end

                stats = fit_model(model, scene, conditions, obs_border_pt_lst_t, θ, outliers=outlier_frames)

                est_ηpList[data_range_] .= stats["η"]
                est_βpList[data_range_] .= stats["β"]

                θ[1] = stats["η"]
                θ[2] = stats["β"]

                update_model!(model)
        
                iterList = stats["iterList"]
                costList = stats["cost_list"]
                ηpList = stats["ηList"]
                βpList = stats["βList"]

                viscosity_type = "bulk_viscosity"
                est_model, est_scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F[data_range_], control, viscosity_type, sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)
                est_model.η = est_ηpList[data_range_] 
                est_scene.β = est_βpList[data_range_]
                est_μ_list, gradList, borderPts2DList, fields_est, _, pos2D_est, pos3D_est, _, _, _ = simulate(est_model, est_scene, conditions)
                est_h_list = get_height(est_μ_list, h)

            else
                println("Number of time windows: $(length(windows))")

                for ti::Int in 1:length(windows)
                    data_range_ = data_ranges_[ti]
                    data_range_prev = ti > 1 ? data_ranges_[ti-1] : 1:data_range_[1]

                    scene.sim_time = time_windows[ti]
                    if data_type == "physical"
                        _F = -F_ext*ones(Float64, round(Int, scene.sim_time*frame_rate)) # force applied to the cylinder in N
                        scene.cParam = _F
                    else
                        scene.cParam = F[data_range_]
                    end
                    obs_border_pt_lst_t = windows[ti] # align the observation points with the simulation time

                    printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
                    println("Data frame : $(data_range_)")
                    println("Time frame : $(scene.sim_time)")
                    @info "Time window $(t_windows[ti])"
                    println("observation Window size: $(size(obs_border_pt_lst_t,1)) seconds")

                    if data_type != "physical" && viscosity_model != "carreau"
                        η_gt_ = model_gt.η[data_range_]
                        av_η = mean(η_gt_)
                        avg_ηList[data_range_] .= av_η
                        printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
                    end
                    
                    if data_type == "physical" && ti == 2
                        println("Second time window, applying correction factor to the parameters for better convergence...")
                        θ[1] = stats["η"]*4
                        θ[2] = stats["β"]
                    end

                    @info "Fitting model in time window $(ti)..."
                    stats = fit_model(model, scene, conditions, obs_border_pt_lst_t, θ, outliers=outlier_frames)

                    est_ηpList[data_range_] .= stats["η"]
                    est_βpList[data_range_] .= stats["β"]
                    push!(cost_list, stats["cost_list"])

                    θ[1] = stats["η"]
                    θ[2] = stats["β"]

                    update_model!(model)
                end

                @info "Completed all time windows."
                
                viscosity_type = "bulk_viscosity"
                est_model, est_scene = def_problem(geom_exp, ne_exp, η_start, _fem..., β_start, F, control, viscosity_type, sim_time_gt, t_steps_gt; viscosity_model=viscosity_model, _mesh_path_kw(exp_params)...)
                est_model.η = est_ηpList
                est_scene.β = est_βpList

                write_csv(joinpath(exp_path,"data", "est_η"), est_ηpList)
                write_csv(joinpath(exp_path,"data", "est_β"), est_βpList)

                write_csv(joinpath(exp_path,"data", "avg_η"), avg_ηList)
                write_csv(joinpath(exp_path,"data", "window_data","time_windows"), time_windows)
                write_csv(joinpath(exp_path,"data", "window_data","t_windows"), t_windows)
                write_csv(joinpath(exp_path,"data", "window_data","data_ranges"), data_ranges_)
                write_csv(joinpath(exp_path,"data", "window_data","windows_sizes"), windows)
                write_csv(joinpath(exp_path,"data", "window_data","cost_windows"), cost_list)
                
                est_μ_list, gradList, simBorderPts, fields_est, surface_pts_3D_est, pos2D_est, pos3D_est, _, _, _ = simulate(est_model, est_scene, conditions)
                est_h_list = get_height(est_μ_list, h)

                if data_type != "physical" && viscosity_model != "carreau"
                    gt_μ_list, gradList, borderPts2DList_gt, fields_gt, surface_pts_3D_gt, pos2D_gt, pos3D_gt, _, _, _ = simulate(model_gt, scene_exp, conditions)
                    gt_h_list = get_height(gt_μ_list, h)
                    write_csv(joinpath(exp_path,"data", "η_gt"), model_gt.η)
                    write_csv(joinpath(exp_path,"data", "β_gt"), β_gt)
                    write_csv(joinpath(exp_path,"data", "gt_h"), gt_h_list)
                    
                    write_2d_data(joinpath(exp_path,"data", "sim_data","2D_surface_points_gt"), pos2D_gt)
                    write_2d_data(joinpath(exp_path,"data", "sim_data","2D_border_points_gt"), borderPts2DList_gt)

                    write_data(joinpath(exp_path,"data", "sim_data","3D_points_gt"), pos3D_gt)
                    write_data(joinpath(exp_path,"data", "sim_data","motion_fields_gt "), fields_gt)
                    write_data(joinpath(exp_path,"data", "sim_data","3D_surface_points_gt"), surface_pts_3D_gt)
                end
            end
            
            write_csv(joinpath(exp_path,"data", "est_h"), est_h_list)

            write_data(joinpath(exp_path,"data", "sim_data","3D_points_est"), pos3D_est)
            write_data(joinpath(exp_path,"data", "sim_data","motion_fields_est"), fields_est)
            write_data(joinpath(exp_path,"data", "sim_data","3D_surface_points_est"), surface_pts_3D_est)
            write_2d_data(joinpath(exp_path,"data", "sim_data","2D_surface_points_est"), pos2D_est)
            write_2d_data(joinpath(exp_path,"data", "sim_data","2D_border_points_est"), simBorderPts)
        end
        end_time = Dates.now()
        write_time_log(start_time, end_time, exp_params; dest_dir=joinpath(exp_path, "logs"))
        @info "Data writing completed. Results saved to $exp_path"
    end
end

"""
    compare_pt_clouds(pred_pts, gt_pts) -> (hausdorff_distances, chamfer_distances, closest_pt_distances)

Compute, frame by frame, the KD-tree-based Hausdorff, Chamfer, and
(one-directional) closest-point RMSE distances between matching pairs of
predicted and ground-truth point clouds.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: matching iterables of
  point-cloud arrays, one pair per frame.

# Returns
- `hausdorff_distances::Vector{Float64}`, `chamfer_distances::Vector{Float64}`,
  `closest_pt_distances::Vector{Float64}`: one distance per frame.
"""
function compare_pt_clouds(pred_pts::AbstractArray, gt_pts::AbstractArray)
    hausdorff_distances = Float64[]
    chamfer_distances = Float64[]
    closest_pt_distances = Float64[]
    for (sim_pts, gt_pts) in zip(pred_pts, gt_pts)
        hausdorff_dist = hausdorff_distance_kdtree(sim_pts, gt_pts)
        chamfer_dist = chamfer_distance_kdtree(sim_pts, gt_pts)
        closest_pt_dist = closest_point_distance_kdtree(sim_pts, gt_pts)
        push!(hausdorff_distances, hausdorff_dist)
        push!(chamfer_distances, chamfer_dist)
        push!(closest_pt_distances, closest_pt_dist)
    end
    return  hausdorff_distances, chamfer_distances, closest_pt_distances
end

"""
    hausdorff_distance(pred_pts, gt_pts) -> Float64

Brute-force (O(N·M)) Hausdorff distance between two point clouds: the largest
of the two directional sup-min distances between `pred_pts` and `gt_pts`.
Superseded by `hausdorff_distance_kdtree`; kept for reference/benchmarking.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `hausdorff_dist::Float64`: the Hausdorff distance between the two point clouds.
"""
function hausdorff_distance(pred_pts::AbstractArray, gt_pts::AbstractArray)
    sup_pred = maximum([minimum(sqrt.(sum((gt_pts .- permutedims(pred_pt)).^2, dims=2))) for pred_pt in eachrow(pred_pts)])
    sup_gt = maximum([minimum(sqrt.(sum((pred_pts .- permutedims(gt_pt)).^2, dims=2))) for gt_pt in eachrow(gt_pts)])
    hausdorff_dist = max(sup_pred, sup_gt)
    return hausdorff_dist
end

"""
    chamfer_distance(pred_pts, gt_pts) -> Float64

Brute-force (O(N·M)) Chamfer distance between two point clouds: the sum of
the two directional mean-min distances between `pred_pts` and `gt_pts`.
Superseded by `chamfer_distance_kdtree`; kept for reference/benchmarking.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `chamfer_dist::Float64`: the Chamfer distance between the two point clouds.
"""
function chamfer_distance(pred_pts::AbstractArray, gt_pts::AbstractArray)
    mean_pred = mean([minimum(sqrt.(sum((gt_pts .- permutedims(pred_pt)).^2, dims=2))) for pred_pt in eachrow(pred_pts)])
    mean_gt = mean([minimum(sqrt.(sum((pred_pts .- permutedims(gt_pt)).^2, dims=2))) for gt_pt in eachrow(gt_pts)])
    chamfer_dist = mean_pred + mean_gt
    return chamfer_dist
end

"""
    _as_points_by_rows(pts) -> AbstractArray

Normalize a point cloud to a `(n_points, n_dims)` layout. Upstream writers are
inconsistent about whether points are stored by row or by column (e.g.
`surface_nodes` vs `3D_points_*`), so orient defensively: the coordinate axis
(2 or 3) is always the smaller of the two dimensions for a real point cloud.

# Arguments
- `pts::AbstractArray`: point cloud, stored as either `(n_points, n_dims)` or
  `(n_dims, n_points)`.

# Returns
- `AbstractArray`: `pts`, transposed if necessary, oriented as `(n_points, n_dims)`.
"""
function _as_points_by_rows(pts::AbstractArray)
    return size(pts, 1) < size(pts, 2) ? permutedims(pts) : pts
end

"""
    _nn_min_dists(query_pts, ref_pts) -> Vector{Float64}

Nearest-neighbor distance from each row of `query_pts` to the closest row of
`ref_pts`, via a `KDTree` instead of brute-force O(N·M) broadcasting.

# Arguments
- `query_pts::AbstractArray`, `ref_pts::AbstractArray`: point clouds (row- or
  column-major; normalized internally via `_as_points_by_rows`).

# Returns
- `dists::Vector{Float64}`: nearest-neighbor distance for each point in `query_pts`.
"""
function _nn_min_dists(query_pts::AbstractArray, ref_pts::AbstractArray)
    tree = KDTree(permutedims(_as_points_by_rows(ref_pts)))
    _, dists = nn(tree, permutedims(_as_points_by_rows(query_pts)))
    return dists
end

"""
    hausdorff_distance_kdtree(pred_pts, gt_pts) -> Float64

KDTree-accelerated Hausdorff distance: the larger of the two directional
max nearest-neighbor distances between `pred_pts` and `gt_pts`.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `Float64`: the Hausdorff distance between the two point clouds.
"""
function hausdorff_distance_kdtree(pred_pts::AbstractArray, gt_pts::AbstractArray)
    d_pred_to_gt = _nn_min_dists(pred_pts, gt_pts)
    d_gt_to_pred = _nn_min_dists(gt_pts, pred_pts)
    return max(maximum(d_pred_to_gt), maximum(d_gt_to_pred))
end

"""
    chamfer_distance_kdtree(pred_pts, gt_pts) -> Float64

KDTree-accelerated Chamfer distance: the sum of the two directional mean
nearest-neighbor distances between `pred_pts` and `gt_pts`.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `Float64`: the Chamfer distance between the two point clouds.
"""
function chamfer_distance_kdtree(pred_pts::AbstractArray, gt_pts::AbstractArray)
    d_pred_to_gt = _nn_min_dists(pred_pts, gt_pts)
    d_gt_to_pred = _nn_min_dists(gt_pts, pred_pts)
    return mean(d_pred_to_gt) + mean(d_gt_to_pred)
end

"""
    closest_point_distance_kdtree(pred_pts, gt_pts) -> Float64

RMSE of the `pred_pts` → `gt_pts` nearest-neighbor distances (one-directional,
unlike the symmetric `chamfer_distance_kdtree`/`hausdorff_distance_kdtree`).

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `Float64`: RMSE of the `pred_pts` → `gt_pts` nearest-neighbor distances.
"""
function closest_point_distance_kdtree(pred_pts::AbstractArray, gt_pts::AbstractArray)
    d_pred_to_gt = _nn_min_dists(pred_pts, gt_pts)
    return sqrt(mean(d_pred_to_gt.^2))
end

"""
    predict(filepath, filepath_gt)

Out-of-sample forward prediction for bulk-viscosity experiments already fit by
`optimize`/`replot`-style time windowing. For every leaf
directory under `filepath` (excluding `view_1`, walking
`elem_size/simtime/noise/dt/view/window`), re-simulate each time window using
the η/β estimated for that window (and, from the second window onward, also
simulate the *previous* window's estimate forward for comparison), then write
the resulting predicted height, 3D points, motion fields, and 2D
border/surface points back into that experiment's `data` directory.

Only applies when `sim_params["viscosity_type"] == "bulk_viscosity"`; a no-op
otherwise.

# Arguments
- `filepath::String`: root of the result tree to walk and update in place.
- `filepath_gt::String`: ground-truth data root (supplies `sim_params`).

# Returns
None. All outputs are side effects written to disk.
"""
function predict(filepath, filepath_gt)

    sim_params = read_json(joinpath(filepath_gt,"data","sim_params")) 

    r::Float64 = sim_params["r"]
    h::Float64 = sim_params["h"]
    
    viscosity_type::String = sim_params["viscosity_type"]
    F = Array(float.(sim_params["cParam"]))

    sim_time::Float64 = sim_params["simulation_time"]
    t_steps::Float64 = sim_params["time_steps"]
    
    camera_matrix::AbstractArray = reshape(Array(float.(sim_params["camera_matrix"])),3,3)
    _obj_pose::AbstractArray = sim_params["obj_pose"]
    obj_pose::Vector{Float64} = if size(_obj_pose, 1) == 4 && size(_obj_pose, 2) == 4
                                    _obj_pose[1:3,4]
                                else
                                    _obj_pose
                                end
    control::String = sim_params["control_type"]
    edge_radius::Union{Float64, Nothing} = get(sim_params, "edge_radius", nothing)

    viscosity_model::String = ""
    if viscosity_type == "bulk_viscosity" && haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
        viscosity_model = sim_params["model_type"]
    end

    ndim::Int = 3
    nDof_p::Int = 1
    nDof_u::Int = ndim

    geometry = get(sim_params, "geometry", :cylinder)
    _mesh_path_kw(d) = haskey(d, "mesh_path") ? (mesh_path=d["mesh_path"],) : (;)
    geom = if geometry === :cube
        @info "Using cuboid geometry with dimensions $r x $r x $h"
        lx = Float64(get(sim_params, "lx", r))
        ly = Float64(get(sim_params, "ly", r))
        lz = Float64(get(sim_params, "lz", h))
        Cuboid(lx, ly, lz, edge_radius)
    elseif geometry === :cylinder
        @info "Using cylinder geometry with radius $r and height $h"
        Cylinder(r, h)
    else
        error("Unsupported geometry type: $geometry")
    end

    gt_h_::Matrix{Float64} = Matrix{Float64}(undef,0,0)

    if viscosity_type == "bulk_viscosity"
        viscosity_type = "constant"

        elem_size_folders = readdir(joinpath(filepath))
        for elem_size_folder in elem_size_folders
            if elem_size_folder == "post_analysis"
                continue
            end

            sim_time_folders = readdir(joinpath(filepath, elem_size_folder))

            for sim_time_folder in sim_time_folders
                if !startswith(sim_time_folder, "simtime_")
                    continue
                end

                noise_folders = readdir(joinpath(filepath, elem_size_folder, sim_time_folder))
                for noise_folder in noise_folders
                    if !startswith(noise_folder, "noise_")
                        continue
                    end
                    step_path = joinpath(filepath, elem_size_folder, sim_time_folder, noise_folder)
                    step_folders = readdir(step_path)

                    for step_folder in step_folders
                        if !startswith(step_folder, "dt_")
                            continue
                        end

                        step_path = joinpath(step_path, step_folder)
                        view_folders = readdir(step_path)
                        for view_folder in view_folders
                            println("Processing view folder: $view_folder")
                            if !startswith(view_folder, "view_")
                                continue
                            end
                            view_path = joinpath(step_path, view_folder)

                            window_dirs = readdir(view_path)
                            for window_dir in window_dirs
    
                                win_exp_path = joinpath(view_path, window_dir)

                                println("Processing window: $win_exp_path")
                                exp_params = read_json(joinpath(win_exp_path ,"data","experiment_parameters"))

                                data_type::String = exp_params["data_type"]
                                ne_exp::Int = exp_params["ne_exp"]
                                z_angle::Float64 = exp_params["z_angle"]
                                element_shape_u = Symbol(exp_params["element_shape_u"])
                                basis_order_u   = Int(exp_params["basis_order_u"])
                                element_shape_p = Symbol(exp_params["element_shape_p"])
                                basis_order_p   = Int(exp_params["basis_order_p"])
                                element_shape_x = Symbol(exp_params["element_shape_x"])
                                basis_order_x   = Int(exp_params["basis_order_x"])

                                est_ηpList = readdlm(joinpath(win_exp_path,"data","est_η.csv"), ',', Float64)
                                est_βpList = readdlm(joinpath(win_exp_path,"data","est_β.csv"), ',', Float64)
                                est_h_list = readdlm(joinpath(win_exp_path,"data","est_h.csv"), ',', Float64)
                                
                                h_pred_list = AbstractArray[]
                                borderPts_pred_list = AbstractArray[]
                                fields_list_pred_list = AbstractArray[]
                                pos2D_pred_list = AbstractArray[]
                                pos3D_pred_list = AbstractArray[]
                                surface_pts_3D_pred_list = AbstractArray[]

                                data_ranges_ = get_time_windows(joinpath(win_exp_path,"data","window_data","data_ranges.csv"))
                                t_windows = readdlm(joinpath(win_exp_path,"data","window_data","t_windows.csv"),',',Float64)
                                time_windows = readdlm(joinpath(win_exp_path,"data","window_data","time_windows.csv"),',',Float64)
                                
                                conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, viewing_angles=[z_angle])
                                _fem = (element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x)
                                model, scene = def_problem(geom, ne_exp, 1.0, _fem..., 1.0, F, control, viscosity_type, sim_time, t_steps; _mesh_path_kw(exp_params)...)
                                _h = h
                                for ti::Int in 1:(size(data_ranges_, 1))

                                    data_range_ = data_ranges_[ti]
                                    scene.sim_time = time_windows[ti]
                                    height_range = collect(Float64, range(start=data_range_[1], stop=data_range_[end]+1, step=1))

                                    η = est_ηpList[data_range_][1] # get the first value in the estimation window the rest are the same for the window
                                    β = est_βpList[data_range_][1] # get the first value in the estimation window the rest are the same for the window
                                    scene.cParam = F[data_range_]

                                    printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
                                    println("Data frame : $(data_range_)")
                                    println("Time frame : $(scene.sim_time)")
                                    @info "Time window $(t_windows[ti])"

                                    if ti > 1
                                        @info "Predicting dynamics in time window $(ti)..."
                                        reset_model!(model)
                                        data_range_prev = data_ranges_[ti-1]
                                        η_pred = est_ηpList[data_range_prev][1] # get the first value in the estimation window the rest are the same for the window
                                        β_pred = est_βpList[data_range_prev][1] # get the first value in the estimation window the rest are the same for the window

                                        model.η = [η_pred]
                                        scene.β = [β_pred]
                                        pred_μ_list, gradList, pred_simBorderPts, fields_pred, surface_pts_3D_pred, pos2D_pred, pos3D_pred, _, _, _  = simulate(model, scene, conditions)
                                        pred_h = get_height(pred_μ_list, _h)

                                        push!(h_pred_list, pred_h)
                                        push!(borderPts_pred_list, pred_simBorderPts)
                                        push!(fields_list_pred_list, fields_pred)
                                        push!(pos2D_pred_list, pos2D_pred)
                                        push!(pos3D_pred_list, pos3D_pred)
                                        push!(surface_pts_3D_pred_list, surface_pts_3D_pred)
                                    end

                                    reset_model!(model)
                                    model.η = [η]
                                    scene.β = [β]
                                    est_μ_list, gradList, est_simBorderPts, fields_est, surface_pts_3D_est, pos2D_est, pos3D_est, _, _, _ = simulate(model, scene, conditions)

                                    h_est = get_height(est_μ_list, _h)
                                    if ti == 1
                                        push!(h_pred_list, h_est)
                                        push!(borderPts_pred_list, est_simBorderPts)
                                        push!(fields_list_pred_list, fields_est)
                                        push!(pos2D_pred_list, pos2D_est)
                                        push!(pos3D_pred_list, pos3D_est)
                                        push!(surface_pts_3D_pred_list, surface_pts_3D_est)
                                    end
                                    update_model!(model)
                                    _h = h_est[end] # update the height for the next window simulation
                                end

                                h_pred_vec               = reduce(vcat, h_pred_list)
                                pos3D_pred_vec           = reduce(vcat, pos3D_pred_list)
                                surface_pts_3D_pred_vec  = reduce(vcat, surface_pts_3D_pred_list)
                                fields_pred_vec          = reduce(vcat, fields_list_pred_list)
                                pos2D_pred_vec           = reduce(vcat, pos2D_pred_list)
                                borderPts_pred_vec       = reduce(vcat, borderPts_pred_list)

                                write_csv(joinpath(win_exp_path,"data", "pred_h"), h_pred_vec)
                                write_data(joinpath(win_exp_path,"data", "sim_data","3D_points_pred"), pos3D_pred_vec)
                                write_data(joinpath(win_exp_path,"data", "sim_data","3D_surface_points_pred"), surface_pts_3D_pred_vec)
                                write_data(joinpath(win_exp_path,"data", "sim_data","motion_fields_pred "), fields_pred_vec)
                                write_2d_data(joinpath(win_exp_path,"data", "sim_data","2D_points_pred"), pos2D_pred_vec)
                                write_2d_data(joinpath(win_exp_path,"data", "sim_data","2D_border_points_pred"), borderPts_pred_vec)
                            end
                        end
                    end
                end
            end
        end
    end
    return
end

"""
    ExpLeaf

One `dt_.../view_...` result directory (`exp_path`) within an `ExpGroup`.
"""
struct ExpLeaf
    step_folder::String
    view_folder::String
    exp_path::String
end

"""
    ExpGroup

All `ExpLeaf`s sharing the same `(elem_size_folder, sim_time_folder,
noise_folder)` combination, as produced by `collect_experiment_groups`.
"""
struct ExpGroup
    elem_size_folder::String
    sim_time_folder::String
    noise_folder::String
    step_path::String
    leaves::Vector{ExpLeaf}
end

"""
    collect_experiment_groups(filepath) -> Vector{ExpGroup}

Walk the `elem_size/simtime/noise/dt/view` result folder tree under `filepath`
once and group leaves by `(elem_size_folder, sim_time_folder, noise_folder)`,
which is the level at which `replot` resets its cross-view comparison
data.

# Arguments
- `filepath::String`: root of the result tree to walk.

# Returns
- `Vector{ExpGroup}`: one group per `(elem_size_folder, sim_time_folder, noise_folder)` combination.
"""
function collect_experiment_groups(filepath)
    groups = ExpGroup[]
    for elem_size_folder in readdir(filepath)
        elem_size_folder == "post_analysis" && continue

        for sim_time_folder in readdir(joinpath(filepath, elem_size_folder))
            startswith(sim_time_folder, "simtime_") || continue
            noise_base = joinpath(filepath, elem_size_folder, sim_time_folder)

            for noise_folder in readdir(noise_base)
                startswith(noise_folder, "noise_") || continue
                step_path = joinpath(noise_base, noise_folder)

                leaves = ExpLeaf[]
                for step_folder in readdir(step_path)
                    startswith(step_folder, "dt_") || continue
                    dt_path = joinpath(step_path, step_folder)

                    for view_folder in readdir(dt_path)
                        startswith(view_folder, "view_") || continue
                        push!(leaves, ExpLeaf(step_folder, view_folder, joinpath(dt_path, view_folder)))
                    end
                end
                push!(groups, ExpGroup(elem_size_folder, sim_time_folder, noise_folder, step_path, leaves))
            end
        end
    end
    return groups
end

"""
    replot(filepath, filepath_gt)

Regenerate all comparison plots (and some derived post-processed data) for
every fitted experiment found under `filepath`, using the ground-truth
`sim_params`/height/contour data in `filepath_gt`. Does not re-run any
optimization; it only reads previously written `optimize`/`predict` outputs
and produces plots/derived data from them. Iterates leaves grouped by
`collect_experiment_groups`.

For `viscosity_type == "constant"`: plots the observed-vs-simulated contour,
closest-point distance error, estimated-vs-ground-truth height (and its
error), η/β vs. iteration, and cost vs. iteration. When `noise_level == 0`,
additionally post-processes the stored ``(η, β)`` cost-surface grid — locating
its minimum, computing the local Hessian and its principal (steepest/flattest)
directions, and plotting cost-surface slices along them — and writes the
result to `data/direction_analysis`/`data/slice_data`. When `noise_level != 0`
(repeated noisy-sample runs), instead plots the covariance ellipse and
error-banded height/relative-error/normalized-error series across samples.

For `viscosity_type == "bulk_viscosity"`: for each time-windowed
sub-experiment, plots per-window estimated vs. ground-truth η, β, and their
ratio/product/sum, height and height error/normalized error, closest-point
distance error, and (via `compare_pt_clouds`) Hausdorff/Chamfer/
closest-point 3D surface distances — each comparing the within-window
estimate against the forward-predicted trajectory from the previous window,
with vertical lines marking window boundaries.

# Arguments
- `filepath::String`: root of the result tree to walk and regenerate plots for.
- `filepath_gt::String`: ground-truth data root (supplies `sim_params`, height,
  and contour data).

# Returns
None. All outputs are side effects (PDFs under each experiment's `plots/`
directory, JSON/CSV under `data/`).
"""
function replot(filepath, filepath_gt)
    sim_params = read_json(joinpath(filepath_gt,"data","sim_params"))

    r::Float64 = sim_params["r"]
    h::Float64 = sim_params["h"]


    viscosity_type::String = sim_params["viscosity_type"]
    F = Array(float.(sim_params["cParam"]))

    sim_time::Float64 = sim_params["simulation_time"]
    t_steps::Float64 = sim_params["time_steps"]
    
    camera_matrix::Matrix{Float64} = reshape(Array(float.(sim_params["camera_matrix"])),3,3)
    _obj_pose::AbstractArray = sim_params["obj_pose"]
    obj_pose::Vector{Float64} = if size(_obj_pose, 1) == 4 && size(_obj_pose, 2) == 4
                                    _obj_pose[1:3,4]
                                else
                                    _obj_pose
                                end

    control::String = sim_params["control_type"]

    viscosity_model::String = ""
    if viscosity_type == "bulk_viscosity" && haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
            viscosity_model = sim_params["model_type"]
    end

    ndim::Int = 3
    nDof_p::Int = 1
    nDof_u::Int = ndim

    gt_h_::Matrix{Float64} = Matrix{Float64}(undef,0,0)

    for group in collect_experiment_groups(filepath)
        for leaf in group.leaves
            view_folder = leaf.view_folder
   
            exp_path = leaf.exp_path
            if viscosity_type == "constant"
                println("Comparing experiments in: $exp_path")
                
                η_gt = float(sim_params["η"][1])
                β_gt = float(sim_params["β"][1])

                exp_params = read_json(joinpath(exp_path,"data","experiment_parameters"))

                noise_level = exp_params["noise_level"]
                sim_time_exp = exp_params["sim_time_exp"]   
                ne_exp = exp_params["ne_exp"]
                data_type = exp_params["data_type"]
                control = exp_params["control"]
                num_exp_points::Int = round(Int, 20/t_steps)

                printstyled(num_exp_points, " experimental points will be used for optimization\n"; color = :blue)

                if data_type != "physical"
                    gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
                end
                if data_type  == "synthetic"
                    ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                    @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                elseif data_type == "simulated"
                    ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
                    @info "Data type $data_type Reading simulated contour data from $(joinpath(filepath_gt,"data","sim_data","contour_data")) of $(length(ObsDataList)) time steps"
                else
                    error("Unknown data type: $data_type")
                end

                time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))
                gt_h = gt_h_[1:(round(Int, sim_time/ t_steps))+1]

                conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
                if noise_level == 0.0
                    obs_border_pt_lst, sim_border_pt_lst, nSplinex, nSpliney, splinex, spliney, _ = _get_borders(data_type, filepath_gt, exp_path, num_exp_points)
                    
                    est_η = readdlm(joinpath(exp_path,"data","η.csv"), ',', Float64)
                    est_β = readdlm(joinpath(exp_path,"data","β.csv"), ',', Float64)
                    est_h = readdlm(joinpath(exp_path,"data","est_h.csv"), ',', Float64)
                    cost_iter = readdlm(joinpath(exp_path,"data","cost_iter.csv"), ',', Float64)
                    
                    stats = read_json(joinpath(exp_path,"data","stats")) 
                    η = stats["η"]
                    β = stats["β"]
                    
                    d_est, pairs = closest_point(sim_border_pt_lst, obs_border_pt_lst)

                    if sim_time_exp == 5.0
                        cont_y_min = 350
                    elseif sim_time_exp == 10.0
                        cont_y_min = 400
                    elseif sim_time_exp == 20.0
                        cont_y_min = 420
                    else
                        cont_y_min = 480
                    end
                    
                    contour_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(contour_plt, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                    Plots.plot!(contour_plt, nSplinex[end], nSpliney[end], label="Ground truth contour", color=:red)
                    Plots.yflip!(contour_plt, true)
                    Plots.xlims!(contour_plt, 480, 1520)
                    Plots.ylims!(contour_plt, cont_y_min, 1200)
                    Plots.xlabel!(contour_plt, L"x\;\mathrm{[px]}")
                    Plots.ylabel!(contour_plt, L"y\;\mathrm{[px]}")
                    Plots.savefig(contour_plt, joinpath(exp_path,"plots","contour_comparison.pdf"))

                    contour_plt_zoom = set_plot(fs, sz=(1200, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(contour_plt_zoom, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                    Plots.plot!(contour_plt_zoom, nSplinex[end], nSpliney[end], label=false, color=:red)
                    Plots.xlims!(contour_plt_zoom, 1000, 1520)
                    Plots.ylims!(contour_plt_zoom, cont_y_min, 1200)
                    Plots.xticks!(contour_plt_zoom, 1100:200:1520)
                    Plots.yflip!(contour_plt_zoom, true)
                    Plots.xlabel!(contour_plt_zoom, L"x\;\mathrm{[px]}")
                    Plots.ylabel!(contour_plt_zoom, L"y\;\mathrm{[px]}")
                    Plots.savefig(contour_plt_zoom, joinpath(exp_path,"plots","contour_comparison_zoom.pdf"))

                
                    costList = stats["cost_list"]
                    iterList = stats["iterList"]

                    plt_cnt_error = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(plt_cnt_error, time[1:(length(d_est)-1)], d_est[2:end], label="Closest point distance error", legend=:outerbottom, legend_column=2)
                    Plots.xlabel!(plt_cnt_error, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_cnt_error, L"\mathrm{Closest\;Point\;Distance\;[px]}")
                    Plots.xlims!(plt_cnt_error, 0, end_obs_win)
                    Plots.ylims!(plt_cnt_error, 0, max(maximum(d_est[2:end])*1.1, 0.5))
                    Plots.savefig(plt_cnt_error, joinpath(exp_path,"plots","closest_point_distance_error.pdf"))
                    
                    est_h = vec(Float64.(collect(est_h)))
                    gt_h = vec(Float64.(collect(gt_h)))
                    
                    n_time = min(length(time), length(est_h), length(gt_h)) #, (end_obs_win/t_steps + 1))
                    if n_time < length(time) || n_time < length(est_h) || n_time < length(gt_h)
                        @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(length(est_h)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                    end
                    
                    h_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(h_plt, time[1:n_time], est_h[1:n_time], label=L"h_{\mathrm{est}}(t)", legend=:outerbottom, legend_column=2)
                    Plots.plot!(h_plt, time[1:n_time], gt_h[1:n_time], label=L"h_{\mathrm{gt}}(t)", legend=:outerbottom, legend_column=2)
                    Plots.xlabel!(h_plt, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(h_plt, L"\mathrm{Height\;[mm]}")
                    Plots.xlims!(h_plt, 0, end_obs_win)
                    Plots.savefig(h_plt, joinpath(exp_path,"plots","h_est.pdf"))

                    error_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(error_plt, time[1:n_time], abs.(est_h[1:n_time] .- gt_h[1:n_time]), label="Height estimation error", legend=:outerbottom, legend_column=1)
                    Plots.xlabel!(error_plt, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(error_plt, "Height Error [mm]")
                    Plots.xlims!(error_plt, 0, end_obs_win)
                    Plots.savefig(error_plt, joinpath(exp_path,"plots","h_est_error.pdf"))

                    η_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(η_plt, est_η, label="Estimated η", marker=1, legend=:outerbottom, legend_column=2)
                    Plots.hline!([η_gt], label="Ground truth η", legend=:outerbottom, legend_column=2)
                    Plots.xlabel!(η_plt, L"\mathrm{Iterations}")
                    Plots.ylabel!(η_plt, latexstring("\$\\eta\$ [kPa s]"))
                    Plots.savefig(η_plt, joinpath(exp_path,"plots","η.pdf"))
                    
                    β_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(β_plt, est_β, label="Estimated β", marker=1, legend=:outerbottom, legend_column=2)
                    Plots.hline!(β_plt, [β_gt], label="Ground truth β",legend=:outerbottom, legend_column=2)
                    Plots.xlabel!(β_plt, L"\mathrm{Iterations}")
                    Plots.ylabel!(β_plt, latexstring("\$\\beta\$ [MPa s m\$^{-1}\$]"))   
                    Plots.savefig(β_plt, joinpath(exp_path,"plots","β.pdf"))

                    cost_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(cost_plt, iterList, costList, label="Cost", marker=1, yscale=:log10, xminorgrid = :false,legend=:outerbottom, legend_column=1)
                    Plots.xlabel!(cost_plt, L"\mathrm{Iterations}")
                    Plots.ylabel!(cost_plt, L"\mathrm{Cost\;[px]}")
                    Plots.savefig(cost_plt, joinpath(exp_path,"plots","cost_steps.pdf"))

                    cost_log_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(cost_log_plt, iterList, costList, label="Cost", marker=1, yscale=:log10, xscale=:log10, legend=:outerbottom, legend_column=1)
                    Plots.xlabel!(cost_log_plt, L"\mathrm{Iterations}")
                    Plots.ylabel!(cost_log_plt, L"\mathrm{Cost\;[px]}")
                    Plots.savefig(cost_log_plt, joinpath(exp_path,"plots","cost_steps_log.pdf"))

                    # Coerce contour parameters to concrete arrays (be defensive against scalar values)
                    function _to_vector(x)
                        if isa(x, Number)
                            return Float64[x]
                        elseif isa(x, AbstractVector)
                            return Float64.(collect(x))
                        else
                            try
                                return Float64.(collect(x))
                            catch err
                                error("Cannot coerce input to vector: $err")
                            end
                        end
                    end
                    
                    function _to_matrix(x)
                        if isa(x, AbstractMatrix)
                            return Float64.(collect(x))
                        else
                            try
                                rows = collect(x)

                                # Handle a collection/set of row vectors by stacking to a matrix
                                if !isempty(rows)
                                    row_arrays = map(rows) do r
                                        if isa(r, AbstractVector)
                                            Float64.(collect(r))
                                        elseif isa(r, AbstractMatrix) && (size(r,1) == 1 || size(r,2) == 1)
                                            vec(Float64.(collect(r)))
                                        else
                                            nothing
                                        end
                                    end

                                    if all(!isnothing, row_arrays)
                                        row_arrays = Vector{Vector{Float64}}(row_arrays)
                                        ncols = length(row_arrays[1])
                                        if any(length(row) != ncols for row in row_arrays)
                                            error("Cannot coerce input to matrix: inconsistent row lengths")
                                        end
                                        return reduce(vcat, permutedims.(row_arrays))
                                    end
                                end

                                return Float64.(rows)
                            catch err
                                error("Cannot coerce input to matrix: $err")
                            end
                        end
                    end
                    
                    try 
                        contour_plot_params = read_json(joinpath(exp_path,"data","contour_plot_params")) 
                        ηList = _to_vector(contour_plot_params["η_list"])
                        βList = _to_vector(contour_plot_params["β_list"])
                        CostMat = _to_matrix(contour_plot_params["cost_mat"])

                        # --- Post-process: compute Hessian-based principal directions from CostMat ---
                            # Ensure CostMat shape matches (length(ηList), length(βList))
                            nx = length(ηList); ny = length(βList)
                            sz = size(CostMat)
                            if ndims(CostMat) == 1 && length(CostMat) == nx*ny
                                CostMat = reshape(CostMat, nx, ny)
                            elseif ndims(CostMat) == 2
                                if sz[1] == ny && sz[2] == nx
                                    # CostMat is (ny, nx) => transpose to (nx, ny)
                                    CostMat = CostMat'
                                elseif sz[1] == nx && sz[2] == ny
                                    # already in expected orientation
                                else
                                    @warn "CostMat shape $(sz) does not match (nx,ny)=($(nx),$(ny)). Attempting to reshape if possible."
                                    if prod(sz) == nx*ny
                                        CostMat = reshape(vec(CostMat), nx, ny)
                                    end
                                end
                            end
                            # Prepare interpolation grid (rows->β, cols->η)
                            Z_for_interp = CostMat'

                            # locate minimum in the cost matrix (η-index, β-index)
                            minval, minidx = findmin(CostMat)
                            i0 = Int(minidx[1]); j0 = Int(minidx[2])
                            # guard against out-of-range indices
                            ni = length(ηList); nj = length(βList)
                            if i0 < 1 || i0 > ni || j0 < 1 || j0 > nj
                                error("Minimum index from CostMat (i0=$i0, j0=$j0) is out of bounds for ηList/βList lengths (", ni, ",", nj, ")")
                            end
                            η0 = ηList[i0]; β0 = βList[j0]

                            # choose interior index for finite-difference Hessian computation
                            ii = clamp(i0, 2, max(2, ni-1))
                            jj = clamp(j0, 2, max(2, nj-1))
                            if ii != i0 || jj != j0
                                @warn "Minimum on boundary; using nearby interior index (ii,jj)=($(ii),$(jj)) for Hessian"
                            end

                            # grid spacings (assume near-uniform)
                            dη = (ni > 1) ? (ηList[min(2,ni)] - ηList[1]) : 1.0
                            dβ = (nj > 1) ? (βList[min(2,nj)] - βList[1]) : 1.0

                            C = CostMat
                            # centered finite differences for second derivatives
                            d2C_dη2 = (C[ii+1,jj] - 2C[ii,jj] + C[ii-1,jj]) / (dη^2)
                            d2C_dβ2 = (C[ii,jj+1] - 2C[ii,jj] + C[ii,jj-1]) / (dβ^2)
                            d2C_dηdβ = (C[ii+1,jj+1] - C[ii+1,jj-1] - C[ii-1,jj+1] + C[ii-1,jj-1]) / (4.0*dη*dβ)

                            H = [d2C_dη2 d2C_dηdβ; d2C_dηdβ d2C_dβ2]
                            ev = eigen(Symmetric(H))
                            evals = ev.values
                            evecs = ev.vectors

                            # steepest = largest eigenvalue, flattest = smallest
                            idx_max = argmax(evals); idx_min = argmin(evals)
                            v_steep = evecs[:, idx_max] / norm(evecs[:, idx_max])
                            v_flat  = evecs[:, idx_min] / norm(evecs[:, idx_min])

                            # param-space sampling along directions (normalized) centered at minimum
                            span_η = maximum(ηList) - minimum(ηList)
                            span_β = maximum(βList) - minimum(βList)
                            diag_span = sqrt(span_η^2 + span_β^2)
                            tvec = collect(range(-0.6*diag_span, stop=0.6*diag_span, length=401))

                            etas_steep = η0 .+ tvec .* v_steep[1]
                            betas_steep = β0 .+ tvec .* v_steep[2]
                            etas_flat = η0 .+ tvec .* v_flat[1]
                            betas_flat = β0 .+ tvec .* v_flat[2]

                            # keep samples inside the original grid extents
                            minη, maxη = minimum(ηList), maximum(ηList)
                            minβ, maxβ = minimum(βList), maximum(βList)
                            mask_steep = (etas_steep .>= minη) .& (etas_steep .<= maxη) .& (betas_steep .>= minβ) .& (betas_steep .<= maxβ)
                            mask_flat  = (etas_flat  .>= minη) .& (etas_flat  .<= maxη) .& (betas_flat  .>= minβ) .& (betas_flat  .<= maxβ)

                            t_steep = tvec[mask_steep]; t_flat = tvec[mask_flat]
                            etas_steep = etas_steep[mask_steep]; betas_steep = betas_steep[mask_steep]
                            etas_flat  = etas_flat[mask_flat];  betas_flat  = betas_flat[mask_flat]

                            # interpolate cost along these samples
                            zs_steep = interp_z_at(etas_steep, betas_steep, ηList, βList, Z_for_interp)
                            zs_flat  = interp_z_at(etas_flat,  betas_flat,  ηList, βList, Z_for_interp)

                            # Contour with overlaid directions and minimum (use est_η/est_β as iteration path)
                            # ensure estimation paths are 1D Float64 vectors
                            est_η = vec(Float64.(collect(est_η)))
                            est_β = vec(Float64.(collect(est_β)))

                            @info "Post-process grid sizes: len(ηList)=$(length(ηList)), len(βList)=$(length(βList)), CostMat_size=$(size(CostMat))"
                            # Compute adaptive clims from CostMat to represent all values
                            cost_min, cost_max = extrema(CostMat)
                            cost_max = min(cost_max, 1e3) # avoid zero max
                            cost_clims = (0, cost_max)
                            plt_dirs = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                            Plots.contour!(plt_dirs, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                            Plots.plot!(plt_dirs, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.plot!(plt_dirs, etas_steep, betas_steep, label = "Steepest dir", color=:black, legend=:outerbottom, legend_column=3)
                            Plots.plot!(plt_dirs, etas_flat,  betas_flat,  label = "Flattest dir",  lw=1, color=:magenta, legend=:outerbottom, legend_column=3)
                            Plots.scatter!(plt_dirs, [η0], [β0], label="Minimum Cost", ms=:5, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_dirs, [η_gt], [β_gt], label="Ground truth", ms=:5, m=:star5, color=def_red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.xlabel!(plt_dirs, L"\eta\;\mathrm{[kPa\, s]}")
                            Plots.ylabel!(plt_dirs, L"\beta\;\mathrm{[Pa\, s \, m]}")
                            Plots.savefig(plt_dirs, joinpath(exp_path, "plots", "cost_surface_with_directions.pdf"))
                            
                            plt_cont = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                            Plots.contour!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                            Plots.plot!(plt_cont, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=:5, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:5, m=:star5, color=def_red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.xlabel!(plt_cont, L"\eta\;\mathrm{[kPa\, s]}")
                            Plots.ylabel!(plt_cont, L"\beta\;\mathrm{[Pa\, s \, m]}")
                            Plots.savefig(plt_cont, joinpath(exp_path, "plots", "cost_surface.pdf"))

                            plt_cont = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                            Plots.contourf!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                            Plots.plot!(plt_cont, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=:5, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:5, m=:star5, color=def_red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.xlabel!(plt_cont, L"\eta\;\mathrm{[kPa\, s]}")
                            Plots.ylabel!(plt_cont, L"\beta\;\mathrm{[Pa\, s \, m]}")
                            Plots.savefig(plt_cont, joinpath(exp_path, "plots", "cost_surface_iter.pdf"))

                            # 2D slices: cost vs distance along the two directions
                            plt_slices = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                            if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                Plots.plot!(plt_slices, t_steep, zs_steep, label = "Steepest direction", color=:black, legend=:outerbottom, legend_column=2)
                            else
                                @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                            end
                            if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                Plots.plot!(plt_slices, t_flat,  zs_flat,  label = "Flattest direction",  lw=1, color=:gray, legend=:outerbottom, legend_column=2)
                            else
                                @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                            end
                            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum")
                            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;[px]}")
                            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
                            Plots.savefig(plt_slices, joinpath(exp_path, "plots", "cost_slices_along_directions.pdf"))

                            Plots.ylims!(plt_slices, 0, 10)
                            Plots.savefig(plt_slices, joinpath(exp_path, "plots", "cost_slices_along_directions_zoomed.pdf"))
                    
                            # save analysis metadata
                            dir_info = Dict("η_min"=>η0, "β_min"=>β0, "Hessian"=>H, "evals"=>evals, "v_steep"=>v_steep, "v_flat"=>v_flat)
                            slice_data = Dict(
                                "steep"=>Dict("t"=>t_steep, "etas"=>etas_steep, "betas"=>betas_steep, "zs"=>zs_steep),
                                "flat"=>Dict("t"=>t_flat, "etas"=>etas_flat, "betas"=>betas_flat, "zs"=>zs_flat)
                            )
                            write_json(joinpath(exp_path, "data", "direction_analysis"), dir_info)
                            write_json(joinpath(exp_path, "data", "slice_data"), slice_data)
                    catch err
                        @warn "Failed to post-process contour parameters and cost surface: $err"
                    end
                else
                    # try
                        η_pred = readdlm(joinpath(exp_path,"data","eta_est.csv"), ',', Float64) # estimated η values per sample [n x n_iter]
                        β_pred = readdlm(joinpath(exp_path,"data","beta_est.csv"), ',', Float64) # estimated β values per sample [n x n_iter]
                        h_pred = readdlm(joinpath(exp_path,"data","h_est.csv"), ',', Float64) # estimated height values per sample [n x sim_time]

                        exp_params = read_json(joinpath(exp_path,"data","experiment_parameters"))
                        
                        sim_time_exp = exp_params["sim_time_exp"]
                        obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noise_level)
                    
                        covarience_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                        plot_covariance!(covarience_plt, η_pred[:,1], β_pred[:,1])
                        Plots.xlabel!(covarience_plt, latexstring("\$\\eta\$ [kPa s]"))
                        Plots.ylabel!(covarience_plt, latexstring("\$\\beta\$ [MPa s m\$^{-1}\$]"))
                        Plots.savefig(covarience_plt, joinpath(exp_path,"plots","covariance.pdf"))

                        n_time = min(length(time), size(h_pred, 2), length(gt_h)) #, (end_obs_win/t_steps + 1))
                        if n_time < length(time) || n_time < size(h_pred, 2) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(size(h_pred, 2)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end

                        h_plot = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=2)
                        Plots.plot!(h_plot, time[1:n_time], gt_h[1:n_time], label=L"h_{\mathrm{gt}}(t)")
                        StatsPlots.errorline!(h_plot, time[1:n_time], h_pred[:,1:n_time], label=L"h_{\mathrm{est}}(t)")
                        Plots.xlabel!(h_plot, L"\mathrm{Time\;[s]}")
                        Plots.ylabel!(h_plot, L"\mathrm{Height\;[mm]}")
                        Plots.xlims!(h_plot, 0, end_obs_win)
                        Plots.savefig(h_plot, joinpath(exp_path,"plots","h_est_noisy.pdf"))
                        
                        # Add inset subplot for zoomed region (10-20 seconds)
                        t_min, t_max = 5.0, 5.1
                        mask_zoom = (time .>= t_min) .& (time .<= t_max)
                        if sum(mask_zoom) > 0
                            idx_zoom = findall(mask_zoom)
                            # Build a single combined figure with two subplots and one shared legend
                            # Use a 2-row layout: top row 2 subplots, bottom row small height for a centered legend
                            plt_layout = @layout [a b; c{0.12h}]

                            # Build a single combined figure and draw into subplots 1,2 and 3
                            plt_combined = set_subplot(fs, sz=(3500, 350), layout=plt_layout)

                            # Full time series on left subplot (subplot=1)
                            Plots.plot!(plt_combined, time, gt_h, subplot=1, label="", color=:blue)
                            try
                                StatsPlots.errorline!(plt_combined, time[1:n_time], h_pred[:,1:n_time], subplot=1, label="", color=def_red)
                            catch e
                                @debug "errorline! failed with h_pred as-is, retrying transposed: $e"
                                StatsPlots.errorline!(plt_combined, time[1:n_time], h_pred[:,1:n_time]', subplot=1, label="", color=def_red)
                            end
                            Plots.xlabel!(plt_combined, L"\mathrm{Time\;[s]}", subplot=1)
                            Plots.ylabel!(plt_combined, L"\mathrm{Height\;[mm]}", subplot=1)

                            # Zoomed region on right subplot (subplot=2)
                            Plots.plot!(plt_combined, time[idx_zoom], gt_h[idx_zoom], subplot=2, label="", color=:blue)
                            try
                                StatsPlots.errorline!(plt_combined, time[idx_zoom], h_pred[:, idx_zoom], subplot=2, label="", color=def_red)
                            catch e
                                @debug "errorline! failed with h_pred as-is, retrying without leading colon: $e"
                                StatsPlots.errorline!(plt_combined, time[idx_zoom], h_pred[idx_zoom], subplot=2, label="", color=def_red)
                            end
                            Plots.xlabel!(plt_combined, L"\mathrm{Time\;[s]}", subplot=2)
                            Plots.xlims!(plt_combined, t_min, t_max, subplot=2)

                            # Bottom subplot (subplot=3) — draw a single centered annotation as deterministic legend
                            Plots.plot!(plt_combined, [], [], subplot=3, label=L"h_{\mathrm{gt}}(t)", framestyle=:none, legend=:outerbottom, color=:blue, legend_column=2, background_color=:transparent)
                            Plots.plot!(plt_combined, [], [], subplot=3, label=L"h_{\mathrm{est}}(t)", framestyle=:none, legend=:outerbottom, color=def_red, background_color=:transparent)
                            Plots.xlims!(plt_combined, 0.0, 1.0, subplot=3)
                            Plots.ylims!(plt_combined, 0.0, 1.0, subplot=3)

                            # Save combined figure
                            Plots.savefig(plt_combined, joinpath(exp_path, "plots", "h_est_noisy_zoomed.pdf"))
                        end

                        error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time]) ./ gt_h[1:n_time]
                        h_error_plot = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                        StatsPlots.errorline!(h_error_plot, time[1:n_time], error', label="Height estimation error", legend=:outerbottom, legend_column=1)
                        Plots.xlabel!(h_error_plot, L"\mathrm{Time\;[s]}")
                        Plots.ylabel!(h_error_plot, L"\mathrm{Relative\;Height\;Error}")
                        Plots.xlims!(h_error_plot, 0, end_obs_win)
                        Plots.savefig(h_error_plot, joinpath(exp_path,"plots","h_rel_error_noisy.pdf"))

                        h_norm = h_pred[:,1:n_time]' ./ gt_h[1:n_time]
                        h_normalized_plot = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                        StatsPlots.errorline!(h_normalized_plot, time[1:n_time], h_norm', label="Normalized height estimation error", legend=:outerbottom, legend_column=1)
                        Plots.xlabel!(h_normalized_plot, L"\mathrm{Time\;[s]}")
                        Plots.ylabel!(h_normalized_plot, L"h_{\mathrm{est}}/h_{\mathrm{gt}}")
                        Plots.xlims!(h_normalized_plot, 0, end_obs_win)
                        Plots.savefig(h_normalized_plot, joinpath(exp_path,"plots","h_normalized_rel_error_noisy.pdf"))
                    # catch err
                    #     @error "Replotting for noisy data with noise $noise_level failed: $err"
                    # end
                end
            elseif viscosity_type == "bulk_viscosity"
                local η_gt = Vector{Float64}(undef, 0)
                local β_gt = Vector{Float64}(undef, 0)
                
                window_dirs = readdir(exp_path)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        println("Skipping directory: $window_dir")
                        continue
                    end
                    win_exp_path = joinpath(exp_path, window_dir)
                    
                    println("Processing window: $win_exp_path")
                    exp_params = read_json(joinpath(win_exp_path ,"data","experiment_parameters"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    ne = exp_params["ne_exp"]
                    
                    gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
                    if data_type != "physical"
                        noise_level = exp_params["noise_level"]
                    else
                        noise_level = 0.0
                    end

                    est_ηpList = readdlm(joinpath(win_exp_path,"data","est_η.csv"), ',', Float64)
                    est_βpList = readdlm(joinpath(win_exp_path,"data","est_β.csv"), ',', Float64)
                    cost_init = readdlm(joinpath(win_exp_path,"data","window_data","cost_windows.csv"), ',', '\n')[1,1]

                    est_h_list = readdlm(joinpath(win_exp_path,"data","est_h.csv"), ',', Float64)
                    pred_h_list = readdlm(joinpath(win_exp_path,"data","pred_h.csv"))
                    ratio_gt = similar(est_ηpList)
                    ratio_opt = similar(est_ηpList)
                    avg_ηList = similar(est_ηpList)

                    data_ranges_ = get_time_windows(joinpath(win_exp_path,"data","window_data","data_ranges.csv"))
                    t_windows = readdlm(joinpath(win_exp_path,"data","window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(joinpath(win_exp_path,"data","window_data","time_windows.csv"),',',Float64)

                    println("Time windows: $(time_windows)")
                    obs_time = sum(time_windows)
                    effective_sim_time = sim_time

                    if obs_time != effective_sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $effective_sim_time, updating time frame"
                        effective_sim_time = min(obs_time, effective_sim_time)
                        obs_time = effective_sim_time
                    end

                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    obs_border_pt_lst, sim_border_pt_lst, nSplinex, nSpliney, splinex, spliney = _get_borders(data_type, filepath_gt, win_exp_path, data_point_len+1; view_folder=view_folder)
                    pred_border_pt_lst, _, _ = read_csv(joinpath(win_exp_path,"data","sim_data","view_1","2D_border_points_pred"))

                    ratio_opt = est_ηpList ./ est_βpList
                    product_opt = est_ηpList .* est_βpList
                    sum_opt = est_ηpList .+ est_βpList

                    if data_type != "physical"
                        η_gt = float.(sim_params["η"])
                        β_gt = float.(sim_params["β"])
                        η_gt = η_gt[1:data_point_len]
                        avg_ηList = readdlm(joinpath(win_exp_path,"data","avg_η.csv"), ',', Float64)
                        ratio_gt = η_gt ./ β_gt
                        product_gt = η_gt .* β_gt
                        sum_gt = η_gt .+ β_gt
                    end

                    gt_h = gt_h_[1:(data_point_len+1)]
                    est_h_list = est_h_list[1:data_point_len+1, :]
                    

                    pred_surface_pt_lst, _, _ = read_csv(joinpath(win_exp_path,"data","sim_data","3D_surface_points_pred"))
                    gt_surface_pt_lst, _, _ = read_csv(joinpath(filepath_gt,"data","sim_data","surface_nodes"))
                    est_surface_pt_lst, _, _ = read_csv(joinpath(win_exp_path,"data","sim_data","3D_surface_points_est"))

                    t_full_h = collect(range(start=0, stop=effective_sim_time, step=t_steps))

                    qoi_est = get_surface_qoi(est_surface_pt_lst, obj_pose, h)
                    qoi_pred = get_surface_qoi(pred_surface_pt_lst, obj_pose, h)
                    qoi_gt = get_surface_qoi(gt_surface_pt_lst, obj_pose, h)
                    plt_surface_qoi = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                            collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        else
                            collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                        end
                        t = t_windows[ti]
                        Plots.vline!(plt_surface_qoi, [t], color=:gray, linestyle=:dash, label=false)
                        if ti === 1
                            Plots.plot!(plt_surface_qoi, t_full_h[range_], qoi_pred[range_pred]./qoi_gt[range_], label=L"\mathrm{QoI}_{\mathrm{pred}}/\mathrm{QoI}_{\mathrm{gt}}", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_surface_qoi, t_full_h[range_], qoi_pred[range_pred]./qoi_gt[range_], label=false, color=def_blue, linestyle=:dash)
                        end
                        pred_range_start = range_pred[end] + 1
                    end
                    Plots.plot!(plt_surface_qoi, t_full_h, qoi_est./qoi_gt, label=L"\mathrm{QoI}_{\mathrm{est}}/\mathrm{QoI}_{\mathrm{gt}}", color=def_red)
                    Plots.xlabel!(plt_surface_qoi, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_surface_qoi, L"\mathrm{QoI}_{\mathrm{est}}/\mathrm{QoI}_{\mathrm{gt}}")
                    Plots.xlims!(plt_surface_qoi, 0, end_obs_win)
                    Plots.ylims!(plt_surface_qoi, y_lims_h_norm)
                    Plots.savefig(plt_surface_qoi, joinpath(win_exp_path,"plots","surface_area_qoi.pdf"))
                    
                    dc_surface_est, dh_surface_est, dcp_surface_est = compare_pt_clouds(est_surface_pt_lst, gt_surface_pt_lst)
                    plt_surface_error_dc = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    plt_surface_error_dh = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    plt_surface_error_dcp = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        dc_surface_pred, dh_surface_pred, dcp_surface_pred = compare_pt_clouds(pred_surface_pt_lst[range_pred], gt_surface_pt_lst[range_])
                        t = t_windows[ti]
                        Plots.vline!(plt_surface_error_dc, [t], color=:gray, linestyle=:dash, label=false)
                        Plots.vline!(plt_surface_error_dh, [t], color=:gray, linestyle=:dash, label=false)
                        Plots.vline!(plt_surface_error_dcp, [t], color=:gray, linestyle=:dash, label=false)
                        if ti === 1
                            Plots.plot!(plt_surface_error_dc, t_full_h[range_], dc_surface_pred, label="Prediction error", color=def_blue, linestyle=:dash)
                            Plots.plot!(plt_surface_error_dh, t_full_h[range_], dh_surface_pred, label="Prediction error", color=def_blue, linestyle=:dash)
                            Plots.plot!(plt_surface_error_dcp, t_full_h[range_], dcp_surface_pred, label="Prediction error", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_surface_error_dc, t_full_h[range_], dc_surface_pred, label=false, color=def_blue, linestyle=:dash)
                            Plots.plot!(plt_surface_error_dh, t_full_h[range_], dh_surface_pred, label=false, color=def_blue, linestyle=:dash)
                            Plots.plot!(plt_surface_error_dcp, t_full_h[range_], dcp_surface_pred, label=false, color=def_blue, linestyle=:dash)
                        end
                        pred_range_start = range_pred[end] + 1
                    end
                    Plots.plot!(plt_surface_error_dc, t_full_h, dc_surface_est, label="Estimation error", color=def_red)
                    Plots.xlabel!(plt_surface_error_dc, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_surface_error_dc, L"\mathrm{Hausdorff\;Distance\;[mm]}")
                    Plots.xlims!(plt_surface_error_dc, 0, end_obs_win)
                    Plots.savefig(plt_surface_error_dc, joinpath(win_exp_path,"plots","surface_point_distance_error_dh.pdf"))

                    Plots.plot!(plt_surface_error_dh, t_full_h, dh_surface_est, label="Estimation error", color=def_red)
                    Plots.xlabel!(plt_surface_error_dh, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_surface_error_dh, L"\mathrm{Chamfer\;Distance\;[mm]}")
                    Plots.xlims!(plt_surface_error_dh, 0, end_obs_win)
                    Plots.savefig(plt_surface_error_dh, joinpath(win_exp_path,"plots","surface_point_distance_error_dc.pdf"))

                    Plots.plot!(plt_surface_error_dcp, t_full_h, dcp_surface_est, label="Estimation error", color=def_red)
                    Plots.xlabel!(plt_surface_error_dcp, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_surface_error_dcp, L"\mathrm{Closest\;Point\;Distance\;[mm]}")
                    Plots.xlims!(plt_surface_error_dcp, 0, end_obs_win)
                    Plots.savefig(plt_surface_error_dcp, joinpath(win_exp_path,"plots","surface_point_distance_error_dcp.pdf"))

                    d_est, _ = closest_point(sim_border_pt_lst, obs_border_pt_lst)

                    plt_cnt_error = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        d_pred, _ = closest_point(pred_border_pt_lst[range_pred], obs_border_pt_lst[range_])
                        t = t_windows[ti]
                        Plots.vline!(plt_cnt_error, [t], color=:gray, linestyle=:dash, label=false)
                        if ti === 1
                            Plots.plot!(plt_cnt_error, t_full_h[range_], d_pred./cost_init, label="Prediction error", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_cnt_error, t_full_h[range_], d_pred./cost_init, label=false, color=def_blue, linestyle=:dash)
                        end
                        pred_range_start = range_pred[end] + 1
                    end
                    Plots.plot!(plt_cnt_error, t_full_h, d_est./cost_init, label="Estimation error", color=def_red)
                    Plots.xlabel!(plt_cnt_error, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_cnt_error, L"d(t)/d^{[0]}")
                    Plots.xlims!(plt_cnt_error, 0, end_obs_win)
                    Plots.ylims!(plt_cnt_error, 0, 4)
                    Plots.savefig(plt_cnt_error, joinpath(win_exp_path,"plots","closest_point_distance_error.pdf"))
                    
                    t_full = collect(range(start=t_steps, stop=effective_sim_time, step=t_steps))

                    plt_η = set_plot(fs, sz=(plt_width, plt_height), legend_column=4, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    Plots.plot!(plt_η, [], label=false, legend=:outerbottom, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        Plots.vline!(plt_η, [t], color=:gray, linestyle=:dash, label=false)
                    end
                    t_prev = 0.1
                    prev_η = 0.0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plt_η, t_win, est_ηpList[data_range_], label=L"\eta_{\mathrm{est}}(t)", color=def_red)
                            Plots.plot!(plt_η, [], label=L"\eta_{\mathrm{pred}}(t)", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_η, t_win, est_ηpList[data_range_], color=def_red, label=false)
                            Plots.plot!(plt_η, t_win, prev_η*ones(length(t_win)), label=false, color=def_blue, linestyle=:dash)
                        end
                        prev_η = est_ηpList[data_range_[end]]
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.plot!(plt_η, t_full, η_gt, label=L"\eta_{\mathrm{gt}}(t)", color=def_green)
                    end
                    Plots.xlabel!(plt_η, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_η, latexstring("\$\\eta(t)\$ [kPa s]"))
                    Plots.xlims!(plt_η, 0, end_obs_win)
                    Plots.savefig(plt_η, joinpath(win_exp_path,"plots","η.pdf"))
                    
                    plt_β = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    t_prev = 0.1
                    prev_β = 0.0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        Plots.vline!(plt_β, [t], color=:gray, linestyle=:dash, label=false)
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plt_β, t_win, est_βpList[data_range_], label=L"\beta_{\mathrm{est}}(t)", color=def_red)
                            Plots.plot!(plt_β, [], label=L"\beta_{\mathrm{pred}}(t)", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_β, t_win, est_βpList[data_range_], color=def_red, label=false)
                            Plots.plot!(plt_β, t_win, prev_β*ones(length(t_win)), label=false, color=def_blue, linestyle=:dash)
                        end
                        prev_β = est_βpList[data_range_[end]]
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.hline!(plt_β, β_gt, label=L"\beta_{\mathrm{gt}}", color=def_green)
                    end
                    Plots.xlabel!(plt_β, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plt_β, latexstring("\$\\beta(t)\$ [MPa s m\$^{-1}\$]"))
                    Plots.xlims!(plt_β, 0, end_obs_win)
                    Plots.savefig(plt_β, joinpath(win_exp_path,"plots","β.pdf"))
                    
                    plot_param_ratio = set_plot(fs, sz=(plt_width, plt_height), legend_column=2, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    t_prev = 0.1
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plot_param_ratio, t_win, ratio_opt[data_range_], label=L"\eta_{\mathrm{est}}/\beta_{\mathrm{est}}", color=def_red)
                        else
                            Plots.plot!(plot_param_ratio, t_win, ratio_opt[data_range_], color=def_red, label=false)
                        end
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.plot!(plot_param_ratio, t_full, ratio_gt, label=L"\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}", color=def_green)
                    end
                    Plots.xlabel!(plot_param_ratio, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plot_param_ratio, latexstring("\$\\eta/\\beta\$ [mm\$^{-1}\$]"))
                    Plots.xlims!(plot_param_ratio, 0, end_obs_win)
                    Plots.savefig(plot_param_ratio, joinpath(win_exp_path,"plots","η_over_β.pdf"))    

                    plot_param_product = set_plot(fs, sz=(plt_width, plt_height), legend_column=2, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    t_prev = 0.1
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plot_param_product, t_win, product_opt[data_range_], label=L"\eta_{\mathrm{est}} \cdot \beta_{\mathrm{est}}", color=def_red)
                        else
                            Plots.plot!(plot_param_product, t_win, product_opt[data_range_], color=def_red, label=false)
                        end
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.plot!(plot_param_product, t_full, product_gt, label=L"\eta_{\mathrm{gt}} \cdot \beta_{\mathrm{gt}}", color=def_green)
                    end
                    Plots.xlabel!(plot_param_product, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(plot_param_product, latexstring("\$\\eta \\cdot \\beta\$"))
                    Plots.xlims!(plot_param_product, 0, end_obs_win)
                    Plots.savefig(plot_param_product, joinpath(win_exp_path,"plots","η_times_β.pdf"))

                    plot_param_sum = set_plot(fs, sz=(plt_width, plt_height), legend_column=2, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    t_prev = 0.1
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plot_param_sum, t_win, sum_opt[data_range_], label=L"\eta_{\mathrm{est}} + \beta_{\mathrm{est}}", color=def_red)
                        else
                            Plots.plot!(plot_param_sum, t_win, sum_opt[data_range_], color=def_red, label=false)
                        end
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.plot!(plot_param_sum, t_full, sum_gt, label=L"\eta_{\mathrm{gt}} + \beta_{\mathrm{gt}}", color=def_green)
                    end
                    Plots.xlabel!(plot_param_sum, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(latexstring("\$\\eta + \\beta\$"))
                    Plots.xlims!(plot_param_sum, 0, end_obs_win)
                    Plots.savefig(plot_param_sum, joinpath(win_exp_path,"plots","η_plus_β.pdf"))
                    
                    h_plt = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    Plots.plot!(h_plt, t_full_h, est_h_list, label=L"h_{\mathrm{est}}", color=def_red)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            Plots.plot!(h_plt, [], label=L"h_{\mathrm{pred}}", linestyle=:dash, color=def_blue)
                        else
                            Plots.plot!(h_plt, t_full_h[range_], pred_h_list[range_pred], linestyle=:dash, color=def_blue, label=false)
                        end
                        Plots.vline!(h_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1
                    end
                    if data_type != "physical"
                        Plots.plot!(h_plt, t_full_h, gt_h, label=L"h_{\mathrm{gt}}", color=def_green)
                    else
                        Plots.plot!(h_plt, t_full_h, gt_h, label=L"h_{\mathrm{m}}", color=def_green)
                    end
                    Plots.xlabel!(h_plt, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(h_plt, L"\mathrm{Height\;[mm]}")
                    Plots.xlims!(h_plt, 0, end_obs_win)
                    Plots.ylims!(h_plt, minimum(vcat(gt_h, est_h_list))*0.8, maximum(vcat(gt_h, est_h_list))*1.2)
                    Plots.savefig(joinpath(win_exp_path,"plots","h.pdf"))

                    h_normalized_plt = set_plot(fs, sz=(plt_width, plt_height), legend_column=2, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            if data_type != "physical"
                                Plots.plot!(h_normalized_plt, t_full_h[range_], est_h_list[range_]./gt_h[range_], label=L"h_{\mathrm{est}}/h_{\mathrm{gt}}", color=def_red)
                                Plots.plot!(h_normalized_plt, [], label=L"h_{\mathrm{pred}}/h_{\mathrm{gt}}", linestyle=:dash, color=def_blue)
                            else
                                Plots.plot!(h_normalized_plt, t_full_h[range_], est_h_list[range_]./gt_h[range_], label=L"h_{\mathrm{est}}/h_{\mathrm{m}}", color=def_red)
                                Plots.plot!(h_normalized_plt, [], label=L"h_{\mathrm{pred}}/h_{\mathrm{m}}", linestyle=:dash, color=def_blue)
                            end
                        else
                            Plots.plot!(h_normalized_plt, t_full_h[range_], pred_h_list[range_pred]./gt_h[range_], linestyle=:dash, color=def_blue, label=false)
                            Plots.plot!(h_normalized_plt, t_full_h[range_], est_h_list[range_]./gt_h[range_], color=def_red, label=false)
                        end
                        Plots.vline!(h_normalized_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1 
                    end
                    Plots.xlabel!(h_normalized_plt, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(h_normalized_plt, L"h_{\mathrm{est}}/h_{\mathrm{gt}}")
                    Plots.xlims!(h_normalized_plt, 0, end_obs_win)
                    Plots.ylims!(h_normalized_plt, y_lims_h_norm)
                    Plots.savefig(joinpath(win_exp_path,"plots","h_normalized.pdf"))

                    error_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    Plots.plot!(error_plt, t_full_h, abs.(est_h_list-gt_h), label="Estimation", color=def_red)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            Plots.plot!(error_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                        else
                            Plots.plot!(error_plt, t_full_h[range_], abs.(pred_h_list[range_pred]-gt_h[range_]), linestyle=:dash, color=def_blue, label=false)
                        end
                        Plots.vline!(error_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1 
                    end
                    Plots.xlabel!(error_plt, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(error_plt, L"\mathrm{Height\;Error\;[mm]}")
                    Plots.xlims!(error_plt, 0, end_obs_win)
                    Plots.savefig(joinpath(win_exp_path,"plots","h_est_error.pdf"))

                    rel_error_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin, legend_column=2)
                    Plots.plot!(rel_error_plt, t_full_h, abs.(est_h_list-gt_h)./gt_h*100, label="Estimation", color=def_red)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            Plots.plot!(rel_error_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                        else
                            Plots.plot!(rel_error_plt, t_full_h[range_], abs.(pred_h_list[range_pred]-gt_h[range_])./gt_h[range_]*100, linestyle=:dash, color=def_blue, label=false)
                        end
                        Plots.vline!(rel_error_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1 
                    end
                    Plots.xlabel!(rel_error_plt, L"\mathrm{Time\;[s]}")
                    Plots.ylabel!(rel_error_plt, latexstring("Relative Error [\$\\%\$]"))
                    Plots.xlims!(rel_error_plt, 0, end_obs_win)
                    Plots.ylims!(rel_error_plt, y_lims_rel_error)
                    Plots.savefig(joinpath(win_exp_path,"plots","rel_error.pdf"))
                    
                end
            end
        end
    end
    return
end

"""
    post_analysis_const(filepath_gt_, filepath, avoid_list)

Cross-experiment post-analysis and figure generation for `viscosity_type ==
"constant"` results (the counterpart of `replot`, which handles a
single experiment). Walks the nested `dir/elem_size_folder/simtime_folder/
noise_folder` result tree under `filepath` (top-level entries in `avoid_list`,
or `"post_analysis_global"`, are skipped), reads each leaf experiment's
already-computed `optimize`/`replot` outputs (η/β estimates, cost history,
height, contours), and overlays them onto shared comparison figures at four
nested levels:

- **global** (across all top-level `dir`s): cost convergence, normalized
  η/β/ratio, height/normalized-height/relative-error, and contour comparison,
  bucketed by simulation-time case (2/5/10/20/30 s).
- **per-`dir`** (across element sizes): mesh-convergence plots of cost, η, β,
  their ratio, and height/relative-height-error.
- **per-`elem_size_folder`** (across simulation-time windows): the same set
  of comparison plots, plus cost-surface slice-direction plots.
- **per-`sim_time_folder`** (across noise levels): covariance ellipse and
  noisy height/error/ratio distribution plots.

# Arguments
- `filepath_gt_::String`: ground-truth data root each experiment is compared against.
- `filepath::String`: root of the result tree to walk.
- `avoid_list`: top-level `dir` names under `filepath` to skip.

# Returns
None. All outputs are PDFs written under `plot_path_global`/`plot_path_elems`/
`plot_path_sim_time`/`plot_path_noise`.
"""
function post_analysis_const(filepath_gt_::String, filepath::String, avoid_list)
    dir_list = readdir(filepath)

    # figure for legend

    # plot for convergence per slip case
    plot_conv_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_2, [],  label=false)
    Plots.xlabel!(plot_conv_2, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_2, L"\mathrm{Cost\;[px]}")
    
    plot_conv_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_5, [],  label=false)
    Plots.xlabel!(plot_conv_5, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_5, L"\mathrm{Cost\;[px]}")
    
    plot_conv_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_10, [],  label=false)
    Plots.xlabel!(plot_conv_10, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_10, L"\mathrm{Cost\;[px]}")

    plot_conv_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_20, [],  label=false)
    Plots.xlabel!(plot_conv_20, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_20, L"\mathrm{Cost\;[px]}")

    plot_conv_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_30, [],  label=false)
    Plots.xlabel!(plot_conv_30, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_30, L"\mathrm{Cost\;[px]}")

    plot_conv_log_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_log_2, [],  label=false)
    Plots.xlabel!(plot_conv_log_2, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_2, L"\mathrm{Cost\;[px]}")

    plot_conv_log_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_log_5, [],  label=false)
    Plots.xlabel!(plot_conv_log_5, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_5, L"\mathrm{Cost\;[px]}")

    plot_conv_log_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_log_10, [],  label=false)
    Plots.xlabel!(plot_conv_log_10, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_10, L"\mathrm{Cost\;[px]}")

    plot_conv_log_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_log_20, [],  label=false)
    Plots.xlabel!(plot_conv_log_20, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_20, L"\mathrm{Cost\;[px]}")

    plot_conv_log_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot!(plot_conv_log_30, [],  label=false)
    Plots.xlabel!(plot_conv_log_30, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_30, L"\mathrm{Cost\;[px]}")

    η_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(η_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_2,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_2,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(β_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_2, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_2, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(η_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_5,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(β_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black, yscale=:log10)
    Plots.xlabel!(β_norm_plot_5, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_5, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")
    
    η_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(η_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_10,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(β_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_10, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_10, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(η_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_20,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_20,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(β_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_20, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_20, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(η_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_30,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_30,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(β_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_30, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_30, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    ratio_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_plot_2, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_2,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_2,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    ratio_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_plot_5, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_5,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")
    
    ratio_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_plot_10, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_10,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    ratio_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_plot_20, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_20, L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_20, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    ratio_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_plot_30, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_30,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_30,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    h_glob_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.plot([5,end_obs_win],[0.998,0.998], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(h_glob_plot_2, 15, 0.997, ("Prediction",:black, :center,10,"computer modern"))
    Plots.xlabel!(h_glob_plot_2,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_glob_plot_2,L"h_{\mathrm{est}}\;\mathrm{[mm]}")
    Plots.xlims!(h_glob_plot_2, 0, end_obs_win)

    h_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.xlabel!(h_glob_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_glob_plot_5,L"h_{\mathrm{est}}\;\mathrm{[mm]}")  
    Plots.xlims!(h_glob_plot_5, 0, end_obs_win) 

    h_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.xlabel!(h_glob_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_glob_plot_10,L"h_{\mathrm{est}}\;\mathrm{[mm]}") 
    Plots.xlims!(h_glob_plot_10, 0, end_obs_win)


    h_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(h_norm_plot_2, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.plot!(h_norm_plot_2,[5,end_obs_win],[0.998,0.998], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(h_norm_plot_2, 15, 0.997, ("Prediction",:black, :center,10,"computer modern"))
    Plots.xlabel!(h_norm_plot_2,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_2,L"h_{\mathrm{est}}/h_{\mathrm{gt}}")
    Plots.xlims!(h_norm_plot_2, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_2, y_lims_h_norm)

    h_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(h_norm_plot_5, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)

    Plots.plot!(h_norm_plot_5,[5,end_obs_win],[1.01,1.01], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(h_norm_plot_5, 15, 1.0125, ("Prediction",:black, :center, 8,"computer modern"))
    Plots.vline!(h_norm_plot_5, [5.0], color=:black, linestyle=:dash, label=false)
    Plots.xlabel!(h_norm_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_5,L"h_{\mathrm{est}}/h_{\mathrm{gt}}")  
    Plots.xlims!(h_norm_plot_5, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_5, y_lims_h_norm)

    h_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(h_norm_plot_10, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.vline!(h_norm_plot_10, [10.0], color=:black, linestyle=:dash, label=false)
    Plots.xlabel!(h_norm_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_10,L"h_{\mathrm{est}}/h_{\mathrm{gt}}") 
    Plots.xlims!(h_norm_plot_10, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_10, y_lims_h_norm)

    rel_height_error_glob_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_2, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_2, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_2, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_2, y_lims_rel_error)

    rel_height_error_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.vline!(rel_height_error_glob_plot_5, [5.0], color=:black, linestyle=:dash, label=false)

    Plots.plot!(rel_height_error_glob_plot_5,[5,end_obs_win],[1.01,1.01], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(rel_height_error_glob_plot_5, 15, 1.0175, ("Prediction",:black, :center, 8,"computer modern"))
    Plots.xlabel!(rel_height_error_glob_plot_5, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_5, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_5, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_5, y_lims_rel_error)

    rel_height_error_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.vline!(rel_height_error_glob_plot_10, [10.0], color=:black, linestyle=:dash, label=false)
    Plots.xlabel!(rel_height_error_glob_plot_10, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_10, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_10, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_10, y_lims_rel_error)

    rel_height_error_glob_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.vline!(rel_height_error_glob_plot_20, [20.0], color=:black, linestyle=:dash, label=false)
    Plots.xlabel!(rel_height_error_glob_plot_20, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_20, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_20, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_20, y_lims_rel_error)

    rel_height_error_glob_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_30, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_30, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_30, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_30, y_lims_rel_error)

    ratio_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_2,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_2,L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

    ratio_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_5,L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

    ratio_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_10,L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

    ratio_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_20, L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_20, L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

    ratio_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_30,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_30,L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")
    
    cont_y_lims = [400, 1200] 
    cont_x_lims = [1380, 1430]
    contour_plt = set_plot(fs, sz=(round(Int,plt_width*1.5), plt_height), legend_column=3)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt, 480, 1520)
    Plots.xlabel!(contour_plt, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt, L"y\;\mathrm{[px]}")

    cnt_plt_width = 179
    contour_plt_zoom = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom, cont_x_lims[1], cont_x_lims[2])
    Plots.xlabel!(contour_plt_zoom, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt_zoom, L"y\;\mathrm{[px]}")

    contour_plt_zoom_05 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_05, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_05, cont_x_lims[1], cont_x_lims[2])
    Plots.xlabel!(contour_plt_zoom_05, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt_zoom_05, L"y\;\mathrm{[px]}")

    contour_plt_zoom_10 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_10, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_10, cont_x_lims[1], cont_x_lims[2])
    Plots.xlabel!(contour_plt_zoom_10, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt_zoom_10, L"y\;\mathrm{[px]}")

    contour_plt_zoom_15 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_15, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_15, cont_x_lims[1], cont_x_lims[2])
    Plots.xlabel!(contour_plt_zoom_15, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt_zoom_15, L"y\;\mathrm{[px]}")

    contour_plt_zoom_20 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_20, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_20, cont_x_lims[1], cont_x_lims[2])
    Plots.xlabel!(contour_plt_zoom_20, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt_zoom_20, L"y\;\mathrm{[px]}")

    cont_plt_legend = set_plot(fs, sz=(cnt_plt_width, plt_height), legend=:outertopright)

    # Create a thin horizontal legend strip for the outer loop (directory/viscosity slip case legend)
    slip_case_legend = set_plot(11, sz=(round(Int, plt_width*1.7), 50), legend=:bottom, legend_column=3, bottom_margin=-35mm, top_margin=2mm, left_margin=-25mm, right_margin=-15mm)
    Plots.plot!(slip_case_legend, [0, 1], [0, 0], label=false, color=:white, linewidth=0)  # Hidden base plot for axes
    Plots.xlims!(slip_case_legend, -0.2, 1.2)
    Plots.ylims!(slip_case_legend, -0.5, 0.5)
    
    global color_palette

    max_iter = 1
    for dir in dir_list
    
        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(joinpath(filepath_gt,"data","sim_params"))
        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        gt_h = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        println("Ground truth η: ", η_gt[1])

        Plots.plot!(slip_case_legend, [0, 0.2], [0, 0], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], linewidth=2, markerstrokewidth=0)

        elem_size_folders = readdir(filepath_dir)

        elem_conv_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.xlabel!(elem_conv_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_conv_plt, L"\mathrm{Cost\;[px]}")

        elem_conv_plt_log = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.xlabel!(elem_conv_plt_log, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_conv_plt_log, L"\mathrm{Cost\;[px]}")

        elem_η_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.hline!(elem_η_plt, [η_gt[1]], label="Ground truth η",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_η_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_plt,L"\eta\;\mathrm{[kPa\, s]}")

        elem_β_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.hline!(elem_β_plt, [β_gt[1]], label="Ground truth β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_β_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_plt, L"\beta\;\mathrm{[Pa\, s \, m]}")

        elem_ratio_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.hline!(elem_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_ratio_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_plt, L"\eta/\beta\;\mathrm{[mm^{-1}]}")
        
        elem_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.hline!(elem_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        Plots.xlabel!(elem_η_norm_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_norm_plt,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

        elem_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.hline!(elem_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        Plots.xlabel!(elem_β_norm_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_norm_plt, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

        elem_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.hline!(elem_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        Plots.xlabel!(elem_ratio_norm_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_norm_plt, L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

        elem_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.plot!(elem_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_rel_height_error_plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(elem_rel_height_error_plt, latexstring("Relative Height Error [\$\\%\$]"))
        Plots.xlims!(elem_rel_height_error_plt, 0, end_obs_win)

        elem_height_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
        Plots.plot!(elem_height_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_height_plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(elem_height_plt, L"\mathrm{Height\;[mm]}")
        Plots.xlims!(elem_height_plt, 0, end_obs_win)
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_16"
                continue
            end
            
            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)
            
            sim_window_η_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_plt, [η_gt[1]], label="Ground truth η")
            Plots.xlabel!(sim_window_η_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_plt,L"\eta\;\mathrm{[kPa\, s]}")
            
            sim_window_β_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_plt, [β_gt[1]], label="Ground truth β")
            Plots.xlabel!(sim_window_β_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_plt, L"\beta\;\mathrm{[Pa\, s \, m]}")

            sim_window_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
            Plots.xlabel!(sim_window_ratio_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_plt, L"\eta/\beta\;\mathrm{[mm^{-1}]}")

            sim_window_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            Plots.xlabel!(sim_window_η_norm_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_norm_plt,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

            sim_window_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            Plots.xlabel!(sim_window_β_norm_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_norm_plt, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

            sim_window_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            Plots.xlabel!(sim_window_ratio_norm_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_norm_plt, L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

            sim_window_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_rel_height_error_plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(sim_window_rel_height_error_plt, latexstring("Relative Height Error [\$\\%\$]"))
            Plots.xlims!(sim_window_rel_height_error_plt, 0, end_obs_win)

            sim_window_height_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_height_plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(sim_window_height_plt, L"\mathrm{Height\;[mm]}")
            Plots.xlims!(sim_window_height_plt, 0, end_obs_win)

            plt_slices = set_plot(fs, sz=(350, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum",  lw=1)
            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;[px]}")
            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
            Plots.ylims!(plt_slices, 0, 50)
                
            for sim_time_folder_ in reverse(sort(sim_time_folders))

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" || sim_time_folder_ == "simtime_20.0" || sim_time_folder_ == "simtime_10.0" || sim_time_folder_ == "simtime_30.0"
                    continue
                end
                
                height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(height_error_plt, [], label=false)
                Plots.xlabel!(height_error_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(height_error_plt, L"\mathrm{Height\;Error\;[mm]}")
                Plots.xlims!(height_error_plt, 0, end_obs_win)
                
                rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(rel_height_error_plt, [], label=false)
                Plots.xlabel!(rel_height_error_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(rel_height_error_plt, latexstring("Relative Height Error [\$\\%\$]"))
                Plots.xlims!(rel_height_error_plt, 0, end_obs_win)
                
                η_β_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(η_β_ratio_plt, [], label=false)
                Plots.hline!(η_β_ratio_plt, [1], label=false, linestyle=:dash)
                Plots.xlabel!(η_β_ratio_plt, L"\mathrm{Iterations}")
                Plots.ylabel!(η_β_ratio_plt, L"\eta/\beta\;\mathrm{[mm^{-1}]}")
                
                noise_cols::Int = 2
                covarience_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.xlabel!(covarience_plt, L"\eta/\eta_{\mathrm{gt}}")
                Plots.ylabel!(covarience_plt, L"\beta/\beta_{\mathrm{gt}}")
                Plots.xlims!(covarience_plt, 0.95, 1.05)
                Plots.ylims!(covarience_plt, 0.0, 1.75)
                
                dist_plot = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.vline!(dist_plot, [1.0], linestyle=:dash, color=:black, label=false)
                Plots.xlabel!(dist_plot, L"\frac{\eta/\beta}{(\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")
                Plots.ylabel!(dist_plot, L"\mathrm{Density}")

                height_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.xlabel!(height_noise_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(height_noise_plt, L"\mathrm{Height\;[mm]}")
                Plots.xlims!(height_noise_plt, 0, end_obs_win)

                height_error_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.vline!(height_error_noise_plt, [5.0], color=:black, linestyle=:dash, label=false)
                Plots.xlabel!(height_error_noise_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(height_error_noise_plt, L"\mathrm{Height\;Error\;[mm]}")
                Plots.xlims!(height_error_noise_plt, 0, end_obs_win)

                rel_height_error_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.vline!(rel_height_error_noise_plt, [5.0], color=:black, linestyle=:dash, label=false)
                Plots.xlabel!(rel_height_error_noise_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(rel_height_error_noise_plt, latexstring("Relative Height Error [\$\\%\$]"))
                Plots.xlims!(rel_height_error_noise_plt, 0, end_obs_win)

                normalized_height_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.vline!(normalized_height_noise_plt, [5.0], color=:black, linestyle=:dash, label=false)
                Plots.hline!(normalized_height_noise_plt, [1.0], linestyle=:dash, color=:black, label=false)
                Plots.ylabel!(normalized_height_noise_plt, L"h_{\mathrm{est}}/h_{\mathrm{gt}}")
                Plots.xlabel!(normalized_height_noise_plt, L"\mathrm{Time\;[s]}")
                Plots.xlims!(normalized_height_noise_plt, 0, end_obs_win)

                ratio_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.hline!(ratio_noise_plt, [1.0], linestyle=:dash, color=:black, label=false)
                Plots.xlabel!(ratio_noise_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(ratio_noise_plt, L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

                η_noise_norm_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.hline!(η_noise_norm_plt, [1.0], linestyle=:dash, color=:black, label=false)
                Plots.xlabel!(η_noise_norm_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(η_noise_norm_plt, L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")
                
                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_)
                noise_folders = readdir(sim_time_folder)

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                noise_iter = 0
                min_η_norm, max_η_norm, min_β_norm, max_β_norm = 0.0, 0.0, 0.0, 0.0
                for noise_folder_ in reverse(sort(noise_folders))
                    
                    if noise_folder_ == "post_analysis_noise" || noise_folder_ == "Results" 
                        continue
                    end
                    
                    exp_path = joinpath(filepath_dir, elem_size_folder_, sim_time_folder_, noise_folder_)
                    printstyled("Processing experiment folder: $(exp_path)\n", color=:magenta)
                    exp_params = read_json(joinpath(exp_path,"data","experiment_parameters"))

                    noise_level = exp_params["noise_level"]
                    ne = exp_params["ne_exp"]
                    sim_time = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    exp_path_n0 = replace(exp_path, "$noise_folder_" => "noise_0.0")
                    num_exp_points::Int = round(Int,sim_time/t_steps)

                    printstyled("Processing for noise level: $(noise_level): $(exp_path)\n", color=:yellow)
                    obs_border_pt_lst, sim_border_pt_lst, gt_Splinex, gt_Spliney, splinex, spliney, _ = _get_borders(data_type, filepath_gt, exp_path_n0, num_exp_points)
                    if noise_level == 0.0
                        
                        est_η = readdlm(joinpath(exp_path,"data","η.csv"), ',', Float64)
                        est_β = readdlm(joinpath(exp_path,"data","β.csv"), ',', Float64)
                        est_h = readdlm(joinpath(exp_path,"data","est_h.csv"), ',', Float64)
                        iter = collect(range(start=1, stop=size(est_η,1), step=1))

                        max_iter = max(max_iter, size(iter,1))

                        if length(est_h) != length(gt_h)
                            @warn "Length mismatch between estimated height ($(length(est_h))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(length(est_h), length(gt_h))
                            est_h = est_h[1:min_length]
                            gt_h = gt_h[1:min_length]
                        end

                        height_error = abs.(est_h - gt_h)
                        rel_height_error = height_error ./ gt_h * 100.0

                        # normalize by ground-truth                        
                        est_η_norm = est_η./η_gt
                        est_β_norm = est_β./β_gt
                        h_norm = est_h ./ gt_h

                        ratio_est = est_η ./ est_β
                        ratio_gt = η_gt ./ β_gt
                        normalized_ratio = (est_η ./ est_β) * (β_gt ./ η_gt)

                        cost_list = readdlm(joinpath(exp_path,"data","cost_iter.csv"), ',', Float64)
                
                        Plots.plot!(sim_window_η_plt, iter, est_η, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_β_plt, iter, est_β, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_ratio_plt, iter, ratio_est, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        
                        time_h = copy(time)
                        n_time = min(length(time), length(rel_height_error))
                        if length(rel_height_error) != length(time) 
                            @warn "Length mismatch between height plot x-axis ($(length(rel_height_error))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(length(rel_height_error), length(time))
                            time_h = time[1:min_length]
                            rel_height_error = rel_height_error[1:min_length]
                            est_h = est_h[1:min_length]
                            gt_h = gt_h[1:min_length]
                        end
                        Plots.plot!(sim_window_rel_height_error_plt, time_h, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_height_plt, time_h, est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_height_plt, time_h, gt_h, label=L"h_{\mathrm{gt}}(t)", legend=:outerbottom, legend_column=2)
                        
                        Plots.plot!(sim_window_η_norm_plt, iter, est_η_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_β_norm_plt, iter, est_β_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_ratio_norm_plt, iter, normalized_ratio, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        
                        try
                            slice_data = read_json(joinpath(exp_path,"data","slice_data"))
                            t_steep = Float64.(collect(slice_data["steep"]["t"]))
                            zs_steep = Float64.(collect(slice_data["steep"]["zs"]))
                            t_flat = Float64.(collect(slice_data["flat"]["t"]))
                            zs_flat = Float64.(collect(slice_data["flat"]["zs"]))
                            
                            # 2D slices: cost vs distance along the two directions
                            if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                Plots.plot!(plt_slices, t_steep, zs_steep, legend = :outerright, label = "Steepest direction; Window = $(sim_time)s")
                            else
                                @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                            end

                            if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                Plots.plot!(plt_slices, t_flat,  zs_flat, legend = :outerright, label = "Flattest direction; Window = $(sim_time)s", legendfontsize=20)
                            else
                                @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                            end
                        catch e
                            @warn "Failed to read or plot slice data: $e"
                        end

                        if sim_time == 2.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_2, iter, est_η_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(η_norm_plot_2, 0:2:(max_iter+1))

                                Plots.plot!(β_norm_plot_2, iter, est_β_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(β_norm_plot_2, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_2, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(plot_conv_2, 0:2:(max_iter+1))
                                
                                Plots.plot!(plot_conv_log_2, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], yscale=:log10)
                                Plots.xticks!(plot_conv_log_2, 1:10:max_iter)
                                
                                Plots.plot!(ratio_plot_2, iter, ratio_est, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], yscale=:log10)
                                Plots.xticks!(ratio_plot_2, 0:2:(max_iter+1))
                                Plots.hline!(ratio_plot_2, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)

                                Plots.plot!(ratio_norm_plot_2, iter, normalized_ratio, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(ratio_norm_plot_2, 0:2:(max_iter+1))
                                
                                Plots.plot!(rel_height_error_glob_plot_2, time_h, rel_height_error, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.plot!(h_norm_plot_2, h_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            end
                        elseif sim_time == 5.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_5, iter, est_η_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(η_norm_plot_5, 0:2:(max_iter+1))
                                
                                Plots.plot!(β_norm_plot_5, iter, est_β_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(β_norm_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_5, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(plot_conv_5, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_log_5, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], yscale=:log10)
                                Plots.xticks!(plot_conv_log_5, 1:2:max_iter)

                                Plots.plot!(ratio_plot_5, iter, ratio_est, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], yscale=:log10)
                                Plots.hline!(ratio_plot_5, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                                Plots.xticks!(ratio_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(ratio_norm_plot_5, iter, normalized_ratio, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.xticks!(ratio_norm_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(rel_height_error_glob_plot_5, time_h[1:51], rel_height_error[1:51], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.plot!(rel_height_error_glob_plot_5, time_h, rel_height_error, label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], linestyle=:dash)
                                
                                Plots.plot!(h_norm_plot_5, time_h[1:51], h_norm[1:51], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.plot!(h_norm_plot_5, time_h, h_norm[1:n_time], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], linestyle=:dash)
                                
                                Plots.plot!(h_glob_plot_5, time_h[1:51], h_norm[1:51], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.plot!(h_glob_plot_5, time_h, h_norm[1:n_time], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], linestyle=:dash)
                                Plots.plot!(h_glob_plot_5, time_h, gt_h[1:n_time], label=false, style=:dash, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                
                                Plots.plot!(contour_plt, gt_Splinex[end], gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                Plots.plot!(contour_plt_zoom, gt_Splinex[end], gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            end
                            println(size(height_error, 1), " vs ", size(gt_h), " vs ", size(time_h), " vs ", length(rel_height_error))
                            Plots.plot!(height_error_plt, time_h, height_error[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)

                            Plots.plot!(rel_height_error_plt, time_h, rel_height_error[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)

                            Plots.plot!(elem_conv_plt, iter, cost_list, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_conv_plt, 0:2:(max_iter+1))
                            
                            Plots.plot!(elem_conv_plt_log, iter, cost_list, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom, yscale=:log10)
                            Plots.xticks!(elem_conv_plt_log, 0:2:(max_iter+1))

                            Plots.plot!(elem_η_plt, iter, est_η, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_η_plt, 0:2:(max_iter+1))
                            
                            Plots.plot!(elem_β_plt, iter, est_β, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_β_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_ratio_plt, iter, ratio_est, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_ratio_plt, 0:2:(max_iter+1))
                            
                            Plots.plot!(elem_rel_height_error_plt, time_h, rel_height_error[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), legend=:outerbottom)
                            Plots.plot!(elem_height_plt, time_h, est_h[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), legend=:outerbottom)
                            Plots.plot!(elem_height_plt, time_h, gt_h[1:n_time], label=L"h_{\mathrm{gt}}(t)", legend=:outerbottom)
                        
                            Plots.plot!(elem_η_norm_plt, iter, est_η_norm, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_η_norm_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_β_norm_plt, iter, est_β_norm, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_β_norm_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_ratio_norm_plt, iter, normalized_ratio[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_ratio_norm_plt, 0:2:(max_iter+1))
                        end
                    else
                        η_pred = readdlm(joinpath(exp_path,"data","eta_est.csv"), ',', Float64) # estimated η values per sample
                        β_pred = readdlm(joinpath(exp_path,"data","beta_est.csv"), ',', Float64) # estimated β values per sample
                        h_pred = readdlm(joinpath(exp_path,"data","h_est.csv"), ',', Float64) # estimated height values per sample

                        n_obs_border_pt_lst, n_gt_Splinex, n_gt_Spliney = add_noise(obs_border_pt_lst, nFactor=noise_level)

                        η_β_pred = η_pred[:,1] ./ β_pred[:,1].* (β_gt ./ η_gt) # normalize the ratio by ground truth to compare across different noise levels
                        dist = fit(Normal, η_β_pred)
                        StatsPlots.plot!(dist_plot, dist, label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), legend_column=noise_cols)

                        η_norm = η_pred[:,1] ./ η_gt
                        β_norm = β_pred[:,1] ./ β_gt

                        if minimum(η_norm) < min_η_norm || min_η_norm == 0.0
                            min_η_norm = minimum(η_norm)
                        end
                        if maximum(η_norm) > max_η_norm || max_η_norm == 0.0
                            max_η_norm = maximum(η_norm)
                        end
                        if minimum(β_norm) < min_β_norm || min_β_norm == 0.0
                            min_β_norm = minimum(β_norm)
                        end
                        if maximum(β_norm) > max_β_norm || max_β_norm == 0.0
                            max_β_norm = maximum(β_norm)
                        end

                        n_samples = size(η_pred, 2)
                        η_β_list = []
                        η_β_norm_list = []
                        iter_list = []
                        for n in 1:n_samples
                            β = readdlm(joinpath(exp_path,"data","opt_data","beta_steps","run_$n.csv"), ',', Float64)
                            η = readdlm(joinpath(exp_path,"data","opt_data","eta_steps","run_$n.csv"), ',', Float64)
                            iter = readdlm(joinpath(exp_path,"data","opt_data","iter","run_$n.csv"), ',', Float64)
                            ratio = η./β 
                            normalized_ratio = ratio * (β_gt ./ η_gt)
                            push!(η_β_list, ratio[end])
                            push!(η_β_norm_list, normalized_ratio[end])
                            push!(iter_list, iter)
                        end
                        plot_covariance!(covarience_plt, η_pred[:,1]./η_gt, β_pred[:,1]./β_gt, label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), legend_column=noise_cols, color_ellipse=color_palette[(noise_iter) % length(color_palette)+1])

                        if noise_level == 0.5
                            Plots.plot!(contour_plt_zoom_05, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(cont_plt_legend, [], [], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                        elseif noise_level == 1.0
                            Plots.plot!(contour_plt_zoom_10, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                        elseif noise_level == 1.5
                            Plots.plot!(contour_plt_zoom_15, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                        elseif noise_level == 2.0
                            Plots.plot!(contour_plt_zoom_20, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                        end

                        n_time = min(length(time), size(h_pred, 2), length(gt_h))
                        if n_time < length(time) || n_time < size(h_pred, 2) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(size(h_pred, 2)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end
                        height_error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time])
                        normalized_height = h_pred[:,1:n_time]' ./ gt_h[1:n_time]
                        rel_height_error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time]) ./ gt_h[1:n_time] * 100.0 # in percentage

                        time_h = copy(time)
                        if size(rel_height_error,1) != length(time) 
                            @warn "Length mismatch between height plot x-axis ($(size(rel_height_error,1))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(size(rel_height_error,1), length(time))
                            time_h = time[1:min_length] 
                            h_pred = h_pred[:,1:min_length]
                            gt_h = gt_h[1:min_length]
                        end
                        Plots.plot!(height_noise_plt, time_h[1:n_time], gt_h[1:n_time], style=:dash, label=false)
                        StatsPlots.errorline!(height_noise_plt, time_h[1:n_time], h_pred[:,1:n_time]', label=false, groupcolor=color_palette[(noise_iter) % length(color_palette)+1], linestyle=:dash)
                        StatsPlots.errorline!(height_noise_plt, time_h[1:51], h_pred'[1:51,:], label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), groupcolor=color_palette[(noise_iter) % length(color_palette)+1])
                     
                        StatsPlots.errorline!(height_error_noise_plt, time_h, height_error', label=false, groupcolor=color_palette[(noise_iter) % length(color_palette)+1], linestyle=:dash)
                        StatsPlots.errorline!(height_error_noise_plt, time_h[1:51], height_error[1:51,:]', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), color=color_palette[(noise_iter) % length(color_palette)+1])
                        
                        StatsPlots.errorline!(rel_height_error_noise_plt, time_h, rel_height_error', label=false, groupcolor=color_palette[(noise_iter) % length(color_palette)+1], linestyle=:dash)
                        StatsPlots.errorline!(rel_height_error_noise_plt, time_h[1:51], rel_height_error[1:51,:]', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), groupcolor=color_palette[(noise_iter) % length(color_palette)+1])
                        
                        StatsPlots.errorline!(normalized_height_noise_plt, time_h, normalized_height', label=false, groupcolor=color_palette[(noise_iter) % length(color_palette)+1], linestyle=:dash)
                        StatsPlots.errorline!(normalized_height_noise_plt, time_h[1:51], normalized_height[1:51,:]', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), groupcolor=color_palette[(noise_iter) % length(color_palette)+1])
                        
                        noise_iter += 1
                    end
                end
                xlims = [min_η_norm, max_η_norm]
                ylims = [min_β_norm, max_β_norm]
                Plots.hline!(covarience_plt, [1], linestyle=:dash, label=false, color=:black)
                Plots.vline!(covarience_plt, [1], linestyle=:dash, label=false, color=:black)
                Plots.plot!(covarience_plt, bottom_margin=-20mm)

                plot_path_noise = joinpath(sim_time_folder,"post_analysis_noise","plots")
                set_file(plot_path_noise)
                Plots.savefig(covarience_plt, joinpath(plot_path_noise,"covariance_$dir.pdf"))
                Plots.savefig(height_noise_plt, joinpath(plot_path_noise,"height_$dir.pdf"))
                Plots.savefig(height_error_noise_plt, joinpath(plot_path_noise,"height_error_$dir.pdf"))
                Plots.savefig(rel_height_error_noise_plt, joinpath(plot_path_noise,"relative_height_error_$dir.pdf"))
                Plots.savefig(normalized_height_noise_plt, joinpath(plot_path_noise,"normalized_height_$dir.pdf"))
                Plots.savefig(ratio_noise_plt, joinpath(plot_path_noise,"eta_beta_ratio_$dir.pdf"))
                Plots.savefig(η_noise_norm_plt, joinpath(plot_path_noise,"eta_normalized_$dir.pdf"))
                Plots.savefig(dist_plot, joinpath(plot_path_noise,"eta_beta_ratio_distribution_$dir.pdf"))
                
                @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = joinpath(elem_size_folder,"post_analysis_time","plots")
            set_file(plot_path_sim_time)

            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η_$dir.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β_$dir.pdf"))
            Plots.savefig(sim_window_ratio_plt, joinpath(plot_path_sim_time,"η_β_ratio_$dir.pdf"))
            Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error_$dir.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height_$dir.pdf"))
            Plots.savefig(sim_window_η_norm_plt, joinpath(plot_path_sim_time,"η_normalized_$dir.pdf"))
            Plots.savefig(sim_window_β_norm_plt, joinpath(plot_path_sim_time,"β_normalized_$dir.pdf"))
            Plots.savefig(sim_window_ratio_norm_plt, joinpath(plot_path_sim_time,"η_β_ratio_normalized_$dir.pdf"))
            Plots.savefig(plt_slices, joinpath(plot_path_sim_time,"cost_slices_along_directions_$dir.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
        
        plot_path_elems = joinpath(filepath_dir,"post_analysis","plots")
        set_file(plot_path_elems)
        Plots.savefig(elem_conv_plt, joinpath(plot_path_elems,"conv_$dir.pdf"))
        Plots.savefig(elem_conv_plt_log, joinpath(plot_path_elems,"conv_log_$dir.pdf"))
        Plots.savefig(elem_ratio_plt, joinpath(plot_path_elems,"η_β_ratio_$dir.pdf"))
        Plots.savefig(elem_η_plt, joinpath(plot_path_elems,"η_$dir.pdf"))
        Plots.savefig(elem_β_plt, joinpath(plot_path_elems,"β_$dir.pdf"))
        Plots.savefig(elem_ratio_norm_plt, joinpath(plot_path_elems,"η_β_ratio_normalized_$dir.pdf"))
        Plots.savefig(elem_η_norm_plt, joinpath(plot_path_elems,"η_normalized_$dir.pdf"))
        Plots.savefig(elem_β_norm_plt, joinpath(plot_path_elems,"β_normalized_$dir.pdf"))
        Plots.savefig(elem_rel_height_error_plt, joinpath(plot_path_elems,"relative_height_error_$dir.pdf"))
        Plots.savefig(elem_height_plt, joinpath(plot_path_elems,"height_$dir.pdf"))
        @info "Saved plots to $plot_path_elems"
    end

    plot_path_global = joinpath(filepath,"post_analysis_global","plots")
    set_file(plot_path_global)

    @info "Saving plots to $plot_path_global"
    
    Plots.savefig(plot_conv_5, joinpath(plot_path_global,"conv_5.pdf"))
    Plots.savefig(plot_conv_log_5, joinpath(plot_path_global,"conv_log_5.pdf"))
    Plots.savefig(η_norm_plot_5, joinpath(plot_path_global,"η_normalized_5.pdf"))
    Plots.savefig(β_norm_plot_5, joinpath(plot_path_global,"β_normalized_5.pdf"))
    Plots.savefig(ratio_plot_5, joinpath(plot_path_global,"η_β_ratio_5.pdf"))
    Plots.savefig(ratio_norm_plot_5, joinpath(plot_path_global,"η_β_ratio_normalized_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"h_normalized_5.pdf"))
    Plots.savefig(h_glob_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(contour_plt, joinpath(plot_path_global,"contour_comparison_5.pdf"))
    Plots.savefig(contour_plt_zoom, joinpath(plot_path_global,"contour_comparison_5_zoomed.pdf"))
    Plots.savefig(contour_plt_zoom_05, joinpath(plot_path_global,"contour_comparison_zoom_05.pdf"))
    Plots.savefig(contour_plt_zoom_10, joinpath(plot_path_global,"contour_comparison_zoom_10.pdf"))
    Plots.savefig(contour_plt_zoom_15, joinpath(plot_path_global,"contour_comparison_zoom_15.pdf"))
    Plots.savefig(contour_plt_zoom_20, joinpath(plot_path_global,"contour_comparison_zoom_20.pdf"))
    Plots.savefig(cont_plt_legend, joinpath(plot_path_global,"slip_legend_vertical.pdf"))
    Plots.savefig(slip_case_legend, joinpath(plot_path_global,"slip_case_legend.pdf"))
    @info "Saved plots to $plot_path_global"
end

"""
    plot_contours(filepath_gt_, time_point=0)

Extract and plot ground-truth contours (one curve per numerically-named slip
case directory under `filepath_gt_`) at a single specified time point,
producing an overview and a zoomed comparison figure plus a shared slip-case
legend.

# Arguments
- `filepath_gt_::String`: path to ground truth data, containing numerically-named
  slip-case subdirectories.
- `time_point::Int`: time index to plot (`0` for the last time point).

# Returns
None. All outputs are PDFs written under `plot_path_global`.
"""
function plot_contours(filepath_gt_::String, time_point::Int=0)
    dir_list = readdir(filepath_gt_)

    cont_y_lims = [400, 1200] 
    cont_x_lims = [1380, 1430]
    
    contour_plt = set_plot(fs, sz=(round(Int,plt_width*1.5), plt_height), legend_column=1, legend=:outertopright, bottom_margin=0mm, top_margin=0mm)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt, 480, 1520)
    Plots.xlabel!(contour_plt, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt, L"y\;\mathrm{[px]}")

    cnt_plt_width = plt_width*2
    contour_plt_zoom = set_plot(round(Int,fs*1.1),sz=(plt_width*2, plt_height), legend_column=2, legend=:outertopright, bottom_margin=1mm, top_margin=0mm)
    Plots.xticks!(contour_plt_zoom, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom, cont_x_lims[1], cont_x_lims[2])
    Plots.xlabel!(contour_plt_zoom, L"x\;\mathrm{[px]}")
    Plots.ylabel!(contour_plt_zoom, L"y\;\mathrm{[px]}")

    slip_case_legend = set_plot(11, sz=(round(Int, plt_width*1.7), 50), legend=:bottom, legend_column=3, bottom_margin=-35mm, top_margin=2mm, left_margin=-25mm, right_margin=-15mm)
    Plots.plot!(slip_case_legend, [0, 1], [0, 0], label=false, color=:white, linewidth=0)
    Plots.xlims!(slip_case_legend, -0.2, 1.2)
    Plots.ylims!(slip_case_legend, -0.5, 0.5)
    
    global color_palette

    # Sort directories numerically
    numeric_dirs = sort([d for d in dir_list if tryparse(Int, d) !== nothing], by=x->parse(Int, x))
    println("Found directories in ground truth path: ", numeric_dirs)
    for (idx, dir) in enumerate(numeric_dirs)
        
        filepath_gt = joinpath(filepath_gt_, dir)
        if !isdir(filepath_gt)
            continue
        end

        printstyled("Processing slip case: $(dir) at index $idx\n", color=:green)
        
        sim_params_path = joinpath(filepath_gt, "data", "sim_params")
        
        sim_params = read_json(sim_params_path)
        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        
        println("Ground truth η: ", η_gt[1], ", β: ", β_gt[1])

        Plots.plot!(slip_case_legend, [0, 0.2], [0, 0], 
                    label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), 
                    color=color_palette[(idx-1) % length(color_palette)+1], linewidth=2, markerstrokewidth=0)

        contour_data_path = joinpath(filepath_gt, "data", "sim_data", "contour_data")
        try
            ObsDataList, splinexObs, splineyObs = read_csv(contour_data_path)
            
            plot_idx = time_point == 0 ? length(splinexObs) : min(time_point, length(splinexObs))
            
            printstyled("Plotting contour at time index $plot_idx out of $(length(splinexObs))\n", color=:blue)
            
            Plots.plot!(contour_plt, splinexObs[plot_idx], splineyObs[plot_idx], 
                        label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), 
                        color=color_palette[(idx-1) % length(color_palette)+1])
            
            Plots.plot!(contour_plt_zoom, splinexObs[plot_idx], splineyObs[plot_idx], 
                        label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), 
                        color=color_palette[(idx-1) % length(color_palette)+1])
            
        catch err
            @warn "Failed to load contour data for $dir: $err"
            continue
        end
    end
    
    plot_path_global = joinpath(filepath_gt_, "post_analysis_global", "plots")
    set_file(plot_path_global)

    @info "Saving contour plots to $plot_path_global"
    
    Plots.savefig(contour_plt, joinpath(plot_path_global, "gt_contour_comparison.pdf"))
    Plots.savefig(contour_plt_zoom, joinpath(plot_path_global, "gt_contour_comparison_zoom.pdf"))
    Plots.savefig(slip_case_legend, joinpath(plot_path_global, "slip_case_legend.pdf"))
    
    @info "Saved contour plots to $plot_path_global"
end

"""
    post_analysis_bulk(filepath_gt_, filepath, avoid_list)

Cross-experiment post-analysis and figure generation for `viscosity_type ==
"bulk_viscosity"` results — the bulk-viscosity counterpart of
`post_analysis_const`. Walks the same nested `dir/elem_size_folder/
simtime_folder/window` result tree under `filepath` (skipping entries in
`avoid_list`, `"post_analysis_global"`, and a few hardcoded element-size
folders), and overlays each leaf's time-windowed η/β estimates, their
ratio/product/sum, and height/relative-height-error against the ground truth
in `filepath_gt_`, at three nested levels:

- **global** (across all top-level `dir`s): normalized η/β/ratio, product,
  sum, and height/normalized-height/relative-error, bucketed by
  simulation-time case (5/10 s).
- **per-`dir`** (across element sizes): mesh-convergence plots of η, β, their
  ratio, and height/relative-height-error.
- **per-`elem_size_folder`** (across simulation-time windows): the same
  comparison plots, plus cost-surface slice-direction plots and (per
  simulation-time folder) height-error plots aggregated across windows.

# Arguments
- `filepath_gt_::String`: ground-truth data root each experiment is compared against.
- `filepath::String`: root of the result tree to walk.
- `avoid_list`: top-level `dir` names under `filepath` to skip.

# Returns
None. All outputs are PDFs written under `plot_path_global`/
`plot_path_elems`/`plot_path_sim_time`/`plot_path_noise`.
"""
function post_analysis_bulk(filepath_gt_::String, filepath::String, avoid_list)
    dir_list = readdir(filepath)


    η_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(η_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_5,L"\eta_{\mathrm{est}}/\eta^{\mathrm{avg}}_{\mathrm{gt}}(t)")

    β_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(β_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_5, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_5, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(η_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_10,L"\eta_{\mathrm{est}}/\eta^{\mathrm{avg}}_{\mathrm{gt}}(t)")

    β_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(β_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_10, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_10, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    h_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_plot_5,L"h_{\mathrm{est}}\;\mathrm{[mm]}")  
    Plots.xlims!(h_plot_5, 0, end_obs_win)

    h_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_plot_10,L"h_{\mathrm{est}}\;\mathrm{[mm]}") 
    Plots.xlims!(h_plot_10, 0, end_obs_win)

    h_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_5,L"h_{\mathrm{est}}/h_{\mathrm{gt}}")  
    Plots.xlims!(h_norm_plot_5, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_5, y_lims_h_norm)

    h_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_10,L"h_{\mathrm{est}}/h_{\mathrm{gt}}") 
    Plots.xlims!(h_norm_plot_10, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_10, y_lims_h_norm)

    rel_height_error_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_5, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_5, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_5, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_5, y_lims_rel_error)

    rel_height_error_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_10, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_10, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.xlims!(rel_height_error_glob_plot_10, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_10, y_lims_rel_error)

    ratio_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(ratio_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_5,L"\eta_{\mathrm{est}}/\beta_{\mathrm{est}}")

    ratio_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(ratio_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(ratio_plot_10,L"\eta_{\mathrm{est}}/\beta_{\mathrm{est}}")

    product_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(product_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(product_plot_5,L"\eta_{\mathrm{est}} \cdot \beta_{\mathrm{est}}")

    sum_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(sum_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(sum_plot_5,L"\eta_{\mathrm{est}} + \beta_{\mathrm{est}}")
    
    ratio_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(ratio_norm_plot_5, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_5,L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

    ratio_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(ratio_norm_plot_10, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_10,L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")
    
    η_y_max = 0
    β_y_max = 0
    ratio_y_max = 0
    for dir in dir_list
        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(joinpath(filepath_gt,"data","sim_params"))

        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)

        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        elem_size_folders = readdir(filepath_dir)

        elem_η_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_η_plt, [η_gt[1]], label="Ground truth η",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_η_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_plt,L"\eta\;\mathrm{[kPa\, s]}")

        elem_β_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_β_plt, [β_gt[1]], label="Ground truth β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_β_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_plt, L"\beta\;\mathrm{[Pa\, s \, m]}")

        elem_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_ratio_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_plt, L"\eta/\beta\;\mathrm{[mm^{-1}]}")
        
        elem_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        Plots.xlabel!(elem_η_norm_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_norm_plt,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

        elem_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        Plots.xlabel!(elem_β_norm_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_norm_plt, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

        elem_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        Plots.xlabel!(elem_ratio_norm_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_norm_plt, L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

        elem_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.plot!(elem_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_rel_height_error_plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(elem_rel_height_error_plt, L"\mathrm{Relative\;Height\;Error}")
        Plots.xlims!(elem_rel_height_error_plt, 0, end_obs_win)

        elem_height_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.plot!(elem_height_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_height_plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(elem_height_plt, L"\mathrm{Height\;[mm]}")
        Plots.xlims!(elem_height_plt, 0, end_obs_win)
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2" || elem_size_folder_ == "Q2_4" || elem_size_folder_ == "Q2_8"
                continue
            end
            
            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)
            
            sim_window_η_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_plt, [η_gt[1]], label="Ground truth η")
            Plots.xlabel!(sim_window_η_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_plt,latexstring("\$\\eta\$ [kPa s]"))
            
            sim_window_β_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_plt, [β_gt[1]], label="Ground truth β")
            Plots.xlabel!(sim_window_β_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_plt, latexstring("\$\\beta\$ [MPa s m\$^{-1}\$]"))

            sim_window_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
            Plots.xlabel!(sim_window_ratio_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_plt, L"\eta/\beta\;\mathrm{[mm^{-1}]}")

            sim_window_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            Plots.xlabel!(sim_window_η_norm_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_norm_plt,L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

            sim_window_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            Plots.xlabel!(sim_window_β_norm_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_norm_plt, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

            sim_window_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            Plots.xlabel!(sim_window_ratio_norm_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_norm_plt, L"\frac{\eta_{\mathrm{est}}/\beta_{\mathrm{est}}}{\eta_{\mathrm{gt}}/\beta_{\mathrm{gt}}}")

            sim_window_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_rel_height_error_plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(sim_window_rel_height_error_plt, L"\mathrm{Relative\;Height\;Error}")
            Plots.xlims!(sim_window_rel_height_error_plt, 0, end_obs_win)

            sim_window_height_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_height_plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(sim_window_height_plt, L"\mathrm{Height\;[mm]}")
            Plots.xlims!(sim_window_height_plt, 0, end_obs_win)

            plt_slices = set_plot(fs, sz=(350, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum",  lw=1)
            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;[px]}")
            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
            Plots.ylims!(plt_slices, 0, 50)

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" || sim_time_folder_ == "simtime_30.0"
                    continue
                end
                
                height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(height_error_plt, [], label=false)
                Plots.xlabel!(height_error_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(height_error_plt, L"\mathrm{Height\;Error\;[mm]}")
                Plots.xlims!(height_error_plt, 0, end_obs_win)
                
                rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(rel_height_error_plt, [], label=false)
                Plots.xlabel!(rel_height_error_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(rel_height_error_plt, L"\mathrm{Height\;Error\;[mm]}")
                Plots.xlims!(rel_height_error_plt, 0, end_obs_win)

                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_, "noise_0.0")

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                window_dirs = readdir(sim_time_folder)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        println("Skipping directory: $window_dir")
                        continue
                    end
                    win_exp_path = joinpath(filepath, elem_size_folder, sim_time_folder, window_dir)
                    
                    println("Processing window: $win_exp_path")
                    exp_params = read_json(joinpath(win_exp_path ,"data","experiment_parameters"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    noise_level = exp_params["noise_level"]
                    ne = exp_params["ne_exp"]

                    if data_type  == "synthetic"
                        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                        @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                    elseif data_type == "simulated"
                        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
                        @info "Data type $data_type Reading simulated contour data from $(joinpath(filepath_gt,"data","sim_data","contour_data")) of $(length(ObsDataList)) time steps"
                    else
                        error("Unknown data type: $data_type")
                    end

                    obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                    local η_gt::Vector{Float64} = float.(sim_params["η"])
                    local β_gt  = sim_params["β"]

                    est_ηpList = readdlm(joinpath(win_exp_path,"data","est_η.csv"), ',', Float64)
                    est_βpList = readdlm(joinpath(win_exp_path,"data","est_β.csv"), ',', Float64)
                    avg_ηList = readdlm(joinpath(win_exp_path,"data","avg_η.csv"), ',', Float64)
                    est_h_list = readdlm(joinpath(win_exp_path,"data","est_h.csv"), ',', Float64)
                    
                    data_ranges_ = get_time_windows(joinpath(win_exp_path,"data","window_data","data_ranges.csv"))
                    t_windows = readdlm(joinpath(win_exp_path,"data","window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(joinpath(win_exp_path,"data","window_data","time_windows.csv"),',',Float64)

                    sim_border_pt_lst, splinex, spliney = read_csv(joinpath(win_exp_path,"data","sim_data","2D_border_points"))
                    
                    println("Time windows: $(time_windows)")
                    obs_time = sum(time_windows)
                    
                    if obs_time != sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                        min_time = min(obs_time, sim_time)
                        sim_time = min_time
                        obs_time = min_time
                    end
                    
                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    η_gt = η_gt[1:data_point_len]
                    gt_h = gt_h_[1:(data_point_len+1)]
                    sim_border_pt_lst = sim_border_pt_lst[1:data_point_len+1, :]
                    obs_border_pt_lst = obs_border_pt_lst[1:data_point_len+1, :]
                    est_h = est_h_list[1:data_point_len+1, :]

                    height_error = abs.(est_h - gt_h)
                    rel_height_error = height_error ./ gt_h * 100.0

                    # normalize by ground-truth                      
                    est_η_norm = est_ηpList./avg_ηList
                    est_β_norm = est_βpList./β_gt

                    ratio_est = est_ηpList ./ est_βpList
                    ratio_gt = η_gt ./ β_gt
                    normalized_ratio = est_η_norm ./ est_β_norm
            
                    Plots.plot!(sim_window_η_plt, est_ηpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_plt, est_βpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_ratio_plt, ratio_est, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    
                    Plots.plot!(sim_window_rel_height_error_plt, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, gt_h, label=L"h_{\mathrm{gt}}(t)", legend=:outerbottom, legend_column=2)
                    
                    Plots.plot!(sim_window_η_norm_plt, est_η_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_norm_plt, est_β_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_ratio_norm_plt, normalized_ratio, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                    if sim_time_exp == 5.0
                        if ne == 6
                            η_y_max = max(η_y_max, maximum(est_η_norm))
                            β_y_max = max(β_y_max, maximum(est_β_norm))
                            ratio_y_max = max(ratio_y_max, maximum(normalized_ratio))
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))

                                Plots.vline!(η_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(β_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(ratio_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(ratio_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)

                                if ti == 1
                                    Plots.plot!(η_norm_plot_5, t_win, est_η_norm[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(β_norm_plot_5, t_win, est_β_norm[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), yscale=:log10, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(product_plot_5, t_win, est_ηpList[data_range_].*est_βpList[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(sum_plot_5, t_win, est_ηpList[data_range_].+est_βpList[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(ratio_norm_plot_5, t_win, normalized_ratio[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                else
                                    Plots.plot!(η_norm_plot_5, t_win, est_η_norm[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(β_norm_plot_5, t_win, est_β_norm[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=false, yscale=:log10, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(product_plot_5, t_win, est_ηpList[data_range_].*est_βpList[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(sum_plot_5, t_win, est_ηpList[data_range_].+est_βpList[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(ratio_norm_plot_5, t_win, normalized_ratio[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                end
                                t_prev = t+t_steps
                            end
                            Plots.ylims!(η_norm_plot_5, 0, max(η_y_max*1.1, 2))
                            Plots.ylims!(β_norm_plot_5, 0, max(β_y_max*1.1, 2))
                            Plots.ylims!(ratio_norm_plot_5, 0, max(ratio_y_max*1.1, 2))

                            Plots.plot!(h_plot_5, time, est_h, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(h_plot_5, time, gt_h, label=false, linestyle=:dash, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])

                            t = collect(range(start=0, stop=sim_time, step=t_steps))
                            Plots.plot!(rel_height_error_glob_plot_5, t, rel_height_error, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(h_norm_plot_5, t, est_h./gt_h, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            if dir == "1"
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t_win = t_windows[ti]
                                    data_range_ = data_ranges_[ti]

                                    Plots.vline!(h_norm_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                end
                            end
                            
                        end

                        
                        
                    
                    
                    elseif sim_time == 10.0
                        if ne == 6
                            Plots.plot!(η_norm_plot_10, est_η_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(β_norm_plot_10, est_β_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            
                            Plots.plot!(ratio_plot_10, ratio_est, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, yscale=:log10, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.hline!(ratio_plot_10, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)

                            Plots.plot!(h_plot_10, time, est_h, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(h_plot_10, time, gt_h, label=false, linestyle=:dash, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(rel_height_error_glob_plot_10, time, rel_height_error, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(ratio_norm_plot_10, normalized_ratio, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                        end
                    end
                end
                plot_path_noise = joinpath(sim_time_folder,"post_analysis_noise","plots")
                set_file(plot_path_noise)
                @info "Saving plots to $plot_path_noise"
                Plots.savefig(height_error_plt, joinpath(plot_path_noise,"height_error.pdf"))
                Plots.savefig(rel_height_error_plt, joinpath(plot_path_noise,"relative_height_error.pdf"))
                @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = joinpath(elem_size_folder,"post_analysis_time","plots")
            set_file(plot_path_sim_time)
            
            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β.pdf"))
            Plots.savefig(sim_window_ratio_plt, joinpath(plot_path_sim_time,"η_β_ratio.pdf"))
            Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height.pdf"))
            Plots.savefig(sim_window_η_norm_plt, joinpath(plot_path_sim_time,"η_normalized.pdf"))
            Plots.savefig(sim_window_β_norm_plt, joinpath(plot_path_sim_time,"β_normalized.pdf"))
            Plots.savefig(sim_window_ratio_norm_plt, joinpath(plot_path_sim_time,"η_β_ratio_normalized.pdf"))
            Plots.savefig(plt_slices, joinpath(plot_path_sim_time,"cost_slices_along_directions.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
        
        plot_path_elems = joinpath(filepath_dir,"post_analysis","plots")
        set_file(plot_path_elems)

        Plots.savefig(elem_ratio_plt, joinpath(plot_path_elems,"η_β_ratio.pdf"))
        Plots.savefig(elem_η_plt, joinpath(plot_path_elems,"η.pdf"))
        Plots.savefig(elem_β_plt, joinpath(plot_path_elems,"β.pdf"))
        Plots.savefig(elem_ratio_norm_plt, joinpath(plot_path_elems,"η_β_ratio_normalized.pdf"))
        Plots.savefig(elem_η_norm_plt, joinpath(plot_path_elems,"η_normalized.pdf"))
        Plots.savefig(elem_β_norm_plt, joinpath(plot_path_elems,"β_normalized.pdf"))
        Plots.savefig(elem_rel_height_error_plt, joinpath(plot_path_elems,"relative_height_error.pdf"))
        Plots.savefig(elem_height_plt, joinpath(plot_path_elems,"height.pdf"))
        @info "Saved plots to $plot_path_elems"
    end
    plot_path_global = joinpath(filepath,"post_analysis_global","plots")
    set_file(plot_path_global)
    
    Plots.savefig(η_norm_plot_5, joinpath(plot_path_global,"η_normalized_5.pdf"))
    Plots.savefig(β_norm_plot_5, joinpath(plot_path_global,"β_normalized_5.pdf"))
    Plots.savefig(ratio_plot_5, joinpath(plot_path_global,"η_β_ratio_5.pdf"))
    Plots.savefig(product_plot_5, joinpath(plot_path_global,"η_β_product_5.pdf"))
    Plots.savefig(sum_plot_5, joinpath(plot_path_global,"η_β_sum_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(ratio_norm_plot_5, joinpath(plot_path_global,"η_β_ratio_normalized_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"height_normalized_5.pdf"))

    Plots.savefig(η_norm_plot_10, joinpath(plot_path_global,"η_normalized_10.pdf"))
    Plots.savefig(β_norm_plot_10, joinpath(plot_path_global,"β_normalized_10.pdf"))
    Plots.savefig(ratio_plot_10, joinpath(plot_path_global,"η_β_ratio_10.pdf"))
    Plots.savefig(rel_height_error_glob_plot_10, joinpath(plot_path_global,"relative_height_error_10.pdf"))
    Plots.savefig(ratio_norm_plot_10, joinpath(plot_path_global,"η_β_ratio_normalized_10.pdf"))
    Plots.savefig(h_norm_plot_10, joinpath(plot_path_global,"height_normalized_10.pdf"))
    @info "Saved plots to $plot_path_global"
end

"""
    post_analysis_real(filepath_gt_, filepath, avoid_list)

Cross-experiment post-analysis and figure generation for physical (real,
measured) experiments — the counterpart of `post_analysis_bulk` for
`data_type == "physical"`, where no ground-truth η/β parameters exist, only a
measured height `h_m(t)`. Walks the same nested `dir/elem_size_folder/
simtime_folder/window` result tree under `filepath` (skipping entries in
`avoid_list`, `"post_analysis_global"`, a few hardcoded element-size folders,
and `simtime_2.0`), and overlays each leaf's time-windowed η(t)/β(t)
estimates, their ratio/product/sum, and height/relative-height-error against
the measured height in `filepath_gt_`, at three nested levels: global (across
all top-level `dir`s, bucketed by simulation-time case), per-`dir` (across
element sizes), and per-`elem_size_folder` (across simulation-time windows,
plus per-window height-error plots).

# Arguments
- `filepath_gt_::String`: ground-truth (measured height) data root each experiment
  is compared against.
- `filepath::String`: root of the result tree to walk.
- `avoid_list`: top-level `dir` names under `filepath` to skip.

# Returns
None. All outputs are PDFs written under `plot_path_global`/
`plot_path_elems`/`plot_path_sim_time`.
"""
function post_analysis_real(filepath_gt_::String, filepath::String, avoid_list)
    dir_list = readdir(filepath)

    η_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(η_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(η_plot_5,latexstring("\$\\eta_{\\mathrm{est}}(t)\$ [kPa s]"))

    β_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(β_plot_5, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(β_plot_5, latexstring("\$\\beta_{\\mathrm{est}}(t)\$ [MPa s m\$^{-1}\$]"))

    η_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(η_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(η_plot_10,latexstring("\$\\eta_{\\mathrm{est}}(t)\$ [kPa s]"))

    β_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(β_plot_10, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(β_plot_10, latexstring("\$\\beta_{\\mathrm{est}}(t)\$ [MPa s m\$^{-1}\$]"))

    gt_h_plot = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(gt_h_plot,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(gt_h_plot,L"h_{\mathrm{m}}\;\mathrm{[mm]}")  
    Plots.xlims!(gt_h_plot, 0, end_obs_win)

    h_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_plot_5,L"h\;\mathrm{[mm]}")  
    Plots.xlims!(h_plot_5, 0, end_obs_win)

    h_plot_est_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_est_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_plot_est_5,L"h_{\mathrm{est}}\;\mathrm{[mm]}")  
    Plots.xlims!(h_plot_est_5, 0, end_obs_win)

    h_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_plot_10,L"h_{\mathrm{est}}\;\mathrm{[mm]}") 
    Plots.xlims!(h_plot_10, 0, end_obs_win) 

    h_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_5,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_5,L"h_{\mathrm{est}}/h_{\mathrm{m}}")  
    Plots.xlims!(h_norm_plot_5, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_5, y_lims_h_norm)

    h_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(h_norm_plot_10,L"h_{\mathrm{est}}/h_{\mathrm{m}}") 
    Plots.xlims!(h_norm_plot_10, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_10, y_lims_h_norm)

    ratio_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(ratio_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_5,L"\eta_{\mathrm{est}}/\beta_{\mathrm{est}}")

    ratio_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(ratio_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(ratio_plot_10,L"\eta_{\mathrm{est}}/\beta_{\mathrm{est}}")

    product_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(product_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(product_plot_5,L"\eta_{\mathrm{est}} \cdot \beta_{\mathrm{est}}")

    product_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(product_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(product_plot_10,L"\eta_{\mathrm{est}} \cdot \beta_{\mathrm{est}}")

    sum_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(sum_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(sum_plot_5,L"\eta_{\mathrm{est}} + \beta_{\mathrm{est}}")

    sum_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(sum_plot_10,L"\mathrm{Time\;[s]}")
    Plots.ylabel!(sum_plot_10,L"\eta_{\mathrm{est}} + \beta_{\mathrm{est}}")

    rel_height_error_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_5, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_5, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.ylims!(rel_height_error_glob_plot_5, y_lims_rel_error)

    rel_height_error_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_10, L"\mathrm{Time\;[s]}")
    Plots.ylabel!(rel_height_error_glob_plot_10, latexstring("Relative Height Error [\$\\%\$]"))
    Plots.ylims!(rel_height_error_glob_plot_10, y_lims_rel_error)
    
    for dir in dir_list
        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(joinpath(filepath_gt,"data","sim_params"))
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        elem_size_folders = readdir(filepath_dir)
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2" || elem_size_folder_ == "Q2_4" || elem_size_folder_ == "Q2_8"
                continue
            end
            
            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)
            
            sim_window_η_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.xlabel!(sim_window_η_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_plt,L"\eta\;\mathrm{[kPa\, s]}")
            
            sim_window_β_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.xlabel!(sim_window_β_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_plt, L"\beta\;\mathrm{[Pa\, s \, m]}")

            sim_window_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.xlabel!(sim_window_ratio_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_plt, L"\eta/\beta\;\mathrm{[mm^{-1}]}")

            sim_window_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_rel_height_error_plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(sim_window_rel_height_error_plt, L"\mathrm{Relative\;Height\;Error}")
            Plots.xlims!(sim_window_rel_height_error_plt, 0, end_obs_win)

            sim_window_height_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_height_plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(sim_window_height_plt, L"\mathrm{Height\;[mm]}")
            Plots.xlims!(sim_window_height_plt, 0, end_obs_win)

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" #|| sim_time_folder_ == "simtime_10.0"
                    continue
                end
                
                height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(height_error_plt, [], label=false)
                Plots.xlabel!(height_error_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(height_error_plt, L"\mathrm{Height\;Error\;[mm]}")
                Plots.xlims!(height_error_plt, 0, end_obs_win)
                
                rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(rel_height_error_plt, [], label=false)
                Plots.xlabel!(rel_height_error_plt, L"\mathrm{Time\;[s]}")
                Plots.ylabel!(rel_height_error_plt, L"\mathrm{Height\;Error\;[mm]}")
                Plots.xlims!(rel_height_error_plt, 0, end_obs_win)

                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_, "noise_0.0")

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                window_dirs = readdir(sim_time_folder)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        println("Skipping directory: $window_dir")
                        continue
                    end
                    win_exp_path = joinpath(filepath, elem_size_folder, sim_time_folder, window_dir)
                    
                    println("Processing window: $win_exp_path")
                    exp_params = read_json(joinpath(win_exp_path ,"data","experiment_parameters"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    ne = exp_params["ne_exp"]

                    ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                    obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                    est_ηpList = readdlm(joinpath(win_exp_path,"data","est_η.csv"), ',', Float64)
                    est_βpList = readdlm(joinpath(win_exp_path,"data","est_β.csv"), ',', Float64)
                    est_h_list = readdlm(joinpath(win_exp_path,"data","est_h.csv"), ',', Float64)
                    
                    data_ranges_ = get_time_windows(joinpath(win_exp_path,"data","window_data","data_ranges.csv"))
                    t_windows = readdlm(joinpath(win_exp_path,"data","window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(joinpath(win_exp_path,"data","window_data","time_windows.csv"),',',Float64)

                    sim_border_pt_lst, splinex, spliney = read_csv(joinpath(win_exp_path,"data","sim_data","2D_border_points"))
                    
                    println("Time windows: $(time_windows)")
                    obs_time = sum(time_windows)
                    
                    if obs_time != sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                        min_time = min(obs_time, sim_time)
                        sim_time = min_time
                        obs_time = min_time
                    end
                    
                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    gt_h = gt_h_[1:(data_point_len+1)]
                    sim_border_pt_lst = sim_border_pt_lst[1:data_point_len+1, :]
                    obs_border_pt_lst = obs_border_pt_lst[1:data_point_len+1, :]
                    est_h = est_h_list[1:data_point_len+1, :]
                    time = collect(range(start=0, stop=sim_time, step=t_steps))

                    height_error = abs.(est_h - gt_h)
                    rel_height_error = height_error ./ gt_h * 100.0

                    ratio_est = est_ηpList ./ est_βpList
            
                    Plots.plot!(sim_window_η_plt, est_ηpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_plt, est_βpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_ratio_plt, ratio_est, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    
                    Plots.plot!(sim_window_rel_height_error_plt, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time[1:length(est_h)], est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, gt_h, label=L"h_{\mathrm{gt}}(t)", legend=:outerbottom, legend_column=2)
                    
                    if sim_time_exp == 0.5
                        if ne == 6
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))

                                Plots.vline!(η_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(β_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(ratio_plot_5, [t], color=:gray, linestyle=:dash, label=false)

                                if ti == 1
                                    Plots.plot!(η_plot_5, t_win, est_ηpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(β_plot_5, t_win, est_βpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(product_plot_5, t_win, est_ηpList[data_range_].*est_βpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(sum_plot_5, t_win, est_ηpList[data_range_].+est_βpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                else
                                    Plots.plot!(η_plot_5, t_win, est_ηpList[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(β_plot_5, t_win, est_βpList[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(product_plot_5, t_win, est_ηpList[data_range_].*est_βpList[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                    Plots.plot!(sum_plot_5, t_win, est_ηpList[data_range_].+est_βpList[data_range_], label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                                end
                                t_prev = t+t_steps
                            end
                            
                            Plots.plot!(rel_height_error_glob_plot_5, time, rel_height_error, label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(h_norm_plot_5, time, est_h./gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(h_plot_est_5, time, est_h, label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(gt_h_plot, time, gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            Plots.plot!(h_plot_5, time, gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1], linestyle=:dash)
                            Plots.plot!(h_plot_5, time, est_h, label=false, color=color_palette[(parse(Int, dir)-1) % length(color_palette)+1])
                            if dir == "1"
                                for ti::Int in 1:(size(data_ranges_, 1))
                                    t_win = t_windows[ti]
                                    data_range_ = data_ranges_[ti]
                                    Plots.vline!(h_norm_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(h_plot_est_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(h_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                end
                            end
                        
                        end

                        Plots.plot!(height_error_plt, time, height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(rel_height_error_plt, time, rel_height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        
                    elseif sim_time == 10.0
                        if ne == 6
                            Plots.plot!(η_norm_plot_10, est_η_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                            Plots.plot!(β_norm_plot_10, est_β_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                            
                            Plots.plot!(ratio_plot_10, ratio_est, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, legend=:outerbottom, legend_column=3, yscale=:log10)
                            Plots.hline!(ratio_plot_10, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)

                            Plots.plot!(rel_height_error_glob_plot_10, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                            Plots.plot!(ratio_norm_plot_10, normalized_ratio, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                        end
                    end
                end
            end
            plot_path_sim_time = joinpath(elem_size_folder,"post_analysis_time","plots")
            set_file(plot_path_sim_time)
            
            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β.pdf"))
            Plots.savefig(sim_window_ratio_plt, joinpath(plot_path_sim_time,"η_β_ratio.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
    end
    plot_path_global = joinpath(filepath,"post_analysis_global","plots")
    set_file(plot_path_global)
    
    Plots.savefig(η_plot_5, joinpath(plot_path_global,"η_plot_5.pdf"))
    Plots.savefig(β_plot_5, joinpath(plot_path_global,"β_plot_5.pdf"))
    Plots.savefig(ratio_plot_5, joinpath(plot_path_global,"η_β_ratio_5.pdf"))
    Plots.savefig(product_plot_5, joinpath(plot_path_global,"η_β_product_5.pdf"))
    Plots.savefig(sum_plot_5, joinpath(plot_path_global,"η_β_sum_5.pdf"))
    Plots.savefig(h_plot_est_5, joinpath(plot_path_global,"height_5.pdf"))
    Plots.savefig(h_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"height_normalized_5.pdf"))
    Plots.savefig(gt_h_plot, joinpath(plot_path_global,"gt_height.pdf"))
    @info "Saved plots to $plot_path_global"
end

"""
    read_unstrcut_csv(file_path) -> Vector

Read a CSV where each row is a single comma-joined string of numbers (a
"ragged"/unstructured layout `readdlm` can't parse directly as a numeric
matrix), and parse each row into a `Vector{Float64}`.

# Arguments
- `file_path::String`: path to the CSV file.

# Returns
- `data_list::Vector`: one `Vector{Float64}` per row.
"""
function read_unstrcut_csv(file_path::String)
    data_list = []

    data = readdlm(file_path)
    data_len = size(data,1)
    data_iter = 1:data_len
    for vec_iter in data_iter
        s = data[vec_iter,1]
        vec = parse.(Float64, split(s, ','))
        push!(data_list, vec)
    end
    return data_list
end

"""
    _get_borders(data_type, filepath_gt, exp_path, num_exp_points; view_folder="view_1") -> (obs_border_pt_lst, sim_border_pt_lst, gt_Splinex, gt_Spliney, splinex, spliney)

Load and align the observed (ground-truth) and simulated 2D border/contour
data for one experiment, truncated to the first `num_exp_points` frames.

# Arguments
- `data_type::String`: `"synthetic"`, `"physical"`, or `"simulated"` — selects
  where the ground-truth contour data is read from under `filepath_gt`.
- `filepath_gt::String`: ground-truth data root.
- `exp_path::String`: experiment result directory containing the simulated
  `2D_border_points`.
- `num_exp_points::Int`: number of leading time frames to keep.

# Keyword Arguments
- `view_folder::String`: camera-view subdirectory to read ground-truth data
  from (default: `"view_1"`).

# Returns
- `obs_border_pt_lst`, `sim_border_pt_lst`: observed and simulated border
  points, truncated to `num_exp_points` frames.
- `gt_Splinex`, `gt_Spliney`, `splinex`, `spliney`: matching spline-sampled
  x/y coordinates for the observed and simulated contours.
"""
function _get_borders(data_type::String, filepath_gt::String, exp_path::String, num_exp_points::Int; view_folder::String="view_1")

    if data_type  == "synthetic"
        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data",view_folder,"contour_data"))
        @info "Data type $data_type Reading synthetic ground truth contour data of $(length(ObsDataList)) time steps"
    elseif data_type == "physical"
        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data",view_folder,"contour_data"))
        @info "Data type $data_type Reading physical contour data of $(length(ObsDataList)) time steps"
    elseif data_type == "simulated"
        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))
        @info "Data type : $data_type Reading simulated ground truth contour data from $(joinpath(filepath_gt,"data","sim_data",view_folder,"contour_data")) of $(length(ObsDataList)) time steps"
    else
        error("Unknown data type: $data_type")
    end

    obs_border_pt_lst, gt_Splinex, gt_Spliney, pd = add_noise(ObsDataList, nFactor=0.0)
    sim_border_pt_lst, splinex, spliney = read_csv(joinpath(exp_path,"data","sim_data","view_1","2D_border_points_est"))
    @info "Reading simulated contour data from $(joinpath(exp_path,"data","sim_data","view_1","2D_border_points_est")) of $(length(sim_border_pt_lst)) time steps"

    @assert length(ObsDataList) >= num_exp_points "Not enough observation border points: have $(length(ObsDataList)), need at least $num_exp_points"
    @assert length(sim_border_pt_lst) >= num_exp_points "Not enough simulation border points: have $(length(sim_border_pt_lst)), need at least $num_exp_points"

    obs_border_pt_lst = obs_border_pt_lst[1:num_exp_points, :]
    sim_border_pt_lst = sim_border_pt_lst[1:num_exp_points, :]

    gt_Splinex = gt_Splinex[1:num_exp_points, :]
    gt_Spliney = gt_Spliney[1:num_exp_points, :]

    splinex = splinex[1:num_exp_points, :]
    spliney = spliney[1:num_exp_points, :]

    return obs_border_pt_lst, sim_border_pt_lst, gt_Splinex, gt_Spliney, splinex, spliney
end

"""
    save_plotly_surface_html(filename, xg, yg, Z; kwargs...) -> Bool

Render `Z` as an interactive Plotly 3D surface (grid `xg` × `yg`) with
optional overlaid estimation-path and ground-truth scatter traces, and write
it to `filename` as a self-contained HTML file. Used to make the
η/β cost-surface plots (see `replot`) explorable interactively,
complementing the static PDF contour plots.

Loads `PlotlyJS` lazily (via `@eval using`) so it isn't a hard dependency of
this script; returns `false` (and warns) instead of throwing if `PlotlyJS` is
unavailable or rendering fails for any other reason, so callers can treat this
as a best-effort side output. Small z offsets are added to overlay traces to
avoid z-fighting with the surface.

# Arguments
- `filename::AbstractString`: destination path for the HTML file.
- `xg`, `yg`: grid coordinates the surface `Z` is defined on.
- `Z`: surface height values on the `xg` × `yg` grid.

# Keyword Arguments
- `xs`, `ys`, `zs`: estimation-path overlay coordinates (default: empty; `zs`
  is back-filled via `interp_z_at` if omitted or mismatched in length).
- `gt_x`, `gt_y`, `gt_z`: single ground-truth marker coordinates (default: `nothing`).
- `title`: plot title (default: `""`).
- `colormap::AbstractString`: surface colormap (default: `"Viridis"`).
- `surface_label::AbstractString`: legend label for the surface (default: `"Cost surface"`).
- `path_label::AbstractString`: legend label for the estimation-path overlay (default: `"Estimations"`).
- `gt_label::AbstractString`: legend label for the ground-truth marker (default: `"Ground truth"`).
- `path_color::AbstractString`: color of the estimation-path markers (default: `"red"`).
- `gt_color::AbstractString`: colormap of the ground-truth marker (default: `"heat"`).
- `path_marker_size::Int`: marker size for the estimation-path overlay (default: `14`).
- `gt_marker_size::Int`: marker size for the ground-truth marker (default: `18`).
- `x_label::AbstractString`, `y_label::AbstractString`: axis labels (default: `""`).
- `font_size::Int`: axis/label font size (default: `14`).
- `latex_labels::Bool`: render labels as LaTeX (default: `false`).
- `z_offset::Real`: z offset applied to overlay traces to avoid z-fighting (default: `0.0`).
- `surface_opacity::Real`: opacity of the surface (default: `1.0`).

# Returns
- `Bool`: `true` if the HTML file was written; `false` if `PlotlyJS` is
  unavailable or rendering failed.
"""
function save_plotly_surface_html(filename::AbstractString, xg, yg, Z; xs=AbstractVector[], ys=AbstractVector[], zs=AbstractVector[],
                                  title="", colormap="Viridis",
                                  surface_label::AbstractString = "Cost surface",
                                  path_label::AbstractString = "Estimations",
                                  gt_label::AbstractString = "Ground truth",
                                  path_color::AbstractString = "red",
                                  gt_color::AbstractString = "heat",
                                  path_marker_size::Int = 14,
                                  gt_marker_size::Int = 18,
                                  gt_x = nothing, gt_y = nothing, gt_z = nothing,
                                  x_label::AbstractString = "",
                                  y_label::AbstractString = "",
                                  font_size::Int = 14,
                                  latex_labels::Bool = false,
                                  z_offset::Real = 0.0,
                                  surface_opacity::Real = 1.0)
    try
        @eval using PlotlyJS
    catch err
        @warn "PlotlyJS not available; cannot save interactive HTML: $err"
        return false
    end
    try
        xg_arr, yg_arr, Zmat = _ensure_grid_matrix(xg, yg, Z)
    # Plotly expects z as 2D array with rows corresponding to y
    # Surface trace: do NOT include the surface in the legend so the HTML
    # matches the PDF (PDF legend only shows Estimations and Ground truth).
    # Set transparency so overlay lines/markers remain visible. Use a lower
    # opacity to reduce occlusion of the estimation path.
    surf_trace = PlotlyJS.surface(x = xg_arr, y = yg_arr, z = Zmat, colorscale = colormap,
                      name = surface_label, showlegend = false, opacity = surface_opacity, showscale = true)
    # compute a tiny z offset (if z_offset==0.0 we pick a conservative
    # default proportional to the Z range) so overlays sit slightly above
    # the surface and avoid z-fighting/occlusion in the browser.
    zrange = maximum(Zmat) - minimum(Zmat)
    # use a slightly larger default offset so overlay traces are visible
    # in WebGL renderers; user can override via z_offset kwarg
    eff_offset = Float64(z_offset == 0.0 ? 1e-3 * (zrange + eps(Float64)) : z_offset)

    # Build traces list and create Plot without relying on add_traces! helper
        traces = [surf_trace]
        # add overlays if provided
        if length(xs) > 0 && length(ys) > 0
            xsv = Float64.(collect(xs))
            ysv = Float64.(collect(ys))
            # If zs not provided or wrong length, try to compute via interpolation
            if length(zs) == length(xs)
                zsv = Float64.(collect(zs))
            else
                try
                    zsv = Float64.(interp_z_at(xsv, ysv, xg_arr, yg_arr, Zmat))
                catch err
                    @warn "Estimator zs length mismatch and interpolation failed: $err; skipping estimation trace"
                    zsv = Float64[]
                end
            end
            if length(zsv) == length(xsv) && length(zsv) > 0
                @info "Adding estimation trace: $(length(xsv)) points"
                # Plot estimation path as line+markers so the trajectory is visible.
                # Use a star marker and strong black outline so it stands out above the surface.
                # lift the points slightly above the surface to avoid z-fighting
                zsv = zsv .+ eff_offset
                scatter_trace = PlotlyJS.scatter3d(x = xsv, y = ysv, z = zsv, mode = "lines+markers",
                                                   name = path_label,
                                                   marker = PlotlyJS.attr(size = path_marker_size, color = path_color, opacity = 1.0, symbol = "star", line = Dict(:width => 1.5, :color => "black")),
                                                   line = PlotlyJS.attr(width = 6, color = path_color))
                push!(traces, scatter_trace)
                # Add a non-legend projected overlay slightly above the surface
                # as a visual aid so the path is always visible in the HTML viewer.
                try
                    proj_offset = 0.02 * (zrange + eps(Float64))
                    zsv_proj = zsv .+ proj_offset
                    proj_marker = PlotlyJS.attr(size = max(path_marker_size+6, 22), color = path_color, opacity = 1.0, symbol = "circle", line = Dict(:width => 2, :color => "black"))
                    proj_trace = PlotlyJS.scatter3d(x = xsv, y = ysv, z = zsv_proj, mode = "markers",
                                                    name = "", showlegend = false,
                                                    marker = proj_marker)
                    push!(traces, proj_trace)
                catch err
                    @warn "Failed to add projected overlay trace: $err"
                end
            else
                @warn "Estimation trace not added: zs length $(length(zsv)) does not match xs length $(length(xsv))"
            end
        end
        # Optional ground-truth point (single marker)
        if gt_x !== nothing && gt_y !== nothing && gt_z !== nothing
            gx = _ensure_scalar_float(gt_x)
            gy = _ensure_scalar_float(gt_y)
            gz = _ensure_scalar_float(gt_z)
            # lift ground-truth marker slightly as well so it is visible above the surface
            gz = gz + eff_offset
            gt_trace = PlotlyJS.scatter3d(x = [gx], y = [gy], z = [gz], mode = "markers",
                                          name = gt_label,
                                          marker = PlotlyJS.attr(size = gt_marker_size, color = gt_color, symbol = "star"))
            push!(traces, gt_trace)
        end
    # Construct the Plot from the traces vector (some PlotlyJS versions
    # expect a single array argument rather than varargs).
        # Utility: convert a small subset of LaTeX expressions to a plain
        # unicode string for axis titles so they render correctly in Plotly's
        # SVG text when MathJax cannot typeset SVG contents.
        function _tex_to_plain(s::AbstractString)
            mp = Dict("\\eta"=>"η", "\\beta"=>"β", "\\,"=>"·")
            out = s
            # unwrap \mathrm{...} -> inner text
            out = replace(out, r"\\mathrm\{([^}]*)\}" => s"\1")
            # replace common macros
            for (k,v) in mp
                out = replace(out, k => v)
            end
            # remove dollar signs and braces
            out = replace(out, "\$" => "")
            out = replace(out, '{' => "")
            out = replace(out, '}' => "")
            return out
        end

        # Layout: place legend inside the plotting area (top-left) so it does
        # not overlap the surface colorbar on the right. Provide a translucent
        # white background for readability and a thin border.
        legend_cfg = Dict(:x => 0.02,
                          :y => 0.98,
                          :xanchor => "left",
                          :yanchor => "top",
                          :bgcolor => "rgba(255,255,255,0.85)",
                          :bordercolor => "black",
                          :borderwidth => 1,
                          :font => Dict(:size => Int(round(font_size*0.6))))
    # reasonable margins to accommodate colorbar and axes labels
    margin_cfg = Dict(:l => 60, :r => 120, :t => 80, :b => 60)

    # Axis label configs to align with PDF labels (2D axes). Use a smaller
    # tick font so numbers don't appear oversized in the HTML.
    tick_font_size = max(8, Int(round(font_size*0.36)))
    xaxis_cfg = Dict(:title => x_label, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size))
    yaxis_cfg = Dict(:title => y_label, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size))
    # Scene (3D) axis labels — ensure 3D plot axes show the same titles.
    # Use plain (unicode) versions for SVG rendering; MathJax may still
    # typeset the LaTeX string if available.
    plain_x = _tex_to_plain(x_label)
    plain_y = _tex_to_plain(y_label)
    scene_cfg = Dict(:xaxis => Dict(:title => plain_x, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size)),
             :yaxis => Dict(:title => plain_y, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size)),
             :zaxis => Dict(:title => "Cost", :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size)))

        # Try to construct the Plot with a Layout first; if the PlotlyJS
        # version here doesn't accept that signature, fall back to creating
        # the Plot and calling relayout!.
        # Provide a default camera view that looks slightly down at the surface
        # so overlaid paths are visible by default.
        camera_cfg = Dict(:eye => Dict(:x => 1.2, :y => 1.2, :z => 0.9))
        scene_with_camera = merge(scene_cfg, Dict(:camera => camera_cfg))

        layout_obj = try
            PlotlyJS.Layout(title = title, legend = legend_cfg, margin = margin_cfg, autosize = true,
                            xaxis = xaxis_cfg, yaxis = yaxis_cfg, scene = scene_with_camera)
        catch _
          # Older versions may expect a Dict for layout
            Dict(:title => title, :legend => legend_cfg, :margin => margin_cfg, :autosize => true,
                 :xaxis => xaxis_cfg, :yaxis => yaxis_cfg, :scene => scene_with_camera)
        end

        plt = try
            PlotlyJS.Plot(traces, layout_obj)
        catch _
            plt = PlotlyJS.Plot(traces)
            try
                PlotlyJS.relayout!(plt, layout_obj)
            catch _
                # If relayout also fails, at least set the title if available
                if !isempty(title)
                    try
                        PlotlyJS.relayout!(plt, Dict(:title => title))
                    catch _
                    end
                end
            end
        end
        # Render the plot to a HTML string using the MIME show method which is
        # available across PlotlyJS/PlotlyBase versions, then write the string.
        html_str = sprint(io -> show(io, MIME("text/html"), plt))
        # If LaTeX labels requested, inject MathJax into the HTML so TeX is rendered
        if latex_labels
            mj_cfg = raw"""<script>window.MathJax = {tex: {inlineMath: [['$','$'], ['\(','\)']]}};</script>""" * "\n"
            mj_src = "<script src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js\"></script>\n"
            mj_run = "<script>if(window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise();</script>\n"
            inject = mj_cfg * mj_src * mj_run
            # Try to insert into <head> if present, else before <body>, else prepend
            if occursin("</head>", html_str)
                html_str = replace(html_str, "</head>" => inject * "</head>")
            elseif occursin("<body>", html_str)
                html_str = replace(html_str, "<body>" => "<head>" * inject * "</head><body>")
            else
                html_str = inject * html_str
            end
        end
        open(filename, "w") do io
            write(io, html_str)
        end
        return true
    catch err
        @warn "Failed to create/save PlotlyJS HTML: $err"
        return false
    end
end

"""
    interp_z_at(x, y, xg, yg, Z) -> Float64

Bilinearly interpolate `Z` (defined on the grid `xg` × `yg`, with rows
indexing `yg` and columns indexing `xg`) at the point `(x, y)`, clamping
`(x, y)` to the grid extents. Tolerates `Z` being 1D, transposed, or of a
non-`Float64` element type (e.g. JSON-backed arrays) by coercing/reshaping it
to the expected `(length(yg), length(xg))` orientation first.

# Arguments
- `x::Real`, `y::Real`: point to interpolate at.
- `xg::AbstractVector`, `yg::AbstractVector`: grid coordinates `Z` is defined on.
- `Z::AbstractMatrix`: surface values on the `xg` × `yg` grid.

# Returns
- `Float64`: the interpolated value at `(x, y)`.
"""
function interp_z_at(x::Real, y::Real, xg::AbstractVector, yg::AbstractVector, Z::AbstractMatrix)
    # Coerce grid vectors and Z to concrete Float64 arrays/matrix and
    # normalize shapes. This avoids errors when Z is an Adjoint, 1D Vec,
    # or contains Union element types (JSON-backed arrays).
    xg_arr = Float64.(collect(xg))
    yg_arr = Float64.(collect(yg))
    nx = length(xg_arr); ny = length(yg_arr)

    Zmat = try
        Array(Z) |> x -> Float64.(x)
    catch err
        error("Failed to coerce Z for interpolation: $err")
    end

    # If Zmat ended up 1D (or different shape), try to reshape/transpose to
    # match (ny, nx) where rows -> y and cols -> x.
    if ndims(Zmat) == 1
        if length(Zmat) == nx*ny
            Zmat = reshape(Zmat, ny, nx)
        else
            error("interp_z_at: Z is 1D with length=$(length(Zmat)) which is not nx*ny=$(nx*ny)")
        end
    elseif ndims(Zmat) == 2
        if size(Zmat,1) == ny && size(Zmat,2) == nx
        elseif size(Zmat,1) == nx && size(Zmat,2) == ny
            Zmat = Zmat'
        elseif prod(size(Zmat)) == nx*ny
            Zmat = reshape(vec(Zmat), ny, nx)
        else
            @warn "interp_z_at: Z matrix shape $(size(Zmat)) does not match (ny,nx)=($(ny),$(nx)). Attempting to continue but interpolation may be wrong."
        end
    else
        error("interp_z_at: Z has unexpected number of dimensions: $(ndims(Zmat))")
    end

    # clamp x,y to grid extents
    xcl = clamp(x, first(xg_arr), last(xg_arr))
    ycl = clamp(y, first(yg_arr), last(yg_arr))

    ix = searchsortedfirst(xg_arr, xcl)
    if ix == 1
        i0 = 1; i1 = 1; tx = 0.0
    elseif ix > nx
        i0 = nx; i1 = nx; tx = 0.0
    else
        i1 = ix; i0 = max(1, ix-1)
        x0 = xg_arr[i0]; x1 = xg_arr[i1]
        tx = x1==x0 ? 0.0 : (xcl - x0)/(x1 - x0)
    end

    iy = searchsortedfirst(yg, ycl)
    if iy == 1
        j0 = 1; j1 = 1; ty = 0.0
    elseif iy > ny
        j0 = ny; j1 = ny; ty = 0.0
    else
        j1 = iy; j0 = max(1, iy-1)
        y0 = yg_arr[j0]; y1 = yg_arr[j1]
        ty = y1==y0 ? 0.0 : (ycl - y0)/(y1 - y0)
    end

    # values at corners: note Z rows->y (j index), cols->x (i index)
    z00 = float(Zmat[j0, i0]); z10 = float(Zmat[j0, i1])
    z01 = float(Zmat[j1, i0]); z11 = float(Zmat[j1, i1])

    # bilinear interpolation
    z0 = (1-tx)*z00 + tx*z10
    z1 = (1-tx)*z01 + tx*z11
    z = (1-ty)*z0 + ty*z1
    return z
end

"""
    interp_z_at(xs, ys, xg, yg, Z) -> Vector{Float64}

Vectorized form of `interp_z_at(x, y, xg, yg, Z)`: bilinearly
interpolate `Z` at each `(xs[i], ys[i])` pair.

# Arguments
- `xs::AbstractVector{<:Real}`, `ys::AbstractVector{<:Real}`: points to interpolate at.
- `xg::AbstractVector`, `yg::AbstractVector`: grid coordinates `Z` is defined on.
- `Z::AbstractMatrix`: surface values on the `xg` × `yg` grid.

# Returns
- `Vector{Float64}`: the interpolated value at each `(xs[i], ys[i])`.
"""
function interp_z_at(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, xg::AbstractVector, yg::AbstractVector, Z::AbstractMatrix)
    if length(xs) == 1 && length(ys) == 1
        return interp_z_at(xs[1], ys[1], xg, yg, Z)
    end
    # If paired vectors, compute element-wise
    if length(xs) == length(ys)
        return [interp_z_at(xi, yi, xg, yg, Z) for (xi, yi) in zip(xs, ys)]
    end
    error("interp_z_at: expecting scalar x,y or vectors of equal length")
end

"""
    _ensure_scalar_float(x) -> Float64

Coerce `x` to a `Float64`, taking the first element if `x` is a collection
rather than a plain `Number` (defensive against JSON-backed values that may
come back as length-1 arrays instead of scalars).

# Arguments
- `x`: a `Number`, or a collection whose first element is a `Number`.

# Returns
- `Float64`: `x` (or its first element) as a `Float64`.
"""
function _ensure_scalar_float(x)
    if isa(x, Number)
        return Float64(x)
    end
    try
        c = collect(x)
        if length(c) >= 1
            return Float64(c[1])
        else
            error("_ensure_scalar_float: container is empty")
        end
    catch err
        error("_ensure_scalar_float: cannot convert input to scalar Float64: $err")
    end
end

"""
    set_time_window(time_step_len, data; method="linear", window_size=10.0, global_window_sz=30.0) -> (time_windows, windows, data_ranges, t_windows)

Split `data` (indexed along its first dimension by time) into successive time
windows whose cumulative end time grows according to `method`
(`"linear"`: `window_size*iter`; `"quadratic"`: `window_size*iter^2`;
`"exponential"`: `window_size*exp(3*(iter-1))`), stopping once the window end
reaches the end of `data` or `global_window_sz`.

# Arguments
- `time_step_len::Float64`: time step length of `data`.
- `data::AbstractArray`: array indexed along its first dimension by time.

# Keyword Arguments
- `method::String`: window growth schedule, `"linear"`, `"quadratic"`, or
  `"exponential"` (default: `"linear"`).
- `window_size::Float64`: base window size (default: `10.0`).
- `global_window_sz::Float64`: overall time budget the windows are truncated
  to (default: `30.0`).

# Returns
- `time_windows::Vector{Float64}`: duration of each window.
- `windows::Vector{AbstractArray}`: the slice of `data` in each window.
- `data_ranges::Vector{AbstractArray}`: the index range (relative to `data`)
  of each window.
- `t_windows::Vector{Float64}`: cumulative end time of each window.
"""
function set_time_window(time_step_len::Float64, data::AbstractArray; method::String="linear", window_size::Float64=10.0, global_window_sz::Float64=30.0)
    windows::Vector{AbstractArray} = Vector{AbstractArray}()
    time_windows::Vector{Float64} = Vector{Float64}()
    data_ranges::Vector{AbstractArray} = Vector{AbstractArray}()
    t_windows::Vector{Float64} = Vector{Float64}()

    function get_t_window(window_size::Float64, step_len::Float64, iter::Int, method::String)::Float64
        t_window = 0.0
        if method == "linear"
            t_window = round(window_size*iter, digits=1)
        elseif method == "quadratic"
            t_window = round(window_size*iter^2, digits=1)
        elseif method == "exponential"
            if iter == 1
                t_window = round(window_size*exp(0.5*(iter-1)), digits=1)
            else 
                δt = window_size*round(exp(3*(iter-1)), digits=1)
                if δt < window_size
                    @warn "Computed time window increment $δt is less than minimum Window $window_size; using minimum."
                    δt = window_size
                end
                t_window = round(window_size*exp(3*(iter-1)), digits=1)
            end
        end
        return t_window
    end

    iter::Int = 1
    start_point::Int = 1
    t_window_prev::Float64 = 0.0
    t_window_end::Float64 = get_t_window(window_size, time_step_len, iter, method)
    end_point::Int = round(Int,t_window_end*time_step_len)+1

    END_FLAG::Bool = false
    while true
            # If the computed end_point reaches or exceeds the data length, adjust.
            # Keep the original requested value for clearer messages.
            if end_point >= size(data, 1) || t_window_end >= global_window_sz
                if end_point == size(data, 1)
                    @info "Reached end of data at point $end_point."
                elseif end_point > size(data, 1)
                    requested_end = end_point
                    # adjust the time window end to the last available sample
                    t_window_end = round((size(data, 1)-1)/time_step_len, digits=1)
                    end_point = size(data, 1)
                    @warn "Requested end point $requested_end exceeds data size $(size(data,1)); adjusting to end of data (end_point=$end_point). Adjusted end time to $t_window_end seconds."
                elseif t_window_end >= global_window_sz
                    requested_end = t_window_end
                    # adjust the end point to match the global window size
                    t_window_end = global_window_sz
                    end_point = round(Int,t_window_end*time_step_len)+1
                    @warn "Requested time window end $requested_end seconds exceeds global window size $global_window_sz seconds; adjusting to global window size (end_point=$end_point)."
                end
                END_FLAG = true
            end

        data_range = start_point:end_point
        data_range_ = start_point:(end_point-1)
        
        println("Data frame : $data_range")
        println("time windows from : $t_window_prev to $t_window_end")
        println("time Window : $(t_window_end - t_window_prev) seconds")
        println("data length : $(size(data[data_range], 1))")
        println("data range size : $(length(data_range_))")
        println("----------")
        t_window_size = round(t_window_end - t_window_prev, digits=1)
        push!(time_windows, t_window_size)
        push!(windows, data[data_range])
        push!(data_ranges, data_range_)
        push!(t_windows, t_window_end)

        if END_FLAG == true
            break
        end
        iter = iter + 1
        start_point = end_point
        t_window_prev = t_window_end
        t_window_end = get_t_window(window_size, time_step_len, iter, method)
        end_point = round(Int,t_window_end*time_step_len)+1
    end
    return time_windows, windows, data_ranges, t_windows
end

"""
    display_batch_info(n_experiments, n_cores=Threads.nthreads()) -> Int

Log how `n_experiments` parallel experiments will be split across
`n_cores`-sized batches, and return the number of batches.

# Arguments
- `n_experiments::Int`: total number of experiments to run.
- `n_cores::Int`: number of parallel cores/batches (default: `Threads.nthreads()`).

# Returns
- `n_batches::Int`: number of batches `n_experiments` is split into.
"""
function display_batch_info(n_experiments::Int, n_cores::Int=Threads.nthreads())
    n_batches = ceil(Int, n_experiments / n_cores)
    experiments_per_batch = n_experiments ÷ n_batches
    remaining = n_experiments % n_batches
    
    @info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    @info "CPU CORES AVAILABLE: $n_cores cores"
    @info "TOTAL EXPERIMENTS: $n_experiments"
    @info "NUMBER OF BATCHES: $n_batches"
    if n_batches > 1
        @info "EXPERIMENTS PER BATCH: $experiments_per_batch (batch 1-$n_batches each have $(experiments_per_batch + (remaining > 0 ? 1 : 0)) experiments max)"
        @info "BATCH DISTRIBUTION: $remaining batch(es) with extra experiment(s)"
    end
    @info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    return n_batches
end

"""
    _handle_worker_error(err, i, params)

Log a failed parallel worker task (index `i`, its `params`, and the caught
exception `err` with backtrace) via `@error`, for use in a `catch` block.

# Arguments
- `err`: the caught exception.
- `i::Int`: index of the failed task.
- `params`: parameters passed to the failed task, logged for diagnosis.

# Returns
None.
"""
function _handle_worker_error(err, i::Int, params)
    bt = catch_backtrace()
    @error "Task $i failed" params exception=(err, bt)
end


"""
    _max_band_mean(z, tol) -> Float64

Mean of all entries of `z` within `tol` of `maximum(z)`.

Used instead of a hard `maximum` so that a flat edge/face (many points tying
for the max within meshing/floating-point tolerance) doesn't collapse to a
single, arbitrarily-chosen point, and doesn't produce a QoI that jumps
discretely between mesh nodes as the geometry deforms.

# Arguments
- `z::AbstractVector`: values to reduce.
- `tol::Float64`: band width around `maximum(z)` over which entries are averaged.

# Returns
- `Float64`: mean of all entries of `z` within `tol` of `maximum(z)`.
"""
function _max_band_mean(z::AbstractVector, tol::Float64)
    z_max = maximum(z)
    band = z[z .≥ z_max - tol]
    return sum(band) / length(band)
end

"""
    get_surface_qoi(surface_list, cam_pose, height; tol=1e-3) -> Vector{Float64}

Quantity of interest for surface comparison: the depth (z, along the camera's
viewing axis) of the furthest surface points, per frame.

Each surface is transformed into the camera frame via `project_to_camera_frame`,
then reduced with `_max_band_mean` so that flat faces/edges facing the
camera (e.g. a cubic mesh) are averaged over rather than reduced to a single
mesh vertex.

# Arguments
- `surface_list::AbstractArray`: list of 3×N surface point matrices, one per frame.
- `cam_pose::AbstractArray`: camera pose used to build the camera frame.
- `height::Float64`: current specimen height, passed through to `project_to_camera_frame`.

# Keyword Arguments
- `tol::Float64`: band width (in the same units as the mesh) around the max
  depth over which points are averaged (default: `1e-3`).

# Returns
- `qoi::Vector{Float64}`: one QoI value per frame.
"""
function get_surface_qoi(surface_list::AbstractArray, cam_pose::AbstractArray, height::Float64; tol::Float64=1e-3)
    qoi = Vector{Float64}(undef, length(surface_list))

    for (i, sim_surface) in enumerate(surface_list)
        transformed_sim_surface = project_to_camera_frame(sim_surface, cam_pose, height)
        qoi[i] = _max_band_mean(transformed_sim_surface[3, :], tol)
    end

    return qoi
end

"""
    optimize_sim(use_parallel=true)

Driver for the mesh-convergence study: builds one `optimize` call
(`exp_params`) per combination of mesh element count (`nz_list`), noise level,
simulation time, and time step for `data_type == "simulated"`, constant
viscosity, reading ground truth from `ground_truth/sim_data/Stokes` and
writing results under `experiments/sim_data/convergence_analysis/
stokes_convergence`. Dispatches the resulting parameter list to
`run_parallel_tasks` when `use_parallel`, otherwise runs `optimize`
sequentially, logging (not raising) any per-task failure via
`_handle_worker_error`.

# Arguments
- `use_parallel::Bool`: dispatch to `run_parallel_tasks` instead of running
  sequentially (default: `true`).

# Returns
None.
"""
function optimize_sim(use_parallel::Bool=true)

    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2
    nz_list = Union{Int,Float64}[6]
    mode::Symbol = :conv_exp_mesh  # :exp

    if mode == :conv_exp_mesh
        nz_list = Union{Int,Float64}[6]
    end
    dt_list = [0.1] 
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]
    window = "multi_window"
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    geometry::Symbol = :cylinder # :cylinder or :cube

    avoid_dirs = ["post_analysis_global", "1", "2", "6", "4", "5" , "7", "8"]
    for viscosity_type in viscosity_type_list
        _filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Stokes"), control, viscosity_type, "Hex2_16.0", string(geometry))
        if !isdir(_filepath_gt)
            @warn "Ground truth directory not found, skipping: $_filepath_gt"
            continue
        end
        dir_list = readdir(_filepath_gt)
        for dir in dir_list
            if dir in avoid_dirs
                continue
                println("Skipping dir $dir")
            end
            filepath_gt = if mode === :conv_exp_mesh
                            joinpath(resolve_data_path("ground_truth/sim_data/Stokes/force/constant/Hex_2/convergence_analysis/experiment_mesh_convergence_analysis/"))
                          else    
                            joinpath(_filepath_gt, dir)
                          end
            println("Ground truth directory: $filepath_gt")
            for nz in nz_list
                if nz == 6 && mode == :conv_exp_mesh
                    dt_list = [0.1]
                else
                    dt_list = [0.1]
                end
                if nz == 6.5 && viscosity_type == "constant"
                    noise_level_list = [0.5, 1.0, 1.5, 2.0]
                else
                    noise_level_list = [0.0]
                end
                for noise_level in noise_level_list 
                    sim_time_exp_list = [5.0]  # simulation time in seconds
                    println("Simulation time experiments to [run: $sim_time_exp_list")
                    for sim_time_exp::Float16 in sim_time_exp_list
                        if viscosity_type == "constant"
                            window = ""
                        end
                        @info "Running optimization with ne = $nz and simulation time = $sim_time_exp with noise level = $noise_level"
                        for dt in dt_list
                            filepath_res = if mode === :conv_exp_mesh
                                                joinpath(resolve_data_path(joinpath("experiments","sim_data","convergence_analysis","stokes_convergence")), "experiment_mesh_conv", string(geometry), dir, "$(element_shape_x)$(basis_order_x)_$(nz)", "dt_$(dt)")
                                            else
                                                joinpath(resolve_data_path(joinpath("experiments","sim_data","optimization","Stokes")), control, viscosity_type, "Hex2_16.0", string(geometry), dir, "$(element_shape_x)$(basis_order_x)_$(nz)", "simtime_$(sim_time_exp)", "noise_$(noise_level)", "dt_$(dt)", window)
                                            end
                                                @info "Running optimization with element_shape_x = $element_shape_x, basis_order_x = $basis_order_x, ne = $nz, dt = $dt"

                            exp_params = Dict(
                                "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                                "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                                "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                                "ne_exp" => nz, "sim_time_exp" => sim_time_exp, "dt" => dt,
                                "filepath_res" => filepath_res, "filepath_gt" => filepath_gt,
                                "control" => control, "data_type" => "simulated",
                                "camera_matrix" => camera_matrix, "WRITE_GT" => false,
                                "noise_level" => noise_level, "mode" => window)

                            push!(param_list, exp_params)
                        end
                    end
                end
            end
        end
        if use_parallel
            run_parallel_tasks(param_list, optimize; max_workers=15, memory_per_task_mb=8192.0)
        else
            for (i, params) in enumerate(param_list)
                @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
                try
                    optimize(params)
                    @info "Completed write_gt_data for index $i"
                catch err
                    _handle_worker_error(err, i, params)
                end
            end
            println("All experiments completed.")
        end
    end
end

"""
    optimize_syn(use_parallel=true)

Driver for synthetic (rendered-image) bulk-viscosity experiments: builds one
`optimize` call per combination of mesh element count (`ne_list`),
noise level, simulation time, and time step for `data_type == "synthetic"`,
`viscosity_type == "bulk_viscosity"`, reading ground truth from
`ground_truth/sim_data/Stokes` and writing results under
`experiments/syn_data/optimization/Stokes`. Dispatches the resulting
parameter list to `run_parallel_tasks` when `use_parallel`, otherwise runs
`optimize` sequentially, logging (not raising) any per-task failure via
`_handle_worker_error`.

# Arguments
- `use_parallel::Bool`: dispatch to `run_parallel_tasks` instead of running
  sequentially (default: `true`).

# Returns
None.
"""
function optimize_syn(use_parallel::Bool=true)

    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2

    ne_list = [6]
    noise_level_list = [0.0] # No noise for synthetic data
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["bulk_viscosity"]
    window = "multi_window" 
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    geometry::Symbol = :cube # :cylinder or :cube
    avoid_dirs = ["post_analysis_global", "1", "2", "3", "4", "6" , "7", "8"]
    dt_list = [0.1]

    for viscosity_type in viscosity_type_list
        _filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Stokes"), control, viscosity_type, "Hex2_16.0", string(geometry))
        if !isdir(_filepath_gt)
            @warn "Ground truth directory not found, skipping: $_filepath_gt"
            continue
        end
        dir_list = readdir(_filepath_gt)
        for dir in dir_list
            if dir in avoid_dirs
                continue
                println("Skipping dir $dir")
            end
            filepath_gt = joinpath(_filepath_gt, dir)
            for ne in ne_list
                for noise_level in noise_level_list
                    sim_time_exp_list = [5.0]
                    println("Simulation time experiments to run: $sim_time_exp_list")
                    for sim_time_exp::Float16 in sim_time_exp_list
                        if viscosity_type == "constant"
                            window = ""
                        end
                        @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level"
                        for dt in dt_list
                            filepath_res = joinpath(resolve_data_path(joinpath("experiments","syn_data","optimization","Stokes")), control, viscosity_type, "Hex2_16.0", string(geometry), dir, "$(element_shape_x)$(basis_order_x)_$(ne)", "simtime_$(sim_time_exp)", "noise_$(noise_level)", "dt_$(dt)", window)
                            @info "Running optimization with element_shape_x = $element_shape_x, basis_order_x = $basis_order_x with $ne elements"

                            exp_params = Dict(
                                "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                                "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                                "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                                "ne_exp" => ne, "sim_time_exp" => sim_time_exp,
                                "dt" => dt,
                                "filepath_res" => filepath_res, "filepath_gt" => filepath_gt, "control" => control, "data_type" => "synthetic", "camera_matrix" => camera_matrix, "WRITE_GT" => false,
                                "noise_level"=>noise_level, "mode"=>window)

                            push!(param_list, exp_params)
                        end
                    end
                end
            end
        end
        if use_parallel
            run_parallel_tasks(param_list, optimize; max_workers=15, memory_per_task_mb=8192.0)
        else
            for (i, params) in enumerate(param_list)
                @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
                try
                    optimize(params)
                    @info "Completed write_gt_data for index $i"
                catch err
                    _handle_worker_error(err, i, params)
                end
            end
            println("All experiments completed.")
        end
    end
end

"""
    optimize_real(use_parallel=true)

Driver for physically-measured (or Carreau-model synthetic) bulk-viscosity
experiments: builds one `optimize` call per ground-truth directory
and mesh element count (`ne_list`). When `model_type == "carreau"`, reads
ground truth from `ground_truth/sim_data/Carreau` and sets `data_type =
"synthetic"`; otherwise reads from `ground_truth/physical_data` and sets
`data_type = "physical"` (supplying the applied force `F_ext`, cylinder
radius `r`, and height `h`). Dispatches the resulting parameter list to
`run_parallel_tasks` when `use_parallel`, otherwise runs `optimize`
sequentially, logging (not raising) any per-task failure via
`_handle_worker_error`.

# Arguments
- `use_parallel::Bool`: dispatch to `run_parallel_tasks` instead of running
  sequentially (default: `true`).

# Returns
None.
"""
function optimize_real(use_parallel::Bool=true)
    F_ext::Float64 = 9.812*1e3 # force applied to the cylinder in N
    sim_time_exp_list::Vector{Float64} = [0.5] # simulation time in seconds

    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2

    ne_list = [6]
    control = "force" # "force" or "velocity"
    η_start = 50.0
    β_start = 50.0
    viscosity_type = "bulk_viscosity"
    
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'

    model_type = "Stokes" # "carreau" or "Stokes"
    window = "multi_window" # "multi_window" or "single_window"
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    avoid_dirs = ["post_analysis_global"] 

    if model_type == "carreau" && viscosity_type == "bulk_viscosity"
        _filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Carreau"))
        sim_time_exp_list = [5.0]
        F_ext = 1e4 # force applied to the cylinder in N
    else
        _filepath_gt = joinpath(resolve_data_path("ground_truth/physical_data"))
        sim_time_exp_list = [0.5] # simulation time in seconds
        F_ext = 9.812*1e3 # force applied to the cylinder in N
        η_start = 10.0
        β_start = 10.0
        h = 38.5  # height of the cylinder in mm
    end

    # Extract base paths to variables using joinpath for cross-platform compatibility
    data_root = joinpath(homedir(), "SMEAR-PhD", "SMEAR-DataFiles", "Data")
    base_experiments_path = joinpath(data_root, "experiments")
    base_gt_path = joinpath(data_root, "ground_truth")
    
    if !isdir(_filepath_gt)
        @warn "Ground truth directory not found, skipping optimize_real: $_filepath_gt"
        return
    end
    dir_list = readdir(_filepath_gt)
    for dir in dir_list
        if dir in avoid_dirs
            continue
            println("Skipping dir $dir")
        end
        @info "Processing ground truth directory: $dir for $viscosity_type viscosity ..."
        filepath_gt = joinpath(_filepath_gt, dir)
        for ne in ne_list
            println("Simulation time experiments to run: $sim_time_exp_list")
            for sim_time_exp in sim_time_exp_list
                println("Simulation time: $sim_time_exp seconds")
                for ne in ne_list
                    if model_type == "carreau" && viscosity_type == "bulk_viscosity"
                        filepath_res = joinpath(base_experiments_path, "syn_data", "optimization", "Carreau", dir, "$(element_shape_x)$(basis_order_x)_$(ne)", "simtime_$(sim_time_exp)", "noise_0.0", window)
                        exp_params = Dict(
                            "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                            "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                            "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                            "ne_exp" => ne, "sim_time_exp" => sim_time_exp,
                            "η_start" => η_start, "β_start" => β_start,
                            "filepath_res" => filepath_res, "filepath_gt" => filepath_gt,
                            "control" => control, "data_type" => "synthetic",
                            "camera_matrix" => camera_matrix, "WRITE_GT" => false, "noise_level" => 0.0, "mode" => window)
                    else
                        filepath_res = joinpath(base_experiments_path, "physical_data", "optimization", dir, "$(element_shape_x)$(basis_order_x)_$(ne)", "simtime_$(sim_time_exp)", "noise_0.0", window)
                        exp_params = Dict(
                            "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                            "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                            "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                            "ne_exp" => ne, "sim_time_exp" => sim_time_exp,
                            "η_start" => η_start, "β_start" => β_start,
                            "filepath_res" => filepath_res, "filepath_gt" => filepath_gt,
                            "control" => control, "viscosity_type" => viscosity_type,
                            "data_type" => "physical", "r" => r, "h" => h,
                            "camera_matrix" => camera_matrix, "F_ext" => F_ext, "mode" => window)
                    end
                    @info "Running optimization with element_shape_x = $element_shape_x, basis_order_x = $basis_order_x with $ne elements"

                    push!(param_list, exp_params)
                end
            end
        end
    end
    if use_parallel
        run_parallel_tasks(param_list, optimize; max_workers=5, memory_per_task_mb=8192.0)
    else
        for (i, params) in enumerate(param_list)
            @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
            try
                optimize(params)
                @info "Completed write_gt_data for index $i"
            catch err
                _handle_worker_error(err, i, params)
            end
        end
        println("All experiments completed.")
    end
end


"""
    set_plot_config(data_type, viscosity_type)

Set the module-global plot styling variables (`fs`, `plt_height`, `plt_width`,
margins, and the `y_lims_h_norm`/`y_lims_rel_error` axis limits used
throughout `replot`/`post_analysis_*`) from `PLOT_PRESETS` (defined near the
top of this file), keyed by `data_type ∈ ("physical", "synthetic",
"simulated")` — `"synthetic"` is further split by `viscosity_type` since it
uses a narrower plot for `"bulk_viscosity"`. Falls back to `PLOT_CONFIG`
(with a warning) for an unrecognized `data_type`/`viscosity_type` combo.

# Arguments
- `data_type::String`: `"physical"`, `"synthetic"`, or `"simulated"`.
- `viscosity_type::String`: `"constant"` or `"bulk_viscosity"`; only
  consulted when `data_type == "synthetic"`.

# Returns
None. The effect is entirely via the mutated globals.
"""
function set_plot_config(data_type::String, viscosity_type::String)
    global fs, plt_height, plt_width, plt_lft_margin, plt_right_margin, plt_top_margin
    global y_lims_h_norm, y_lims_rel_error

    preset_key = data_type == "synthetic" ? (data_type, viscosity_type) : data_type
    preset = get(PLOT_PRESETS, preset_key, nothing)
    if preset === nothing
        @warn "Unknown data_type/viscosity_type combo '$data_type'/'$viscosity_type', using PLOT_CONFIG default"
        return
    end

    fs = preset[:font_size]
    plt_height = preset[:plot_height]
    plt_width = preset[:plot_width]
    plt_lft_margin = preset[:left_margin]
    plt_right_margin = preset[:right_margin]
    plt_top_margin = preset[:top_margin]

    y_lims_h_norm = preset[:y_lims_h_norm]
    y_lims_rel_error = preset[:y_lims_rel_error]
end

"""
    plot_results()

Top-level driver that regenerates plots for every result tree under
`~/SMEAR-PhD/SMEAR-DataFiles/Data/experiments`, across the configured
`data_type_list` (`"simulated"`/`"synthetic"`/`"physical"`) and their
applicable `model_type` (`"Stokes"`/`"carreau"`) and `viscosity_type`
combinations. For each ground-truth directory (skipping `avoid_dirs`), applies
`set_plot_config` for that data/viscosity type and calls
`replot` on the matching experiment directory (`predict` and
the `post_analysis_*` cross-experiment summaries are wired in but currently
commented out below).

# Returns
None. All outputs are the side effects of `replot`.
"""
function plot_results()
    control::String = "force"
    viscosity_type_list = [] # "constant" or "bulk_viscosity"
    model_type = [] # "carreau" or "Stokes"
    avoid_dirs = ["post_analysis_global","1","2","3","4","6","7","8"] # directories to skip in post-analysis and plotting
    data_type_list = ["synthetic"] # ["simulated", "synthetic", "physical"]
    base_path = ""
    geometry::Symbol = :cube # :cylinder or :cube

    for data_type in data_type_list
        if data_type == "physical"
            model_type = ["Stokes"] # for physical data, we only have Stokes model results for now
        elseif data_type == "synthetic"
            model_type = ["carreau", "Stokes"]
        elseif data_type == "simulated"
            model_type = ["Stokes"]
        end

        data_root = joinpath(homedir(), "SMEAR-PhD", "SMEAR-DataFiles", "Data")
        base_experiments_path = joinpath(data_root, "experiments")
        base_gt_path = joinpath(data_root, "ground_truth")
        for model_type::String in model_type
            if data_type == "synthetic" 
                if model_type == "carreau"
                    base_path = joinpath(base_experiments_path, "syn_data", "optimization", "Carreau")
                    viscosity_type_list = ["bulk_viscosity"] # for Carreau model, we only have bulk viscosity results for now
                else
                    base_path = joinpath(base_experiments_path, "syn_data", "optimization", "Stokes")
                    viscosity_type_list = ["bulk_viscosity"]
                end
            elseif data_type == "simulated"
                base_path = joinpath(base_experiments_path, "sim_data", "optimization", "Stokes")
                viscosity_type_list = ["constant"]
            elseif data_type == "physical"
                base_path = joinpath(base_experiments_path, "physical_data", "optimization")
                viscosity_type_list = ["bulk_viscosity"]
            end
            for viscosity_type::String in viscosity_type_list
                set_plot_config(data_type, viscosity_type)
                if data_type == "physical"
                    filepath_gt = joinpath(base_gt_path, "physical_data")
                    filepath_res = base_path
                else
                    if model_type == "carreau" && viscosity_type == "bulk_viscosity"
                        filepath_gt = joinpath(base_gt_path, "sim_data", "Carreau")
                        filepath_res = base_path
                    else
                        filepath_gt = joinpath(base_gt_path, "sim_data", "Stokes", control, viscosity_type, "Hex2_16.0" , string(geometry))
                        filepath_res = joinpath(base_path, control, viscosity_type, "Hex2_16.0" , string(geometry))
                    end
                end
                
                if !isdir(filepath_gt)
                    @warn "Ground truth directory not found, skipping: $filepath_gt"
                    continue
                end
                dirs = readdir(filepath_gt)
                for dir::String in dirs
                    if dir in avoid_dirs
                        continue
                        println("Skipping dir $dir")
                    end
                    filepath_gt_dir = joinpath(filepath_gt, dir)
                    filepath_res_dir = joinpath(filepath_res, dir)
                    println("Processing ground truth directory: $filepath_gt_dir for $viscosity_type viscosity ...")
                    replot(filepath_res_dir, filepath_gt_dir)
                end
            end
        end
    end
end

optimize_sim(false)
# plot_results()
