using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function main()
    # parameters for the optimization
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne_gt::Int = 6 # number of elements in the mesh for the ground truth

    β_gt_list = [10.0, 50.0, 100.0, 1e3]
    η_gt_list = [100.0]

    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"

    viscosity_type_list = ["constant"] # "constant" or "bulk_viscosity"
    FunctionClass_x_gt_list = ["Q2"] # Function space for the ground truth

    F_ext_::Float64 = 0.15*9.812*1e3 # force applied to the cylinder in N
    sim_time_gt::Float64 = 20.0 # simulation time in seconds
    steps_gt::Float64 = 200.0 # number of time steps

    obj_pose = zeros(Float64, 4,4)
    obj_pose[1,1] = -1.0
    obj_pose[2,3] = -1.0
    obj_pose[3,2] = -1.0
    obj_pose[1:3,4] = [0.0, h/2, 150.0]
    camera_matrix::AbstractMatrix{Float64} = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'

    for viscosity_type in viscosity_type_list
        for FunctionClass_x in FunctionClass_x_gt_list
            run_id = 1
            for β_gt in β_gt_list
                for η_gt in η_gt_list
                    if β_gt <= 1.0
                        F_ext::Float64 = 9.813e3
                    elseif β_gt == 10.0
                        F_ext = 9.813e3*1.75 # force applied to the cylinder in kg.mm/s^2 (N)
                    elseif β_gt == 50.0
                        F_ext = 9.813e3*5.5 # force applied to the cylinder in kg.mm/s^2 (N)
                    elseif β_gt == 100.0
                        F_ext = 9.813e3*10 # force applied to the cylinder in kg.mm/s^2 (N)
                    elseif β_gt == 1e3
                        F_ext = 9.813e3*70 # force applied to the cylinder in kg.mm/s^2 (N)
                    end
                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$(FunctionClass_x_gt)_$(ne_gt)/$run_id")

                    exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "β_gt" => β_gt, 
                                "η_gt" => η_gt, "filepath_gt"=>filepath_gt, "control" => control, "viscosity_type" => viscosity_type, "obj_pose_gt" => obj_pose, 
                                "F_ext" => F_ext, "sim_time_gt" => sim_time_gt, "steps_gt" => steps_gt, "r" => r, "h" => h, "camera_matrix" => camera_matrix)

                    write_data(exp_params)
                    run_id = run_id + 1
                end
            end
        end
    end
end

main()