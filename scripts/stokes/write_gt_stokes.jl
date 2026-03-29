using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function main()
    # parameters for the optimization
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne_gt::Int = 16 # number of elements in the mesh for the ground truth

    β_gt_list = [0.01, 10.0, 50.0, 1e2, 500.0, 1e3, 5e3] # penalty parameters for the ground truth [2e3, 5e3, 1e4, 1e5, 1e10]
    η_gt_list = [1e2] # viscosity values for the ground truth in kg/(mm⋅s)

    control = "force" # "force" or "velocity"

    viscosity_type_list = ["constant"] # "constant" or "bulk_viscosity"
    FunctionClass_x_gt_list = ["Q2"] # Function space for the ground truth

    F_ext::Float64 = 9.813e3*0.85 
    
    sim_time_gt::Float64 = 30.0 # simulation time in seconds
    steps_gt::Int = 300 # number of time steps

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
                if β_gt <= 1.0
                    F_ext = 9.813e3*0.85
                elseif β_gt == 10.0
                    F_ext = 9.813e3*0.93
                elseif β_gt == 50.0
                    F_ext = 9.813e3*0.955
                elseif β_gt == 100.0
                    F_ext = 9.813e3*0.97
                elseif β_gt == 500.0
                    F_ext = 9.813e3*0.99
                elseif β_gt == 1e3
                    F_ext = 9.813e3*1.0
                elseif β_gt == 5e3
                    F_ext = 9.813e3*1.01 # force applied to the cylinder in kg.mm/s^2 (N)

                elseif β_gt == 2e3
                    F_ext = 9.813e3*0.995 # force applied to the cylinder in kg.mm/s^2 (N)
                elseif β_gt == 1e4
                    F_ext = 9.813e3*0.85*700 # force applied to the cylinder in kg.mm/s^2 (N)
                elseif β_gt == 1e5
                    F_ext = 9.813e3*0.85*7e3 # force applied to the cylinder in kg.mm/s^2 (N)
                elseif β_gt == 1e10
                    F_ext = 9.813e3*0.85*7e8 # force applied to the cylinder in kg.mm/s^2 (N)
                end
                for η_gt in η_gt_list
                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$(FunctionClass_x)_$(ne_gt)/$run_id")

                    exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "β_gt" => β_gt, 
                                "η_gt" => η_gt, "filepath_gt"=>filepath_gt, "control" => control, "viscosity_type" => viscosity_type, "obj_pose_gt" => obj_pose, 
                                "F_ext" => F_ext, "sim_time_gt" => sim_time_gt, "steps_gt" => steps_gt, "r" => r, "h" => h, "camera_matrix" => camera_matrix)

                    write_gt_data(exp_params)
                    run_id = run_id + 1
                end
            end
        end
    end
end

main()