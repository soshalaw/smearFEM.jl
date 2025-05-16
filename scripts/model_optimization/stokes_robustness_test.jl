using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
using Measures
using Dates
using DelimitedFiles

function main(βLst, noiseLevelLst, ηLst)

    # test case 
    r = 0.5
    h = 1.0
    ne = 10
    ndim = 3
    FunctionClass_u = "Q2"
    nDof_u = ndim  # number of degree of freedom per node
    FunctionClass_p = "Q1"
    nDof_p = 1  # number of degree of freedom per node

    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

    dateTime = Dates.now()
    sampleNo = 21
    dev = 0.1

    control = "force" # "force" or "velocity"
    viscosity_type = "constant" # "constant" or "bulk_viscosity"

    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",control,"/test3")

    # simulation parameters for the ground truth
    sim_time_gt = 120.0
    steps_gt = 120.0
    t_steps_gt = sim_time_gt/steps_gt

    sim_time = 10.0
    steps = 10.0
    t_steps = sim_time/steps

    sideList = [false]

    F = 3.0

    for η_gt::Float64 in ηLst
        dev_η::Float64 = dev*η_gt
        ηStart::Float64 = η_gt - dev_η
        β_gt::Int = 1
        for β_gt::Float64 in βLst[1:end]
            dev_β::Float64 = dev*β_gt
            βStart::Float64 = β_gt - dev_β

            println("Ground truth: η :", η_gt, " β :", β_gt)
            model_gt, scene_gt = def_problem(r, h, ne, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β_gt, F, control, viscosity_type, sim_time_gt, t_steps_gt)
            filepath_gt = string(filepath,"/experiment_$(η_gt)_$(β_gt)/ground_truth")
            write_sim_data(model_gt, scene_gt, CameraMatrix, filepath_gt)
            
            gt_h = readdlm(string(filepath_gt,"/Results/data/h.csv"), ',', Float64, '\n', header=false)
            # Write the ground truth
            params = Dict("gt_η" => η_gt, "gt_β" => β_gt)
            write_json(string(filepath_gt,"/params"), params)

            ne_exp = 4
            model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β_gt, F, control, viscosity_type, sim_time, t_steps)
            for noiseLevel::Float64 in noiseLevelLst
                ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/Results/contour_data"))  
                ObsDataList = ObsDataList[1:(round(Int,sim_time)+1)]
                if noiseLevel == 0.0
                    # Read the gt data 
                    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
                    
                    for sides::Bool in sideList
                        filepathi = string(filepath,"/experiment_$(η_gt)_$(β_gt)/trials/noise_$(noiseLevel)")
                        conditions = Conditions(CameraMatrix=CameraMatrix,SIDES=sides,filepath=filepathi,ANIMATE=false)

                        θ = [ηStart, βStart]
                        stats = fit_model(model, scene, conditions, obsBorderPts, θ)
                                            
                        iterList = stats["iterList"]
                        costList = stats["costList"]
                        ηpList = stats["ηList"]
                        βpList = stats["βList"]
                        η = stats["η"]
                        β = stats["β"]

                        η_accuracy = abs(1 - abs(η-η_gt)/η_gt)*100
                        β_accuracy = abs(1 - abs(β-β_gt)/β_gt)*100
                        printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                        printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)
                        
                        model, scene = def_problem(r, h, ne_exp, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time_gt, t_steps_gt)
                        # simulate the model with the estimated parameters
                        est_μ_list, others = simulate(model, scene, conditions)

                        est_h = get_height(est_μ_list, h)

                        if maximum(ηpList) > η_gt+dev_η
                            ηStop = maximum(ηpList)*1.1
                        else
                            ηStop = η_gt+dev_η
                        end

                        if minimum(ηpList) < η_gt-dev_η
                            ηStart = minimum(ηpList)*0.9
                        else
                            ηStart = η_gt-dev_η
                        end

                        if maximum(βpList) > β_gt+dev_β
                            βStop = maximum(βpList)*1.1
                        else
                            βStop = β_gt+dev_β
                        end

                        if minimum(βpList) < β_gt-dev_β
                            βStart = minimum(βpList)*0.9
                        else
                            βStart = β_gt-dev_β
                        end

                        sample_space = 41
                        ηList = collect(range(ηStart, stop=ηStop, length=sample_space))
                        βList = collect(range(βStart, stop=βStop, length=sample_space))
                        CostMat = zeros(size(ηList,1),size(βList,1))

                        η_iter = 1:size(ηList,1)
                        β_iter = 1:size(βList,1)

                        for i::Int in η_iter
                            η = ηList[i]
                            for j::Int in β_iter
                                β = βList[j]
                                reset_model!(model)
                                model.η = [η]
                                scene.β = [β]
                                μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

                                # test the closest point function
                                d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 

                                CostMat[i,j] = sum(d_cp)
                            end
                        end

                        write_csv(string(filepathi,"/Results/data/cost_surface"), CostMat)
                        write_csv(string(filepathi,"/Results/data/cost_surface_η"), ηList)
                        write_csv(string(filepathi,"/Results/data/cost_surface_β"), βList)

                        write_csv(string(filepathi,"/Results/data/cost_surface_ηp"), ηpList)
                        write_csv(string(filepathi,"/Results/data/cost_surface_βp"), βpList)

                        write_csv(string(filepathi,"/Results/data/cost_steps"), costList)
                        write_csv(string(filepathi,"/Results/data/cost_steps_iter"), iterList)

                        write_csv(string(filepathi,"/Results/data/est_h"), est_h)

                        write_csv(string(filepathi,"/Results/data/border"), obsBorderPts[round(Int,sim_time/2)])

                        params = Dict("gt_η" => η_gt,
                        "gt_β" => β_gt,
                        "η" => η,
                        "β" => β,
                        "η_accuracy" => η_accuracy,
                        "β_accuracy" => β_accuracy)

                        write_json(string(filepathi,"/Results/data/stats"), params)

                        set_file(string(filepathi,"/Results/plots"))
                        Plots.plot(est_h, label="Estimated height", dpi=400)
                        Plots.plot!(gt_h, label="Ground truth height", dpi=400)
                        Plots.xlabel!("Time (s)")
                        Plots.ylabel!("Height")
                        Plots.savefig(string(filepathi,"/Results/plots/h_est.pdf"))

                        Plots.plot(abs.(est_h-gt_h), label="Height estimation error", dpi=400)
                        Plots.xlabel!("Time (s)")
                        Plots.ylabel!("Error")
                        Plots.savefig(string(filepathi,"/Results/plots/h_est_error.pdf"))

                        # Plot the cost function with iterations
                        Plots.plot(iterList, costList, label="Cost", marker=1, dpi=400, yscale=:log10)
                        Plots.xlabel!("Iterations")
                        Plots.ylabel!("Error")
                        Plots.savefig(string(filepathi,"/Results/plots/cost_steps.pdf"))
                        
                        # Plot the cost function surface
                        Plots.contour(ηList, βList, CostMat, color=:turbo, fill=false, levels=200, xlabel="η", ylabel="β", dpi=400)
                        Plots.plot!(ηpList, βpList, label="Estimations", marker=1)
                        Plots.plot!([η_gt], [β_gt], label="Ground truth", marker=:star5, markersize=5, markercolor=:green)
                        Plots.xlabel!("η")
                        Plots.ylabel!("β")
                        Plots.savefig(string(filepathi,"/Results/plots/cost_surface_iter.pdf"))
                    end
                else
                    n_samples = 10
                    η_pred = zeros(Float64, n_samples)
                    β_pred = zeros(Float64, n_samples)
                    costnList = Vector{AbstractVector}(undef, n_samples)
                    iternList = Vector{AbstractVector}(undef, n_samples)
                    filepathi = string(filepath,"/experiment_$(η_gt)_$(β_gt)/trials/noise_$(noiseLevel)")
                    est_h_list = zeros(Float64, round(Int,(steps+1)), n_samples)

                    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
                    filepathi = string(filepath,"/experiment_$(η_gt)_$(β_gt)/trials/noise_$(noiseLevel)")
                    for n::Int in 1:n_samples
                        obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

                        conditions = Conditions(CameraMatrix=CameraMatrix,filepath=filepathi)

                        θ = [ηStart, βStart]
                        stats = fit_model(model, scene, conditions, obsBorderPts, θ)
                                            
                        iterList = stats["iterList"]
                        costList = stats["costList"]
                        ηpList = stats["ηList"]
                        βpList = stats["βList"]
                        η = stats["η"]
                        β = stats["β"]

                        model, scene = def_problem(r, h, ne_exp, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time_gt, t_steps_gt)
                        
                        # simulate the model with the estimated parameters
                        est_μ_list, others = simulate(model, scene, conditions)

                        est_h_list[:,n] = get_height(est_μ_list, h)

                        η_accuracy = abs(1 - abs(η-η_gt)/η_gt)*100
                        β_accuracy = abs(1 - abs(β-β_gt)/β_gt)*100
                        printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                        printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

                        η_pred[n] = η
                        β_pred[n] = β
                        push!(costnList, costList)
                        push!(iternList, iterList)

                        params = Dict("gt_η" => η_gt,
                        "gt_β" => β_gt,
                        "η" => η,
                        "β" => β,
                        "η_accuracy" => η_accuracy,
                        "β_accuracy" => β_accuracy)

                        write_json(string(filepathi,"/Results/data/stats"), params)

                        write_csv(string(filepathi,"/Results/data/cost_steps/run_$n"), costList)
                        write_csv(string(filepathi,"/Results/data/cost_steps_iter/run_$n"), iterList)
                        write_csv(string(filepathi,"/Results/data/border"), obsBorderPts[round(Int,sim_time/2)])
                    end

                    write_csv(string(filepathi,"/Results/data/η_est"), η_pred)
                    write_csv(string(filepathi,"/Results/data/β_est"), β_pred)
                    write_csv(string(filepathi,"/Results/data/h_est"), est_h_list)

                    plot_covariance(η_pred, β_pred, string(filepathi,"/Results/plots/"))
                    
                    StatsPlots.errorline(est_h_list, label="Estimated height", dpi=400)
                    Plots.plot!(gt_h, label="Ground truth height", dpi=400)
                    Plots.xlabel!("Time (s)")
                    Plots.ylabel!("Height")
                    Plots.savefig(string(filepathi,"/Results/plots/h_est.pdf"))

                    set_file(string(filepathi,"/Results/plots"))
                    plot(x->pdf(pd, x),  size=(200,150), label="")
                    savefig(string(filepathi,"/Results/plots/obs_pdf_$(noiseLevel).pdf"))

                    Plots.plot(nSplinex[10], nSpliney[10], label="", size=(250,300), xaxis=false, yaxis=false) 
                    Plots.xlims!(1100,1350)
                    Plots.ylims!(420,1100)
                    Plots.savefig(string(filepathi,"/Results/plots/noise_contour_$(noiseLevel).pdf"))
                end
            end
        end
    end
end

# save the data to a file and post process

function plot_noise_covariance(ηLst, βLst, noiseLevel)
    
    for η in ηLst
        for β in βLst
            file_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test3/experiment_$(η)_$(β)/trials/")
            Plots.plot([],[],label="")
            for i::Float64 in noiseLevel
                if i != 0.0
                    filepath = string(file_path,"noise_$(i)/Results/data/")
                    η_list = readdlm(string(filepath,"η_est.csv"), ',', Float64, '\n', header=false)
                    β_list = readdlm(string(filepath,"β_est.csv"), ',', Float64, '\n', header=false)
    
                    mean_η = mean(η_list)
                    mean_β = mean(β_list)
                
                    cov_η = cov(η_list)
                    cov_β = cov(β_list)
                    cov_ηβ = cov(η_list, β_list)
                
                    cov_mat = [cov_η cov_ηβ; cov_ηβ cov_β]
                    mean_vec = [mean_η; mean_β]
                
                    # Plot the covariance matrix
                    StatsPlots.covellipse!(mean_vec, cov_mat, label="Noise variance $(i) pixels", aspect_ratio=0.25, xtickfont=font(10), ytickfont=font(10), legendfont=font(10), dpi=400)
                end
            end
            xlabel!("η")
            ylabel!("β")    
            Plots.savefig(string(file_path,"covariance.pdf"))
        end
    end
    filePath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test3")

    for η in ηLst
        Plots.plot([],[],label="")
        for β in βLst
            for i::Float64 in noiseLevel
                file_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test3/experiment_$(η)_$(β)")
                filepath = string(file_path,"/trials/noise_1.5/Results/data/")

                η_list = readdlm(string(filepath,"η_est.csv"), ',', Float64, '\n', header=false)
                β_list = readdlm(string(filepath,"β_est.csv"), ',', Float64, '\n', header=false)

                mean_η = mean(η_list)
                mean_β = mean(β_list)
            
                cov_η = cov(η_list)
                cov_β = cov(β_list)
                cov_ηβ = cov(η_list, β_list)
            
                cov_mat = [cov_η cov_ηβ; cov_ηβ cov_β]
                mean_vec = [mean_η; mean_β]
            
                # Plot the covariance matrix
                StatsPlots.covellipse(mean_vec, cov_mat, label="β = $(β)", aspect_ratio=0.25, xtickfont=font(10), ytickfont=font(10), legendfont=font(10), dpi=400)
                xlabel!("η")
                ylabel!("β")   
                Plots.savefig(string(file_path,"/covariance_vs_slip.pdf"))
            end
        end
    end

end

function plot_height_vs_slip(ηLst, βLst)
    βsz = length(βLst)
    file_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test3/")
    Plots.plot([],[],label="")

    for j::Float64 in ηLst
        for i::Float64 in βLst
            filepath = string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/h.csv")
            h = readdlm(filepath, ',', Float64, '\n', header=false)
            Plots.plot!(h, label=string("β = ", i), dpi=400)
        end
    end
    xlabel!("Time (s)")
    ylabel!("Height")
    Plots.savefig(string(file_path,"height_vs_slip.pdf"))
end

function plot_field_at_height(ηLst, βLst)
    h_Vector = Vector{AbstractArray}(undef, 4)
    β_Vector = Vector{Float64}(undef, 4)
    field_Vector = Vector{AbstractArray}(undef, 4)
    cont_vector = Vector{AbstractArray}(undef, 4)
    t_indexes = Vector{Int}(undef, 4)

    Plots.plot([],[],label="")
    file_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test3/")

    ηiter::Int = 1
    for j::Float64 in ηLst
        βiter::Int = 1
        for i::Float64 in βLst
            filepath = string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/h.csv")
            filepathJSON = string(file_path,"/experiment_$(j)_$(i)/ground_truth/params.json")
            data = read_json(filepathJSON)
            β = data["gt_β"]
            h = readdlm(filepath, ',', Float64, '\n', header=false)
            contours, tmp = read_csv(string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/contour_data"))
            fields, tmp = read_csv(string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/2D_surface_points/contour_data"))
            h_Vector[βiter] = h
            cont_vector[βiter] = contours
            β_Vector[βiter] = β
            field_Vector[βiter] = fields

            βiter += 1
        end
        ηiter += 1
    end
    # define a reference height
    h_ref = h_Vector[end][end]
    # find the closest time step to the reference height
    for i::Int in 1:4
        cost = 1e18
        h = h_Vector[i]
        sz_h = length(h)
        for id_t::Int in 1:sz_h
            if abs(h[id_t]-h_ref) < cost
                cost = abs(h[id_t]-h_ref)
                t_indexes[i] = id_t
            end
        end
        cont_t = cont_vector[i][t_indexes[i]]
        Plots.plot!(cont_t[1,:], cont_t[2,:], label=string("β = ",β_Vector[i]), dpi=400)
    end

    Plots.xlabel!("x")
    Plots.ylabel!("y")
    Plots.xlims!(720,740)
    Plots.ylims!(400,1100) 
    Plots.savefig(string(file_path,"field_at_height.pdf"))
    p = Vector{Plots.Plot}(undef, 8)

    for i::Int in 1:length(βLst)
        cont_t = cont_vector[i][t_indexes[i]]
        cont_tm1 = cont_vector[i][t_indexes[i]-1]
        field = cont_t - cont_tm1
        
        # normalize the field
        norm_field = get_norm(field)

        # select the elements to plot
        length_field = size(field, 2)
        elements = collect(range(start=1,stop=length_field, step=15))
        
        # println(minimum(field[1,elements]), " ", maximum(field[1,elements]))
        # println(minimum(field[2,elements]), " ", maximum(field[2,elements]))
        # println(minimum(abs.(field[1,elements])), " ", maximum(abs.(field[1,elements])))
        # println(minimum(abs.(field[2,elements])), " ", maximum(abs.(field[2,elements])))

        p[2*i-1] = Plots.plot(cont_t[1,:], cont_t[2,:], label=string("β = ",β_Vector[i]), dpi=800, aspect_ratio=0.35, xaxis=false, yaxis=false)
        arrow0!.(cont_t[1,elements], cont_t[2,elements], norm_field[1,elements], .0, field[1,elements], .0, as=0.1, lw=1, la=2)

        p[2*i] = Plots.plot(cont_t[1,:], cont_t[2,:], label=string("β = ",β_Vector[i]), dpi=800, aspect_ratio=0.35, xaxis=false, yaxis=false)
        arrow0!.(cont_t[1,elements], cont_t[2,elements], .0, norm_field[2,elements], .0, field[2,elements]; as=0.1, lw=1, la=2)

        Plots.plot(p[2*i-1], p[2*i], layout=grid(1,2, widths=(4/8,4/8)), size=(1600,560))
        Plots.ylims!(350,1300)
        Plots.xlims!(1000,1500)
        Plots.savefig(string(file_path,"field_at_height_quiver_$(i).pdf"))
    end

end

function arrow0!(x, y, u, v, color_u, color_v; as=0.07, lw=1, lc=:black, la=1, color=0.0)
    # resize the arrow
    nuv = sqrt(u^2 + v^2)
    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
    v5 = v4 - 2*(v4'*v2)*v2
    v4, v5 = as*nuv*v4, as*nuv*v5
    color = sqrt(color_u^2 + color_v^2)/800
    lc = color
    color_palette = :berlin
    Plots.plot!([x,x+u], [y,y+v], lw=lw, line_z=lc, la=la, label= "", c=color_palette)
    Plots.plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lw=lw, line_z=lc, la=la, label= "", c=color_palette)
    Plots.plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lw=lw, line_z=lc, la=la, label= "", c=color_palette)
end

function get_norm(x::AbstractMatrix)
    norm_vec = zeros(size(x,1))
    for i::Int in 1:size(x,1)
        norm_vec[i] = norm(x[i,:])
    end
    min_val = minimum(norm_vec)
    norm_x = x/min_val
    return norm_x.*800
end

function plot_data(ηLst, βLst, noiseLevelLst; n=0)
    file_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test3/")
    plot_colors = [:steelblue :indianred :seagreen :darkorange :mediumpurple :cadetblue :lightcoral :dimgray]
    for k::Float64 in noiseLevelLst
        if k == 0.0
            Plots.plot([],[],label="")
            for j::Float64 in ηLst
                for i::Float64 in βLst
                    filepath = string(file_path,"/experiment_$(j)_$(i)/trials/noise_$(k)/Results/data")
                    conv = readdlm(string(filepath,"/cost_steps.csv"), ',', Float64, '\n', header=false)
                    Plots.plot!(conv, label=string("β = ", i), dpi=400, yscale=:log10)
                end
            end
            xlabel!("Iterations")
            ylabel!("Cost")
            Plots.savefig(string(file_path,"cost_vs_slip.pdf"))
        else
            Plots.plot([],[],label="")
            n = 10  
            convMat = Vector{AbstractMatrix}(undef, n)
            color_iter = 1
            for j::Float64 in ηLst
                for i::Float64 in βLst
                    filepath = string(file_path,"/experiment_$(j)_$(i)/trials/noise_$(k)/Results/data")
                    for n::Int in 1:n
                        conv = readdlm(string(filepath,"/cost_steps/run_$(n).csv"), ',', Float64, '\n', header=false)
                        Plots.plot!(conv, label="", dpi=400, yscale=:log10, color=plot_colors[color_iter], lw=0.5)
                    end
                    Plots.plot!([],[],label=string("β = ", i),color=plot_colors[color_iter])
                    color_iter += 1
                end
            end
            xlabel!("Iterations")
            ylabel!("Cost")
            Plots.savefig(string(file_path,"cost_vs_slip_noise_$(k).pdf"))
        end
    end

    for k::Float64 in noiseLevelLst
        if k == 0.0
            Plots.plot([],[],label="")
            for j::Float64 in ηLst
                for i::Float64 in βLst
                    filepath_gt = string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/h.csv")
                    filepath = string(file_path,"/experiment_$(j)_$(i)/trials/noise_0.0/Results/data")

                    h_gt = readdlm(filepath_gt, ',', Float64, '\n', header=false)
                    h_est = readdlm(string(filepath,"/est_h.csv"), ',', Float64, '\n', header=false)
                    h_error = abs.(h_est - h_gt)

                    Plots.plot!(h_error, label=string("β = ", i), dpi=400)
                end
            end
            xlabel!("Time (s)")
            ylabel!("Height error")
            
            Plots.savefig(string(file_path,"height_error_vs_slip.pdf"))

            Plots.plot([],[],label="")
            for j::Float64 in ηLst
                for i::Float64 in βLst
                    filepath_gt = string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/h.csv")
                    filepath = string(file_path,"/experiment_$(j)_$(i)/trials/noise_0.0/Results/data")

                    h_est = readdlm(string(filepath,"/est_h.csv"), ',', Float64, '\n', header=false)

                    Plots.plot!(h_est, label=string("β = ", i), dpi=400)
                end
            end
            xlabel!("Time (s)")
            ylabel!("Height error")
            Plots.savefig(string(file_path,"height_est_vs_slip.pdf"))
        else
            Plots.plot([],[],label="")
            for j::Float64 in ηLst
                for i::Float64 in βLst
                    filepath_gt = string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/h.csv")
                    filepath = string(file_path,"/experiment_$(j)_$(i)/trials/noise_$(k)/Results/data")

                    h_gt = readdlm(filepath_gt, ',', Float64, '\n', header=false)
                    h_est = readdlm(string(filepath,"/h_est.csv"), ',', Float64, '\n', header=false)
                    h_error = abs.(h_est .- h_gt)

                    StatsPlots.errorline!(h_error, label=string("β = ", i), dpi=400, errorstyle=:stick)
                end
            end
            xlabel!("Time (s)")
            ylabel!("Height error")
            Plots.savefig(string(file_path,"height_error_vs_slip_noise.pdf"))

            Plots.plot([],[],label="")
            for j::Float64 in ηLst
                for i::Float64 in βLst
                    filepath_gt = string(file_path,"/experiment_$(j)_$(i)/ground_truth/Results/data/h.csv")
                    filepath = string(file_path,"/experiment_$(j)_$(i)/trials/noise_$(k)/Results/data")

                    h_est = readdlm(string(filepath,"/h_est.csv"), ',', Float64, '\n', header=false)

                    StatsPlots.errorline!(h_est, label=string("β = ", i), dpi=400, errorstyle=:stick)
                end
            end
            xlabel!("Time (s)")
            ylabel!("Height error")
            Plots.savefig(string(file_path,"height_est_vs_slip_noise_$(k).pdf"))
        end
    end
end

ηLst = [40.0]
βLst = [1 100.0 1000.0 1e4]
noiseLevelLst = [0.0 0.5 1.0]
# noiseLevelLst = [0.0]

main(βLst, noiseLevelLst, ηLst)
plot_height_vs_slip(ηLst, βLst)
plot_field_at_height(ηLst, βLst)
plot_noise_covariance(ηLst, βLst, noiseLevelLst)
plot_data(ηLst, βLst, noiseLevelLst, n=10)

