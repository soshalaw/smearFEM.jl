function plot_noise_covariance(ηLst, βLst, noiseLevel, file_path)
    
    for η in ηLst
        for β in βLst
            file_path = string(file_path,"/experiment_$(η)_$(β)/trials/")
            set_plot(22)
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
    filePath = resolve_data_path("sim_experiments/cost_function_test/optimization/Stokes/force/test3")

    for η in ηLst
        set_plot(22)
        for β in βLst
            for i::Float64 in noiseLevel
                file_path = resolve_data_path("sim_experiments/cost_function_test/optimization/Stokes/force/test3/experiment_$(η)_$(β)")
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

function plot_height_vs_slip(ηLst, βLst, file_path)
    βsz = length(βLst)
    set_plot(22)
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

function plot_field_at_height(ηLst, βLst, file_path)
    scene_size = length(ηLst)* length(βLst)
    h_Vector = Vector{AbstractArray}(undef, scene_size)
    β_Vector = Vector{Float64}(undef, scene_size)
    field_Vector = Vector{AbstractArray}(undef, scene_size)
    cont_vector = Vector{AbstractArray}(undef, scene_size)
    t_indexes = Vector{Int}(undef, scene_size)

    set_plot(22)
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
    println("Reference height: ", h_Vector)
    h_ref = h_Vector[end][1]
    # find the closest time step to the reference height
    for i::Int in 1:scene_size
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
        cont_tp1 = cont_vector[i][t_indexes[i]+1]
        field = cont_tp1-cont_t
        
        # normalize the field
        norm_field = get_norm(field)

        # select the elements to plot
        length_field = size(field, 2)
        elements = collect(range(start=1,stop=length_field, step=15))

        # Plots.plot(cont_t[1,:], cont_t[2,:], label=string("β = ",β_Vector[i]), dpi=800, aspect_ratio=0.35, xaxis=false, yaxis=false)
        set_plot(22)
        Plots.plot(cont_t[1,:], cont_t[2,:], label="", dpi=800, aspect_ratio=0.35, xaxis=false, yaxis=false)
        arrow0!.(cont_t[1,elements], cont_t[2,elements], norm_field[1,elements], norm_field[2,elements], .0, .0, as=0.1, lw=1, la=2)
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
    # color = sqrt(color_u^2 + color_v^2)/800
    # lc = color
    # color_palette = :berlin
    set_plot(22)
    # Plots.plot!([x,x+u], [y,y+v], lw=lw, line_z=lc, la=la, label= "", c=color_palette)
    # Plots.plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lw=lw, line_z=lc, la=la, label= "", c=color_palette)
    # Plots.plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lw=lw, line_z=lc, la=la, label= "", c=color_palette)
    Plots.plot!([x,x+u], [y,y+v], lw=lw, lc=:black, la=la, label= "")
    Plots.plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lw=lw, lc=:black, la=la, label= "")
    Plots.plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lw=lw, lc=:black, la=la, label= "")
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

function plot_data(ηLst, βLst, noiseLevelLst, file_path; n=0)

    plot_colors = [:steelblue :indianred :seagreen :darkorange :mediumpurple :cadetblue :lightcoral :dimgray]
    for k::Float64 in noiseLevelLst
        if k == 0.0
            set_plot(22)
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
            set_plot(22)
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
            set_plot(22)
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

            set_plot(22)
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
            set_plot(22)
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

            set_plot(22)
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