# using ProfileView
using smearFEM
using BenchmarkTools
using DelimitedFiles

    # test case 
    r = 0.5
    h = 1
    ne = 4
    ndim = 3
    FunctionClass_u = "Q2"
    nDof_u = ndim  # number of degree of freedom per node
    FunctionClass_p = "Q1"
    nDof_p = 1  # number of degree of freedom per node
  
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  
    sim_time = 15.0
    steps = 15.0
    t_steps = sim_time/steps
  
    β = 100.0
    η = 40.0
    F = 1.0
    control = "force" # "force" or "displacement"
    viscosity_type = "constant" # "constant" or "linear"
    filepath = string("/home/soshala/tmp/")

    bsparseT = Float64[]
    bdenseT = Float64[]
    bsparseCGT = Float64[]
    bdenseCGT = Float64[]
    
    model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time, t_steps)

    @profview write_sim_data(model, scene, CameraMatrix, filepath)
    # @benchmark write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)
    
    # write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)
    # @profview test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filepathi, mode = mode)
    # @benchmark test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filepathi, mode = mode)

    # @pro @profview test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filename, mode = mode)
    
    # for ne in neLst
    #     bresSparse = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm) setup=(ne=$ne)
    #     bresDense = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, DENSE=true) setup=(ne=$ne)

    #     bresSparseCG = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, CG=true) setup=(ne=$ne)
    #     bresDenseCG = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, DENSE=true, CG=true) setup=(ne=$ne)

    #     push!(bsparseT, median(bresSparse).time/1e9)
    #     push!(bdenseT, median(bresDense).time/1e9)
    #     push!(bsparseCGT, median(bresSparseCG).time/1e9)
    #     push!(bdenseCGT, median(bresDenseCG).time/1e9)

    #     println("ne = ", ne)
    #     println("Sparse time = ", median(bresSparse).time/1e9)
    #     println("Dense time = ", median(bresDense).time/1e9)
    #     println("Sparse CG time = ", median(bresSparseCG).time/1e9)
    #     println("Dense CG time = ", median(bresDenseCG).time/1e9)
    # end 
    # @profview simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, DENSE=dense)

    # bench_results = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, DENSE=dense)

    # @profview simulate_single_tstep_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass, FunctionClass_, nDof, nDof_, β, μ_tp, μ_btm, μ_side)

    # bench_results = @benchmark simulate_single_tstep_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass, FunctionClass_, nDof, nDof_, β, μ_tp, μ_tp, μ_side)

    # rm(filepathi,recursive=true)


    # filenameSparse = "/home/soshala/dense_times"
    # filenameDense = "/home/soshala/sparse_times"
    # filenameSparseCG = "/home/soshala/dense_times_cg"
    # filenameDenseCG = "/home/soshala/sparse_times_cg"
    
    # open(string(filenameSparse,".csv"), "w") do io
    #     writedlm(io, bdenseT, ',')
    # end

    # open(string(filenameDense,".csv"), "w") do io
    #     writedlm(io, bsparseT, ',')
    # end

    # open(string(filenameSparseCG,".csv"), "w") do io
    #     writedlm(io, bdenseCGT, ',')
    # end

    # open(string(filenameDenseCG,".csv"), "w") do io
    #     writedlm(io, bsparseCGT, ',')
    # end
    
    rm(filepathi,recursive=true)