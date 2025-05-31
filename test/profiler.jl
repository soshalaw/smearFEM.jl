# using ProfileView
using smearFEM
using BenchmarkTools
using DelimitedFiles

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    neLst = [2 4 6 8 10 12 14 16 18 20]
    # neLst = [2]
    ndim = 3
    FunctionClass = "Q2"
    FunctionClass_ = "Q1"
    nDof = ndim  # number of degree of freedom per node
    nDof_ = 1
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 30
    tSteps = 45
    
    Control = "force" # "force" or "displacement"

    mode = "standard" # "standard" or "lame"

    writeData = false
    Young = 12.0
    ν = 12.0
    μ_tp = -0.1
    μ_btm = 0
    μ_side = 0
    filepathi = string("/home/soshala/tmp/")

    bsparseT = Float64[]
    bdenseT = Float64[]
    bsparseCGT = Float64[]
    bdenseCGT = Float64[]
    
    # @profview write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)
    # @benchmark write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)
    
    # write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)

    # @profview test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filepathi, mode = mode)
    # @benchmark test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filepathi, mode = mode)

    # @pro @profview test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filename, mode = mode)
    
    for ne in neLst
        bresSparse = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm) setup=(ne=$ne)
        bresDense = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, DENSE=true) setup=(ne=$ne)

        bresSparseCG = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, CG=true) setup=(ne=$ne)
        bresDenseCG = @benchmark simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, DENSE=true, CG=true) setup=(ne=$ne)

        push!(bsparseT, median(bresSparse).time/1e9)
        push!(bdenseT, median(bresDense).time/1e9)
        push!(bsparseCGT, median(bresSparseCG).time/1e9)
        push!(bdenseCGT, median(bresDenseCG).time/1e9)

        println("ne = ", ne)
        println("Sparse time = ", median(bresSparse).time/1e9)
        println("Dense time = ", median(bresDense).time/1e9)
        println("Sparse CG time = ", median(bresSparseCG).time/1e9)
        println("Dense CG time = ", median(bresDenseCG).time/1e9)
    end 
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
    