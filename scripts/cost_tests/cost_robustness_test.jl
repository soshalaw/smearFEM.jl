
using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using Plots

using Dates
function main()

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 3
    ndim = 3
    FunctionClass = "Q2"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 15
    tSteps = 45
    dateTime = Dates.now()
    sampleNo = 13
    dev = 0.3
    noiseProfileList = [1 2]

    noiseLevelLst = [0 0.25 0.5 1]
    YoungtstLst = [30 35 40]
    νtstLst = [0.25 0.3 0.35]
    βLst = [100 1000 10000]
    ControlList = ["force", "displacement"]

    # Derived Lame constants from Young's modulus and Poisson ratio
    lambdatstLst_ = YoungtstLst.*νtstLst./((νtstLst.+1).*(-2*νtstLst.+1))
    mutstLst_ = YoungtstLst./(2*(νtstLst.+1))

    # Round the values
    lambdatstLst = [round(i) for i in lambdatstLst_]
    mutstLst = [round(j) for j in mutstLst_]

    Standard = false # Test for Poisson ratio and Young's modulus
    Lame = true # Test for Lame constants

    folder = string("/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/cost_function_test/robusteness/exp_",dateTime, "/")
    
    iter = 1
    for noiseLevel in noiseLevelLst
        for noiseProfile in noiseProfileList
            for Youngtst in YoungtstLst
                for νtst in νtstLst
                    for βtst in βLst
                        dev_β = βtst*dev
                        βList = collect(range(βtst-dev_β, stop=βtst+dev_β, length=sampleNo))
                        for Control in ControlList
                            if Standard == true
                                ############################################################################################################################################
                                ## Testing with youngs modulus and poisson ratio
                                ############################################################################################################################################
                                filepath = string(folder,"/standard/")

                                for Youngtst in YoungtstLst
                                    for νtst in νtstLst
                                        println("For testing with Young's modulus and Poisson ratio...")

                                        filepathi = string(filepath,"experiment_",iter)
                                        write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi)
                                        
                                        ObsDataList, splinex, spliney = read_csv(filepathi, NOISE=true, nProfile=noiseProfile, nFactor=noiseLevel)

                                        dev_ν = νtst*dev
                                        dev_Young = Youngtst*dev

                                        YngList = collect(range(Youngtst-dev_Young, stop=Youngtst+dev_Young, length=sampleNo))
                                        νList = collect(range(νtst-dev_ν, stop=νtst+dev_ν, length=sampleNo))

                                        hcostListν = Float64[]
                                        cpCostListν = Float64[]
                                        hcostListYoung = Float64[]
                                        cpcostListYoung = Float64[]
                                        hCostListβ = Float64[]
                                        cpCostListβ = Float64[]

                                        for Young in YngList
                                            hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Young, νtst, ndim, FunctionClass, nDof, βtst, CameraMatrix, 
                                                                endTime, tSteps, Control, filepath=filepathi,NOISE=true, noiseProfile=noiseProfile,
                                                                noiseLevel=noiseLevel)                                       
                                            push!(hcostListν, sum(hcost)/length(hcost))
                                            push!(cpCostListν, sum(cpCost)/length(cpCost))
                                        end

                                        for ν in νList
                                            hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Youngtst, ν, ndim, FunctionClass, nDof, βtst, CameraMatrix, 
                                                                endTime, tSteps, Control, filepath=filepathi,NOISE=true, noiseProfile=noiseProfile,
                                                                noiseLevel=noiseLevel)
                                            push!(hcostListYoung, sum(hcost)/length(hcost))
                                            push!(cpcostListYoung, sum(cpCost)/length(cpCost))
                                        end

                                        for β in βList
                                            hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                endTime, tSteps, Control, filepath=filepathi,NOISE=true, noiseProfile=noiseProfile,
                                                                noiseLevel=noiseLevel)
                                            push!(hCostListβ, sum(hcost)/length(hcost))
                                            push!(cpCostListβ, sum(cpCost)/length(cpCost))
                                        end
 
                                        # plot the cost function with respect to the poisson ratio
                                        Plots.plot(νList, hcostListν, label="Height sample Cost", marker=1, dpi=400)
                                        Plots.vline!([νtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("Poisson Ratio (ν)")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_nu.png"))

                                        Plots.plot(νList, cpCostListν, label="Closest Point Cost", marker=1, dpi=400)
                                        Plots.vline!([νtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("Poisson Ratio (ν)")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_nu.png"))

                                        # plot the cost function with respect to the young's modulus
                                        Plots.plot(YngList, hcostListYoung, label="Height sample Cost", marker=1, dpi=400)
                                        Plots.vline!([Youngtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("Young's Modulus (E)")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_young.png"))

                                        Plots.plot(YngList, cpcostListYoung, label="Closest Point Cost", marker=1, dpi=400)
                                        Plots.vline!([Youngtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("Young's Modulus (E)")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_young.png"))

                                        # plot the cost function with respect to friction parameter
                                        Plots.plot(βList, hCostListβ, label="Height sample Cost", marker=1, dpi=400)
                                        Plots.vline!([βtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("β")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_β.png"))

                                        Plots.plot(βList, cpCostListβ, label="Closest Point Cost", marker=1, dpi=400)
                                        Plots.vline!([βtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("β")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_β.png"))

                                        params = Dict("Paramter type" => "Standard", "E" => Youngtst, "ν" => νtst, "β" => βtst, "Control" => Control,
                                                      "Noise Level" => noiseLevel, "Noise Profile" => noiseProfile)

                                        write_json(string(filepathi,"/Results/params"), params)
                                        iter += 1
                                    end
                                end
                            end
                            
                            if Lame == true
                            ############################################################################################################################################
                            ## Testing with lame constants
                            ############################################################################################################################################
                                filepathlame = string(folder,"/Lame/")
                                for lambdatst in lambdatstLst
                                    for mutst in mutstLst
                                        println("For testing with Lame constants...")

                                        dev_λ = lambdatst*dev
                                        dev_μ = mutst*dev

                                        filepathi = string(filepathlame,"experiment_",iter)
                                        simBorderPts, simBorderNodes, splinex, spliney = write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, βtst, CameraMatrix, endTime, tSteps, Control,filepathi, mode="lame")
                                        ObsDataList, splinexObs, splineyObs = readData(filepathi, NOISE=true, nProfile=noiseProfile, nFactor=noiseLevel)
                                        
                                        ObsData = [ObsDataList, splinexObs, splineyObs]

                                        println(size(splineyObs))
                                        # animate_fields(filepath = string(filepathi,"/Results"), p=splinex, q=spliney)
                                        animate_fields(filepath = string(filepathi,"/Results"), BorderNodes2D=simBorderPts, pObs=splinexObs, qObs=splineyObs)
                                        lambdaList = collect(range(lambdatst-dev_λ, stop=lambdatst+dev_λ, length=sampleNo))
                                        muList = collect(range(mutst-dev_μ, stop=mutst+dev_μ, length=sampleNo))

                                        hCostListμ = Float64[]
                                        cpCostListμ = Float64[]
                                        hCostListλ = Float64[]
                                        cpCostListλ = Float64[]
                                        hCostListβ = Float64[]
                                        cpCostListβ = Float64[]

                                        for λ in lambdaList
                                            hcost, cpCost = compare(x0, x1, y0, y1, z0, z1, ne, λ, mutst, ndim, FunctionClass, nDof, βtst, CameraMatrix, 
                                                                endTime, tSteps, Control, "lame", ObsData)    
                                            push!(hCostListλ, sum(hcost)/length(hcost))
                                            push!(cpCostListλ, sum(cpCost)/length(cpCost))
                                        end

                                        for μ in muList
                                            hcost, cpCost = compare(x0, x1, y0, y1, z0, z1, ne, lambdatst, μ, ndim, FunctionClass, nDof, βtst, CameraMatrix, 
                                                                    endTime, tSteps, Control, "lame", ObsData)
                                            push!(hCostListμ, sum(hcost)/length(hcost))
                                            push!(cpCostListμ, sum(cpCost)/length(cpCost))
                                        end

                                        for β in βList 
                                            hcost, cpCost = compare(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                    endTime, tSteps, Control, "lame", ObsData)
                                            push!(hCostListβ, sum(hcost)/length(hcost))
                                            push!(cpCostListβ, sum(cpCost)/length(cpCost))
                                        end

                                        # plot the cost function with respect to λ
                                        Plots.plot(lambdaList, hCostListλ, label="Height sample Cost", marker=1, dpi=400)
                                        Plots.vline!([lambdatst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("λ")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_lbd.png"))

                                        Plots.plot(lambdaList, cpCostListλ, label="Closest Point Cost", marker=1, dpi=400)
                                        Plots.vline!([lambdatst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("λ")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_lbd.png"))

                                        # plot the cost function with respect to μ
                                        Plots.plot(muList, hCostListμ, label="Height sample Cost", marker=1, dpi=400)
                                        Plots.vline!([mutst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("μ")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_mu.png"))

                                        Plots.plot(muList, cpCostListμ, label="Closest Point Cost", marker=1, dpi=400)
                                        Plots.vline!([mutst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("μ")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_mu.png"))

                                        # plot the cost function with respect to friction parameter
                                        Plots.plot(βList, hCostListβ, label="Height sample Cost", marker=1, dpi=400)
                                        Plots.vline!([βtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("β")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_β.png"))

                                        Plots.plot(βList, cpCostListβ, label="Closest Point Cost", marker=1, dpi=400)
                                        Plots.vline!([βtst],linestyle=:dash,linecolor=:grey)
                                        Plots.xlabel!("β")
                                        Plots.ylabel!("Mean Error")
                                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_β.png"))

                                        params = Dict("Paramter type" => "Lame", "λ" => lambdatst, "μ" => mutst, "β" => βtst, "Control" => Control,
                                                      "Noise Level" => noiseLevel, "Noise Profile" => noiseProfile)

                                        write_json(string(filepathi,"/Results/params"), params)
                                        iter += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

main()