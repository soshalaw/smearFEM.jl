
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
    ne = 4
    ndim = 3
    FunctionClass = "Q2"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    endTime = 15
    tSteps = 45
    dateTime = Dates.now()
    sampleNo = 13
    
    YoungtstLst = [30]
    νtstLst = [0.3]
    βLst = [100]
    ControlList = ["displacement"]

    # Derived Lame constants from Young's modulus and Poisson ratio
    lambdatstLst_ = YoungtstLst.*νtstLst./((νtstLst.+1).*(-2*νtstLst.+1))
    mutstLst_ = YoungtstLst./(2*(νtstLst.+1))

    # Round the values
    lambdatstLst = [round(i) for i in lambdatstLst_]
    mutstLst = [round(j) for j in mutstLst_]

    Standard = false # Test for Poisson ratio and Young's modulus
    Lame = true # Test for Lame constants

    homeDir = homedir()
    filepath = string(homeDir,"/SMEAR-PhD/SMEAR/Data/sim_experiments/cost_function_test/robusteness/",Date(dateTime),"/",Time(dateTime),"/")
    for βtst in βLst
        dev_β = βtst*0.3
        βList = collect(range(β-dev_β, stop=β+dev_β, length=sampleNo))
        iter = 1
        for Control in ControlList
            if Standard == true
                ############################################################################################################################################
                ## Testing with youngs modulus and poisson ratio
                ############################################################################################################################################
                for Youngtst in YoungtstLst
                    for νtst in νtstLst
                        println("For testing with Young's modulus and Poisson ratio...")

                        filepathi = string(filepath,"experiment_",iter)
                        write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi)
                        
                        dev_ν = νtst*0.3
                        dev_Young = Youngtst*0.3

                        YngList = collect(range(Youngtst-dev_Young, stop=Youngtst+dev_Young, length=sampleNo))
                        νList = collect(range(νtst-dev_ν, stop=νtst+dev_ν, length=sampleNo))

                        hcostListν = Float64[]
                        cpCostListν = Float64[]
                        hcostListYoung = Float64[]
                        cpcostListYoung = Float64[]
                        hcostListβ = Float64[]
                        cpCostListβ = Float64[]

                        costMatcp = zeros(length(YngList),length(νList))
                        costMath = zeros(length(YngList),length(νList))

                        hcost = zeros(Float64,tSteps)
                        cpCost = zeros(Float64,tSteps)

                        i = 1
                        for Young in YngList
                            j = 1
                            for ν in νList
                                for β in βList
                                    hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi)
                                    if Young == Youngtst
                                        push!(hcostListν, sum(hcost)/length(hcost))
                                        push!(cpCostListν, sum(cpCost)/length(cpCost))
                                    end
                                    if ν == νtst
                                        push!(hcostListYoung, sum(hcost)/length(hcost))
                                        push!(cpcostListYoung, sum(cpCost)/length(cpCost))
                                    end
                                    if β == βtst
                                        push!(hcostListβ, sum(hcost)/length(hcost))
                                        push!(cpCostListβ, sum(cpCost)/length(cpCost))
                                    end
                                end
                                costMath[i,j] = sum(hcost)/length(hcost)
                                costMatcp[i,j] = sum(cpCost)/length(cpCost)
                                j += 1
                            end
                            i += 1
                        end
                        # plot the cost function with respect to the poisson ratio
                        Plots.plot(νList, hcostListν, label="Height sample Cost", marker=1, dpi=400)
                        Plots.xlabel!("Poisson Ratio (ν)")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_nu.png"))

                        Plots.plot(νList, cpCostListν, label="Closest Point Cost", marker=1, dpi=400)
                        Plots.xlabel!("Poisson Ratio (ν)")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_nu.png"))

                        # plot the cost function with respect to the young's modulus
                        Plots.plot(YngList, hcostListYoung, label="Height sample Cost", marker=1, dpi=400)
                        Plots.xlabel!("Young's Modulus (E)")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_young.png"))

                        Plots.plot(YngList, cpcostListYoung, label="Closest Point Cost", marker=1, dpi=400)
                        Plots.xlabel!("Young's Modulus (E)")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_young.png"))

                        # plot the cost function with respect to friction parameter
                        Plots.plot(βList, hcostListβ, label="Height sample Cost", marker=1, dpi=400)
                        Plots.xlabel!("β")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_β.png"))

                        Plots.plot(βList, cpCostListβ, label="Closest Point Cost", marker=1, dpi=400)
                        Plots.xlabel!("β")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_β.png"))

                        # contour plot of the cost function with respect to the young's modulus and poisson ratio
                        Plots.contourf(νList, YngList, costMath, c=:viridis, ylabel="Young's Modulus (E)", xlabel="Poisson Ratio (ν)", title="Height Sample Cost", levels=40)
                        Plots.savefig(string(filepathi,"/Results/cost/contour_height.png"))

                        Plots.contourf(νList, YngList, costMatcp, c=:viridis, ylabel="Young's Modulus (E)", xlabel="Poisson Ratio (ν)", title="Closest Point Cost", levels=40)
                        Plots.savefig(string(filepathi,"/Results/cost/contour_cp.png"))

                        params = Dict("Paramter type" => "standard", "Young" => Youngtst, "ν" => νtst, "β" => βtst)

                        write_json(string(filepathi,"/Results/params"), params)
                        iter += 1
                    end
                end
            end
            
            if Lame == true
            ############################################################################################################################################
            ## Testing with lame constants
            ############################################################################################################################################
                filepathlame = string("/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/cost_function_test/exp_",dateTime,"/Lame/")
                for lambdatst in lambdatstLst
                    for mutst in mutstLst
                        println("For testing with Lame constants...")

                        dev_λ = lambdatst*0.3
                        dev_μ = mutst*0.3

                        filepathi = string(filepathlame,"experiment_",iter)
                        write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode="lame")

                        lambdaList = collect(range(lambdatst-dev_λ, stop=lambdatst+dev_λ, length=sampleNo))
                        muList = collect(range(mutst-dev_μ, stop=mutst+dev_μ, length=sampleNo))

                        hcostListμ = Float64[]
                        cpCostListμ = Float64[]
                        hcostListλ = Float64[]
                        cpcostListλ = Float64[]
                        hcostListβ = Float64[]
                        cpCostListβ = Float64[]

                        costMatcp = zeros(length(lambdaList),length(muList))
                        costMath = zeros(length(lambdaList),length(muList))

                        i = 1
                        for λ in lambdaList
                            j = 1
                            for μ in muList
                                hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi, mode="lame")
                                if λ == lambdatst
                                    push!(hcostListμ, sum(hcost)/length(hcost))
                                    push!(cpCostListμ, sum(cpCost)/length(cpCost))
                                end
                                if μ == mutst
                                    push!(hcostListλ, sum(hcost)/length(hcost))
                                    push!(cpcostListλ, sum(cpCost)/length(cpCost))
                                end
                                costMath[i,j] = sum(hcost)/length(hcost)
                                costMatcp[i,j] = sum(cpCost)/length(cpCost)
                                j += 1
                            end
                            i += 1
                        end

                        for β in βList
                            hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi, mode="lame")
                            push!(hcostListβ, sum(hcost)/length(hcost))
                            push!(cpCostListβ, sum(cpCost)/length(cpCost))
                        end

                        # plot the cost function with respect to λ
                        Plots.plot(lambdaList, hcostListλ, label="Height sample Cost", marker=1, dpi=400)
                        Plots.xlabel!("λ")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_lbd.png"))

                        Plots.plot(lambdaList, cpcostListλ, label="Closest Point Cost", marker=1, dpi=400)
                        Plots.xlabel!("λ")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_lbd.png"))

                        # plot the cost function with respect to μ
                        Plots.plot(muList, hcostListμ, label="Height sample Cost", marker=1, dpi=400)
                        Plots.xlabel!("μ")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_mu.png"))

                        Plots.plot(muList, cpCostListμ, label="Closest Point Cost", marker=1, dpi=400)
                        Plots.xlabel!("μ")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_mu.png"))

                        # plot the cost function with respect to friction parameter
                        Plots.plot(βList, hcostListβ, label="Height sample Cost", marker=1, dpi=400)
                        Plots.xlabel!("β")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_height_β.png"))

                        Plots.plot(βList, cpCostListβ, label="Closest Point Cost", marker=1, dpi=400)
                        Plots.xlabel!("β")
                        Plots.ylabel!("Mean Error")
                        Plots.savefig(string(filepathi,"/Results/cost/cost_cp_β.png"))

                        # contour plot of the cost function with respect to λ and μ
                        Plots.contourf(lambdaList, muList, costMath, levels=40, c=:viridis, xlabel="λ", ylabel="μ", title="Height Sample Cost")
                        Plots.savefig(string(filepathi,"/Results/cost/contour_height.png"))

                        Plots.contourf(lambdaList, muList, costMatcp, levels=40, c=:viridis, xlabel="λ", ylabel="μ", title="Closest Point Cost")
                        Plots.savefig(string(filepathi,"/Results/cost/contour_cp.png"))

                        params = Dict("Paramter type" => "Lame", "Young" => lambdatst, "ν" => mutst, "β" => βtst)

                        write_json(string(filepathi,"/Results/params"), params)
                        iter += 1
                    end
                end
            end
        end
    end
end

main()