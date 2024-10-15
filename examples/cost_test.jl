
using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using Plots

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
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 15
    tSteps = 45
    
    Control = "displacement" # "force" or "displacement"

    filepath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/standard/"
    filepathlame = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/Lame/"

    ############################################################################################################################################
    ## Testing with youngs modulus and poisson ratio
    ############################################################################################################################################
    YoungtstLst = [40]
    νtstLst = [0.35]

    iter = 1
    for Youngtst in YoungtstLst
        for νtst in νtstLst

            filepathi = string(filepath,"experiment_",iter)
            write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi)
            
            println("For testing with Young's modulus and Poisson ratio...")
            YngList = (Youngtst-1):0.1:(Youngtst+1)
            νList = (νtst-0.05):0.005:(νtst+0.05)

            hcostListν = []
            cpCostListν = []
            hcostListYoung = []
            cpcostListYoung = []

            costMatcp = zeros(length(YngList),length(νList))
            costMath = zeros(length(YngList),length(νList))

            println("For plotting the contours...")
            i = 1
            for Young in YngList
                j = 1
                for ν in νList
                    hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi)
                    if Young == Youngtst
                        push!(hcostListν, sum(hcost)/length(hcost))
                        push!(cpCostListν, sum(cpCost)/length(cpCost))
                    end
                    if ν == νtst
                        push!(hcostListYoung, sum(hcost)/length(hcost))
                        push!(cpcostListYoung, sum(cpCost)/length(cpCost))
                    end
                    costMath[i,j] = sum(hcost)/length(hcost)
                    costMatcp[i,j] = sum(cpCost)/length(cpCost)
                    j += 1
                end
                i += 1
            end

            # plot the cost function with respect to the poisson ratio
            Plots.plot(νList, hcostListν, label="Height sample Cost", marker=1, dpi=400)
            Plots.xlabel!("Poisson Ratio")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost/cost_height_nu",Int8(νtst*100),"young",Youngtst,".png"))

            Plots.plot(νList, cpCostListν, label="Closest Point Cost", marker=1, dpi=400)
            Plots.xlabel!("Poisson Ratio")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost/cost_cp_nu",Int8(νtst*100),"young",Youngtst,".png"))

            # plot the cost function with respect to the young's modulus
            Plots.plot(YngList, hcostListYoung, label="Height sample Cost", marker=1, dpi=400)
            Plots.xlabel!("Young's Modulus")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost/cost_height_young",Int8(νtst*100),"young",Youngtst,".png"))

            Plots.plot(YngList, cpcostListYoung, label="Closest Point Cost", marker=1, dpi=400)
            Plots.xlabel!("Young's Modulus")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost/cost_cp_young",Int8(νtst*100),"young",Youngtst,".png"))

            # contour plot of the cost function with respect to the young's modulus and poisson ratio
            Plots.contourf(νList, YngList, costMath, c=:viridis, ylabel="Young's Modulus (E)", xlabel="Poisson Ratio (ν)", title="Height Sample Cost", levels=40)
            Plots.savefig(string(filepathi,"/Results/cost/contour_height_young",Int8(νtst*100),"young",Youngtst,".png"))

            Plots.contourf(νList, YngList, costMatcp, c=:viridis, ylabel="Young's Modulus (E)", xlabel="Poisson Ratio (ν)", title="Closest Point Cost", levels=40)
            Plots.savefig(string(filepathi,"/Results/cost/contour_cp_young",Int8(νtst*100),"young",Youngtst,".png"))

            iter += 1
        end
    end

    ############################################################################################################################################
    ## Testing with lame constants
    ############################################################################################################################################
    lambdatstLst = [34.0]
    mutstLst = [14.0]

    iter = 1
    for lambdatst in lambdatstLst
        for mutst in mutstLst
            println("For testing with Lame constants...")

            filepathi = string(filepathlame,"experiment_",iter)
            write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode="lame")

            lambdaList = (lambdatst-1):0.1:(lambdatst+1)
            muList = (mutst-1):0.1:(mutst+1)

            hcostListμ = []
            cpCostListμ = []
            hcostListλ = []
            cpcostListλ = []

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

            # plot the cost function with respect to the poisson ratio
            Plots.plot(lambdaList, hcostListλ, label="Height sample Cost", marker=1, dpi=400)
            Plots.xlabel!("λ")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost_lame/cost_height_lbd",lambdatst,"mu",mutst,".png"))

            Plots.plot(lambdaList, cpcostListλ, label="Closest Point Cost", marker=1, dpi=400)
            Plots.xlabel!("λ")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost_lame/cost_cp_lbd",lambdatst,"mu",mutst,".png"))

            # plot the cost function with respect to the young's modulus
            Plots.plot(muList, hcostListμ, label="Height sample Cost", marker=1, dpi=400)
            Plots.xlabel!("μ")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost_lame/cost_height_mu",lambdatst,"mu",mutst,".png"))

            Plots.plot(muList, cpCostListμ, label="Closest Point Cost", marker=1, dpi=400)
            Plots.xlabel!("μ")
            Plots.ylabel!("Mean Error")
            Plots.savefig(string(filepathi,"/Results/cost_lame/cost_cp_mu",lambdatst,"mu",mutst,".png"))

            # contour plot of the cost function with respect to the young's modulus and poisson ratio
            Plots.contourf(lambdaList, muList, costMath, levels=40, c=:viridis, xlabel="λ", ylabel="μ", title="Height Sample Cost")
            Plots.savefig(string(filepathi,"/Results/cost_lame/contour_height_lbd",lambdatst,"mu",mutst,".png"))

            Plots.contourf(lambdaList, muList, costMatcp, levels=40, c=:viridis, xlabel="λ", ylabel="μ", title="Closest Point Cost")
            Plots.savefig(string(filepathi,"/Results/cost_lame/contour_cp_lbd",lambdatst,"mu",mutst,".png"))

            iter += 1
        end
    end
end

main()