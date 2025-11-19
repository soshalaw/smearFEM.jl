using smearFEM
using LinearAlgebra
using Plots
using ProgressMeter
using DelimitedFiles

function test_linear_elas()

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 8
    ndim = 3
    FunctionClass = "Q1"
    nDof = ndim  # number of degree of freedom per node
    Young = 30
    ν = 0.3

    β_lst = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5]
    ne_lst = [2, 4, 8, 16, 32, 64]

    condition_number_lst_β = Float64[]
    condition_number_lst_ne = Float64[]

    for β in β_lst
        push!(condition_number_lst_β, condition_number(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β))
    end

    β = 100
    Plots.plot(β_lst, condition_number_lst_β, yscale=:log10, xscale=:log10,label="Condition number vs β", xlabel="β", ylabel="Condition number", title="Condition number vs β")

    for ne in ne_lst
        push!(condition_number_lst_ne, condition_number(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β))
    end

end

function test_stokes()

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 5
    ndim = 3
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    ν = 1

    β_lst = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5]

    condition_number_lst_β = Float64[]
    condition_number_lst_ne = Float64[]

    stepsβ = (length(β_lst))

    β = 100

    j = 1
    pr = Progress(stepsβ; desc="Convergance analysis with β ...", showspeed=true)
    for β in β_lst
        cond_num_β = condition_number_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β)
        push!(condition_number_lst_β, cond_num_β)
        next!(pr, showvalues = [(:iterations,j)])
        try
            yield()
        catch
        end
        j += 1
    end
    
    Plots.plot(β_lst, condition_number_lst_β, yscale=:log10, xscale=:log10,label="Condition number vs β", xlabel="β", ylabel="Condition number", title="Condition number vs β")
    Plots.savefig("/home/soshala/condition_number_beta.png")
    
    ne_lst = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
    stepsne = (length(ne_lst))

    i = 1
    pr = Progress(stepsne; desc="Convergance analysis with mesh size ...", showspeed=true)
    for ne in ne_lst
        cond_num_ne = condition_number_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β)
        i += 1
        push!(condition_number_lst_ne, cond_num_ne)
        next!(pr, showvalues = [(:iterations,i)])
        try
            yield()
        catch
        end
    end

    Plots.plot(ne_lst, condition_number_lst_ne, yscale=:log10, xscale=:log10,label="Condition number vs Elements", xlabel="Elements", ylabel="Condition number", title="Condition number vs Elements")
    Plots.savefig("/home/soshala/condition_number_stokes.png")


    return condition_number_lst_ne
end

function condition_number_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β)    
    NodeList_u, IEN_u, ID_u, IEN_u_top, IEN_u_btm, BorderNodesList_u = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_u)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList_u, x0, x1, y0, y1)
    # q_tp, q_side, q_btm, C_uc = set_boundary_cond_stokes(NodeList_u, ne, ndim, FunctionClass_u, nDof_u)
    
    NodeList_p, IEN_p, ID_p, IEN_p_top, IEN_p_btm, BorderNodesList_p = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_p)  # generate the mesh grid
    # NodeListCylinderp = inflate_cylinder(NodeList_p, x0, x1, y0, y1)
    
    mdl = def_model("stokes", ne=ne, NodeList=NodeListCylinder, IEN=IEN_u, IEN_top=IEN_u_top, IEN_btm=IEN_u_btm, ndim=ndim, nDof=nDof_u, 
                        FunctionClass=FunctionClass_u, ID=ID_u, IEN_2=IEN_p, IEN_2_top=IEN_p_top, IEN_2_btm=IEN_p_btm, 
                            ndim_2=ndim, nDof_2=nDof_p, FunctionClass_2=FunctionClass_p ) # define the model

    A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
    B = assemble_system_B(mdl)                       # assemble the stiffness matrix
    b = apply_boundary_conditions_stokes(mdl)        # apply the neumann boundary conditions

    A = A_bar + β*b

    K = [A B; B' zeros(size(B,2),size(B,2))]     # assemble the system of equations

    cNum = cond(K,1)

    return cNum
end

function condition_number_linear_elas(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β)
    mode = "standard"

    cMat = get_cMat(mode, Young, ν)

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere
    # q_tp, q_btm, C_uc = setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof)

    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeListCylinder, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                        FunctionClass=FunctionClass, Young=Float64(Young), ν=ν, cMat=cMat)
        
    K = assemble_system(mdl)                   # assemble the stiffness matrix
    b = apply_boundary_conditions(mdl)         # apply the neumann boundary conditions

    K_bar = K + β*b

    return cond(K_bar,1)
end

res = test_stokes()

# cStr = string(counter,pad=3)

filename = "/home/soshala/data"

open(string(filename,".csv"), "w") do io
    writedlm(io, res, ',')
end