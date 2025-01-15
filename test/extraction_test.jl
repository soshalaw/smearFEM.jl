using smearFEM
using LinearAlgebra
using Test

filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"
# CPointList, W, C, IEN, IEN_top, IEN_btm = read_h5(string(filePath,"/cylinder.h5"))
CPointList, W, C, IEN = read_h5(string(filePath,"/cylinder.h5"))
NodeList = zeros(Float64, 3, size(IEN,2)*size(IEN,1))
IEN_list = zeros(Int64, size(IEN))

ndim = size(CPointList,1)
ne = Int64(round((size(IEN,2))^(1/ndim)))

# dummy solution
q = ones(3,size(CPointList,2))

ID = zeros(Int64, ndim, size(CPointList,2))
cpiter = 1:size(CPointList,2)
for m in cpiter
    for l in 1:ndim
        ID[l,m] = ndim*(m-1) + l
    end
end 

# construct the parent element
lPoints = zeros(27,3)
lPoints[1,:] = [-1 , -1, -1]
lPoints[2,:] = [1 , -1, -1]
lPoints[3,:] = [1 , 1, -1]
lPoints[4,:] = [-1 , 1, -1]

lPoints[5,:] = [-1 , -1, 1]
lPoints[6,:] = [1 , -1, 1]
lPoints[7,:] = [1 , 1, 1]
lPoints[8,:] = [-1 , 1, 1]
 
lPoints[9,:] = [0 , -1, -1]
lPoints[10,:] = [1 , 0, -1]
lPoints[11,:] = [0 , 1, -1]
lPoints[12,:] = [-1 , 0, -1]

lPoints[13,:] = [0 , -1, 1]
lPoints[14,:] = [1 , 0, 1]
lPoints[15,:] = [0 , 1, 1]
lPoints[16,:] = [-1 , 0, 1]

lPoints[17,:] = [-1 , -1, 0]
lPoints[18,:] = [1 , -1, 0]
lPoints[19,:] = [1 , 1, 0]
lPoints[20,:] = [-1 , 1, 0]

lPoints[21,:] = [0 , -1, 0]
lPoints[22,:] = [1 , 0, 0]
lPoints[23,:] = [0 , 1, 0]
lPoints[24,:] = [-1 , 0, 0]

lPoints[25,:] = [0, 0, -1]
lPoints[26,:] = [0, 0, 1]
lPoints[27,:] = [0, 0, 0]

eiter = 1:size(IEN,2)
# eiter = 1:1
nodeiter = 1:size(lPoints,1)

for e in eiter 
    cIter = 1:size(C,2)
    for i in cIter
        c = C[:,:,e]
        @test 1.0 ≈ sum(c[:,i]) atol=1e-5
    end
    for j in nodeiter
        L, ΔL = basis_function(lPoints[j,1],lPoints[j,2],lPoints[j,3], "Q2")

        # Compute the Bspline basis
        B = C[:,:,e]*L

        for b in B 
            if abs(b) > 1e-10
                @test (b > 0)
            end
        end 
        @test 1.0 ≈ sum(B) atol=1e-5

        NodeList[:,(e-1)*size(IEN,1)+j] = B'*CPointList[:,IEN[:,e]]'
        IEN_list[j,e] = (e-1)*size(IEN,1)+j
    end
end
display(IEN_list)
write_vtk(filePath, "q", NodeList, IEN_list, ne, ndim, q, ID=ID, FunctionClass="Q2")
