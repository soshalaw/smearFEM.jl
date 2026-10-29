# usage: julia diff.jl a.jls b.jls
using Serialization, LinearAlgebra, SparseArrays

flat(x::Number) = [Float64(x)]
flat(x::AbstractArray{<:Number}) = vec(Float64.(collect(x)))
flat(x::AbstractArray) = isempty(x) ? Float64[] : reduce(vcat, flat.(x))

a, b = deserialize(ARGS[1]), deserialize(ARGS[2])
names = ["output", "gradList", "borderPts2D", "displacement", "surface_pts_3D",
         "pos2D", "pos3D", "velocity", "pressure", "gradList_3d"]
worst = 0.0
for (n, x, y) in zip(names, a, b)
    fx, fy = flat(x), flat(y)
    if length(fx) != length(fy)
        println(rpad(n, 16), "LENGTH MISMATCH ", length(fx), " vs ", length(fy))
        global worst = Inf
        continue
    end
    rel = norm(fx - fy) / max(norm(fx), eps())
    global worst = max(worst, rel)
    println(rpad(n, 16), "rel. diff = ", rel)
end
println("worst: ", worst)
