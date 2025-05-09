
function mat_nan_inf_check(v::AbstractArray)
    row, col = size(v,1), size(v,2)
    for i in 1:row
        for j in 1:col
            if isnan(v[i,j])
                throw(ArgumentError("NaN value found in the matrix"))
            elseif isinf(v[i,j])
                throw(ArgumentError("Inf value found in the matrix"))
            end
        end
    end
end