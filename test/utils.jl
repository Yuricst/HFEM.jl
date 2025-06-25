"""Utils for test"""

using Printf

function print_matrix(A)
    for i in 1:size(A,1)
        for j in 1:size(A,2)
            @printf("% 1.6e  ", A[i,j])
        end
        println()
    end
end