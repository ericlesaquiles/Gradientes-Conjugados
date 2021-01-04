using LinearAlgebra # para podermos usar I, a identidade, por conveniência e sem perda de eficiência computacional
function CG(A, b, x0, imax = size(A)[1], Minv = I, ϵ = 1e-5)
    i  = 0
    x  = x0
    r  = b - A*x
    d  = Minv*r
    δn = r'd     
    δ0 = δn
    while i < imax && δn > ϵ^2*δ0
        q  = A*d
        α  = δn/(d'q)
        x  = x + α*d
        if (i + 1) ÷ 50 == 0  # por questão de estabilidade numérica, recalculamos r a cada 50 iterações
            r = b - A*x
        else
            r  = r - α*q
        end
        s  = Minv*r
        δv = δn       # δ velho
        δn = r's      # δ novo
        β  = δn/δv
        d  = s + β*d
        i += 1
    end
    return x
end

A  = [3 2; 2 6]

b  = [2; -8]
x0 = [14; -20]

x = CG(A, b, x0)

A\b

using BenchmarkTools, MatrixDepot, IterativeSolvers, LinearAlgebra, SparseArrays

# Matriz de Wathen de dimensões 30401 x 30401
A = matrixdepot("wathen", 100)

using UnicodePlots
spy(A)

# Nível de esparsidade
count(!iszero, A) / length(A)

b = ones(size(A, 1))
# Resolve Ax=b by CG
xcg = cg(A, b);
@benchmark cg($A, $b)

using Preconditioners
@time p = CholeskyPreconditioner(A, 2)

xpcg = cg(A, b, Pl=p)
# same answer?
norm(xcg - xpcg)

@benchmark cg($A, $b, Pl=$p)
