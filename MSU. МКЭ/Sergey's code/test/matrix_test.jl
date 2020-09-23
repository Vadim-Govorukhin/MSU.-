
function K(triSet::Vector, p::Vector)

	N = size(p,1)
    
	k = spzeros(N,N)
    
	for tri in triSet
        
		tri = collect(tri)
        
		s = Sw(tri..., p)
        
		for (i,α) in enumerate(tri)
            aⁱ, bⁱ, _ = abc(tri..., i, p)
            for (j,β) in enumerate(tri)
                aʲ, bʲ, _ = abc(tri..., j, p)
                k[α,β] += (aⁱ*aʲ+bⁱ*bʲ)*s
            end
        end
    end
    k
end



# из диссертации Коняева
function M(triSet::Vector, p::Vector)

    N = size(p,1)

    m = spzeros(N,N)

    X, Y = zeros(3), zeros(3)

    for tri in triSet

        tri = collect(tri)

        s = Sw(tri..., p)

        for (i,α) in enumerate(tri)

            X[i], Y[i] = p[α].x, p[α].y

        end

        x₀ = sum(X) / 3

        y₀ = sum(Y) / 3

        X .-= x₀

        Y .-= y₀

        for (i,α) in enumerate(tri)

            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)

            for (j,β) in enumerate(tri)

                aʲ, bʲ, cʲ = abc(tri..., j, p)

                α₁, α₂, α₃ = aⁱ*aʲ, bⁱ*bʲ, cⁱ*cʲ

                α₄, α₅, α₆ = aⁱ*bʲ+aʲ*bⁱ, aⁱ*cʲ+aʲ*cⁱ, bⁱ*cʲ+bʲ*cⁱ

                β₃ = α₃ + α₁*x₀^2 + α₂*y₀^2 + α₄*x₀*y₀ + α₅*x₀ + α₆*y₀

                m[α,β] += s/12*(α₁*sum(X.^2) + α₂*sum(Y.^2) + α₄*dot(X,Y)) + s*β₃

            end

        end

    end

    m

end



# в декартовых с нулем в (0,0)
function M2(triSet::Vector, p::Vector)
    N = size(p,1)
    m = spzeros(N,N)
    for tri in triSet
        tri = collect(tri)
        s = Sw(tri..., p)
        t = collect(tri)
        x₁, y₁ = p[t[1]].x, p[t[1]].y
        x₂, y₂ = p[t[2]].x, p[t[2]].y
        x₃, y₃ = p[t[3]].x, p[t[3]].y
        for (i,α) in enumerate(tri)
            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)
            for (j,β) in enumerate(tri)
                aʲ, bʲ, cʲ = abc(tri..., j, p)
                α₁ = aⁱ*bʲ + bⁱ*aʲ
                α₂ = aⁱ*aʲ
                α₃ = bⁱ*bʲ
                α₄ = aⁱ*cʲ+cⁱ*aʲ
                α₅ = bⁱ*cʲ+cⁱ*bʲ
                I = α₁*(x₁*y₁ + x₂*y₂ + x₃*y₃) -
                   2α₂*(x₁*x₂ + x₁*x₃ + x₂*x₃) -
                   2α₃*(y₁*y₂ + y₁*y₃ + y₂*y₃) +
                    α₁*(x₁ + x₂ + x₃)*(y₁ + y₂ + y₃) +
                   2α₂*(x₁ + x₂ + x₃)^2 +
                   2α₃*(y₁ + y₂ + y₃)^2 +
                   4α₄*(x₁ + x₂ + x₃) +
                   4α₅*(y₁ + y₂ + y₃) +
                  12cⁱ*cʲ
                m[α,β] += s/12 * I
            end
        end
    end
    m
end

# в декартовых с нулем в центре масс треугольника
function M3(triSet::Vector, p::Vector)
    N = size(p,1)
    m = spzeros(N,N)
    X, Y = zeros(3), zeros(3)
    for tri in triSet
        tri = collect(tri)
        s = Sw(tri..., p)
        for (i,α) in enumerate(tri)
            X[i], Y[i] = p[α].x, p[α].y
        end
        x₀ = sum(X) / 3
        y₀ = sum(Y) / 3
        X .-= x₀
        Y .-= y₀
        for (i,α) in enumerate(tri)
            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)
            for (j,β) in enumerate(tri)
                aʲ, bʲ, cʲ = abc(tri..., j, p)
                α₁ = aⁱ*aʲ
                α₂ = bⁱ*bʲ
                α₃ = aⁱ*bʲ + bⁱ*aʲ
                α₄ = (aⁱ*x₀ + bⁱ*y₀ + cⁱ)*(aʲ*x₀ + bʲ*y₀ + cʲ)
                I = α₁*dot(X,X) + α₂*dot(Y,Y) + α₃*dot(X,Y) + 12α₄
                m[α,β] += s/12 * I
            end
        end
    end
    m
end



# из диссертации Коняева

function T(triSet::Vector, p::Vector, q::Function)

    N = size(p,1)

    t = spzeros(N,N)

    X, Y = zeros(3), zeros(3)

    for tri in triSet

        tri = collect(tri)

        s = Sw(tri..., p)

        for (i,α) in enumerate(tri)

            X[i], Y[i] = p[α].x, p[α].y

        end

        x₀ = sum(X) / 3

        y₀ = sum(Y) / 3

        X .-= x₀

        Y .-= y₀

        for (i,α) in enumerate(tri)

            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)

            for (j,β) in enumerate(tri)

                aʲ, bʲ, cʲ = abc(tri..., j, p)

                α₁, α₂, α₃ = aⁱ*aʲ, bⁱ*bʲ, cⁱ*cʲ

                α₄, α₅, α₆ = aⁱ*bʲ+aʲ*bⁱ, aⁱ*cʲ+aʲ*cⁱ, bⁱ*cʲ+bʲ*cⁱ

                for (l,γ) in enumerate(tri)

                    aˡ, bˡ, cˡ = abc(tri..., l, p)

                    γ₁, γ₂, γ₃ = α₁*cˡ+α₅*aˡ, α₂*cˡ+α₆*bˡ, α₃*cˡ

                    γ₄ = α₄*cˡ+α₅*bˡ+α₆*aˡ

                    γ₅, γ₆, γ₇, γ₈ = α₅*cˡ+α₃*aˡ, α₆*cˡ+α₃*bˡ, α₁*aˡ, α₂*bˡ

                    γ₉, γ₁₀ = α₄*aˡ+α₁*bˡ, α₄*bˡ+α₂*aˡ

                    η₁, η₂ = γ₁+3γ₇*x₀+γ₉*y₀, γ₂+3γ₈*y₀+γ₁₀*x₀

                    η₃ = γ₃ + γ₁*x₀^2 + γ₂*y₀^2 + γ₄*x₀*y₀ + γ₅*x₀ + γ₆*y₀ +

                         γ₇*x₀^3 + γ₈*y₀^3 + γ₉*x₀^2*y₀ + γ₁₀*x₀*y₀^2

                    η₄ = γ₄ + 2γ₉*x₀ + 2γ₁₀*y₀

                    η₅ = γ₅ + 2γ₁*x₀ + γ₄*y₀ + 3γ₇*x₀^2 + 2γ₉*x₀*y₀ + γ₁₀*y₀^2

                    η₆ = γ₆ + 2γ₂*y₀ + γ₄*x₀ + 3γ₈*y₀^2 + 2γ₁₀*x₀*y₀ + γ₉*y₀^2

                    t[α,β] += q(p[γ].x,p[γ].y) * ( s/30 *

                        (γ₇*sum(X) + γ₈*sum(Y) + γ₉*dot(X.^2,Y) + γ₁₀*dot(X,Y.^2)) +

                        s/12 * (η₁*sum(X.^2) + η₂*sum(Y.^2) + η₄*dot(X,Y)) + s*η₃)

                end

            end

        end

    end

    t

end




# в декартовых с нулем в (0,0)
function T2(triSet::Vector, p::Vector, q::Function)
    N = size(p,1)
    t = spzeros(N,N)
    X, Y = zeros(3), zeros(3)
    for tri in triSet
        tri = collect(tri)
        s = Sw(tri..., p)
        for (i,α) in enumerate(tri)
            X[i], Y[i] = p[α].x, p[α].y
        end
        for (i,α) in enumerate(tri)
            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)
            for (j,β) in enumerate(tri)
                aʲ, bʲ, cʲ = abc(tri..., j, p)
                for (l,γ) in enumerate(tri)
                    aˡ, bˡ, cˡ = abc(tri..., l, p)
                    α₁ = bⁱ*aʲ*aˡ + aⁱ*bʲ*aˡ + aⁱ*aʲ*bˡ
                    α₂ = aⁱ*bʲ*bˡ + bⁱ*aʲ*bˡ + bⁱ*bʲ*aˡ
                    α₃ = aⁱ*bʲ*cˡ + bⁱ*aʲ*cˡ + bⁱ*cʲ*aˡ +
                         aⁱ*cʲ*bˡ + cⁱ*aʲ*bˡ + cⁱ*bʲ*aˡ
                    α₄ = aⁱ*aʲ*aˡ
                    α₅ = bⁱ*bʲ*bˡ
                    α₆ = cⁱ*aʲ*aˡ + aⁱ*cʲ*aˡ + aⁱ*aʲ*cˡ
                    α₇ = cⁱ*bʲ*bˡ + bⁱ*cʲ*bˡ + bⁱ*bʲ*cˡ
                    α₈ = aⁱ*cʲ*cˡ + cⁱ*aʲ*cˡ + cⁱ*cʲ*aˡ
                    α₉ = bⁱ*cʲ*cˡ + cⁱ*bʲ*cˡ + cⁱ*cʲ*bˡ

                    I = 2α₁*(sum(X)^2*sum(Y) + 2*dot(X.^2,Y) - X[1]*X[2]*Y[3] - X[1]*X[3]*Y[2] - X[2]*X[3]*Y[1]) +
                        2α₂*(sum(X)*sum(Y)^2 + 2*dot(X,Y.^2) - X[3]*Y[1]*Y[3] - X[2]*Y[1]*Y[3] - X[1]*Y[2]*Y[3]) +
                        5α₃*(sum(X)*sum(Y) + dot(X,Y)) +
                        6α₄*(sum(X)^3 - 2*dot(X,[X[2]+X[3],X[1]+X[3],X[1]+X[2]].^2) + 7X[1]*X[2]*X[3]) +
                        6α₅*(sum(Y)^3 - 2*dot(Y,[Y[2]+Y[3],Y[1]+Y[3],Y[1]+Y[2]].^2) + 7Y[1]*Y[2]*Y[3]) +
                       10α₆*(sum(X)^2 - (X[1]*X[2] + X[1]*X[3] + X[2]*X[3])) +
                       10α₇*(sum(Y)^2 - (Y[1]*Y[2] + Y[1]*Y[3] + Y[2]*Y[3])) +
                       20α₈*sum(X) +
                       20α₉*sum(Y) +
                       60cⁱ*cʲ*cˡ
                    t[α,β] += q(p[γ].x,p[γ].y) * s/60 * I
                end
            end
        end
    end
    t
end

# в декартовых с нулем в центре масс треугольника

function T3(triSet::Vector, p::Vector, q::Function)
    N = size(p,1)
    t = spzeros(N,N)
    X, Y = zeros(3), zeros(3)
    for tri in triSet
        tri = collect(tri)
        s = Sw(tri..., p)

        for (i,α) in enumerate(tri)

            X[i], Y[i] = p[α].x, p[α].y

        end
        x₀ = sum(X) / 3

        y₀ = sum(Y) / 3

        X .-= x₀

        Y .-= y₀

        for (i,α) in enumerate(tri)

            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)

            for (j,β) in enumerate(tri)

                aʲ, bʲ, cʲ = abc(tri..., j, p)

                α₁, α₂, α₃ = aⁱ*aʲ, bⁱ*bʲ, cⁱ*cʲ
                α₄, α₅, α₆ = aⁱ*bʲ+aʲ*bⁱ, aⁱ*cʲ+aʲ*cⁱ, bⁱ*cʲ+bʲ*cⁱ
                for (l,γ) in enumerate(tri)
                    aˡ, bˡ, cˡ = abc(tri..., l, p)
                    α₁ = bⁱ*aʲ*aˡ + aⁱ*bʲ*aˡ + aⁱ*aʲ*bˡ
                    α₂ = aⁱ*bʲ*bˡ + bⁱ*aʲ*bˡ + bⁱ*bʲ*aˡ
                    α₃ = aⁱ*bʲ*cˡ + bⁱ*aʲ*cˡ + bⁱ*cʲ*aˡ +
                         aⁱ*cʲ*bˡ + cⁱ*aʲ*bˡ + cⁱ*bʲ*aˡ +
                         2α₁*x₀ + 2α₂*y₀
                    α₄ = aⁱ*aʲ*aˡ
                    α₅ = bⁱ*bʲ*bˡ
                    α₆ = cⁱ*aʲ*aˡ + aⁱ*cʲ*aˡ + aⁱ*aʲ*cˡ + 3α₄*x₀ + α₁*y₀
                    α₇ = cⁱ*bʲ*bˡ + bⁱ*cʲ*bˡ + bⁱ*bʲ*cˡ + 3α₅*y₀ + α₂*x₀
                    α₈ = aⁱ*cʲ*cˡ + cⁱ*aʲ*cˡ + cⁱ*cʲ*aˡ
                    α₉ = bⁱ*cʲ*cˡ + cⁱ*bʲ*cˡ + cⁱ*cʲ*bˡ
                    α₁₀ = (aⁱ*x₀ + bⁱ*y₀ + cⁱ) *
                          (aʲ*x₀ + bʲ*y₀ + cʲ) *
                          (aˡ*x₀ + bˡ*y₀ + cˡ)
                    I = 2α₁*dot(X.^2,Y) + 2α₂*dot(X,Y.^2) +
                        5α₃*dot(X,Y) +
                        6α₄*X[1]*X[2]*X[3] +
                        6α₅*Y[1]*Y[2]*Y[3] +
                        5α₆*dot(X,X) + 5α₇*dot(Y,Y) + 60*α₁₀
                    t[α,β] += q(p[γ].x,p[γ].y) * s/60 * I
                end
            end
        end
    end
    t

end




L(x1, y1, x2, y2) = √((x2-x1)^2 + (y2-y1)^2)


# в декартовых с нулем в (0,0)

function P(bSet::Matrix, p::Vector, pf::Function)

    N = size(p,1)

    Pij = spzeros(N,N)

    for tri in eachcol(bSet)

        tri = collect(tri)

        x₁, y₁ = p[tri[1]].x, p[tri[1]].y
        x₂, y₂ = p[tri[2]].x, p[tri[2]].y

        len = L(x₁, y₁, x₂, y₂)

        for (i,α) in enumerate(tri[1:2])

            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)

            for (j,β) in enumerate(tri[1:2])

                aʲ, bʲ, cʲ = abc(tri..., j, p)

                for (l,γ) in enumerate(tri[1:2])

                    aˡ, bˡ, cˡ = abc(tri..., l, p)

                    α₁ = aⁱ*aʲ*bˡ + aⁱ*bʲ*aˡ + bⁱ*aʲ*aˡ

                    α₂ = aⁱ*bʲ*bˡ + bⁱ*aʲ*bˡ + bⁱ*bʲ*aˡ

                    α₃ = aⁱ*bʲ*cˡ + bⁱ*aʲ*cˡ + bⁱ*cʲ*aˡ +

                         aⁱ*cʲ*bˡ + cⁱ*aʲ*bˡ + cⁱ*bʲ*aˡ

                    α₄ = aⁱ*aʲ*aˡ

                    α₅ = bⁱ*bʲ*bˡ

                    α₆ = aⁱ*aʲ*cˡ + aⁱ*cʲ*aˡ + cⁱ*aʲ*aˡ

                    α₇ = bⁱ*bʲ*cˡ + bⁱ*cʲ*bˡ + cⁱ*bʲ*bˡ

                    α₈ = aⁱ*cʲ*cˡ + cⁱ*aʲ*cˡ + cⁱ*cʲ*aˡ

                    α₉ = bⁱ*cʲ*cˡ + cⁱ*bʲ*cˡ + cⁱ*cʲ*bˡ

                    I = α₁*(2x₁^2*y₁ + (y₁+y₂)*(x₁+x₂)^2 + 2x₂^2*y₂) +

                        α₂*(2x₁*y₁^2 + (y₁+y₂)^2*(x₁+x₂) + 2x₂*y₂^2) +

                       2α₃*(2x₁*y₁ + x₁*y₂ + x₂*y₁ + 2x₂*y₂) +

                       3α₄*(x₁+x₂)*(x₁^2 + x₂^2) +

                       3α₅*(y₁+y₂)*(y₁^2 + y₂^2) +

                       4α₆*(x₁^2 + x₁*x₂ + x₂^2) +

                       4α₇*(y₁^2 + y₁*y₂ + y₂^2) +

                       6α₈*(x₁+x₂) +

                       6α₉*(y₁+y₂) +

                      12cⁱ*cʲ*cˡ


                    Pij[α,β] += pf(p[γ].x,p[γ].y) * len/12 * I

                end

            end

        end

    end

    Pij

end



# в декартовых с нулем в центре масс треугольника

function P2(bSet::Matrix, p::Vector, pf::Function)
    N = size(p,1)
    Pij = spzeros(N,N)
    for tri in eachcol(bSet)
        tri = collect(tri)
        x₁, y₁ = p[tri[1]].x, p[tri[1]].y
        x₂, y₂ = p[tri[2]].x, p[tri[2]].y
        len = L(x₁, y₁, x₂, y₂)
        x₀, y₀ = (x₁+x₂)/2, (y₁+y₂)/2
        x₁, x₂ = x₁-x₀, x₂-x₀
        y₁, y₂ = y₁-y₀, y₂-y₀
        for (i,α) in enumerate(tri[1:2])
            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)
            for (j,β) in enumerate(tri[1:2])
                aʲ, bʲ, cʲ = abc(tri..., j, p)
                for (l,γ) in enumerate(tri[1:2])
                    aˡ, bˡ, cˡ = abc(tri..., l, p)
                    α₁ = bⁱ*aʲ*aˡ + aⁱ*bʲ*aˡ + aⁱ*aʲ*bˡ
                    α₂ = aⁱ*bʲ*bˡ + bⁱ*aʲ*bˡ + bⁱ*bʲ*aˡ
                    α₃ = cⁱ*aʲ*aˡ + aⁱ*cʲ*aˡ + aⁱ*aʲ*cˡ
                    α₄ = cⁱ*bʲ*bˡ + bⁱ*cʲ*bˡ + bⁱ*bʲ*cˡ
                    α₅ = aⁱ*bʲ*cˡ + aⁱ*cʲ*bˡ + bⁱ*aʲ*cˡ +
                         bⁱ*cʲ*aˡ + cⁱ*aʲ*bˡ + cⁱ*bʲ*aˡ
                    α₆ = aⁱ*cʲ*cˡ + cⁱ*aʲ*cˡ + cⁱ*cʲ*aˡ
                    α₇ = cⁱ*bʲ*cˡ + bⁱ*cʲ*cˡ + cⁱ*cʲ*bˡ
                    α₈ = aⁱ*aʲ*aˡ
                    α₉ = bⁱ*bʲ*bˡ
                    β₁ = α₁*y₀ + α₃ + 3α₈*x₀
                    β₂ = α₂*x₀ + α₄ + 3α₉*y₀
                    β₃ = 2α₁*x₀ + 2α₂*y₀ + α₅
                    β₄ = (aⁱ*x₀ + bⁱ*y₀ + cⁱ) *
                         (aʲ*x₀ + bʲ*y₀ + cʲ) *
                         (aˡ*x₀ + bˡ*y₀ + cˡ)
                    I = β₁*x₁^2 + β₂*y₁^2 + β₃*x₁*y₁ + 3β₄
                    Pij[α,β] += pf(p[γ].x,p[γ].y) * len/3 * I
                end
            end
        end
    end
    Pij
end

# в декартовых с нулем в (0,0)
function H(bSet::Matrix, p::Vector, hf::Function)
    N = size(p,1)
    h = zeros(N)
    for tri in eachcol(bSet)
        tri = collect(tri)
        x₁, y₁ = p[tri[1]].x, p[tri[1]].y
        x₂, y₂ = p[tri[2]].x, p[tri[2]].y
        l = L(x₁, y₁, x₂, y₂)
        for (i,α) in enumerate(tri[1:2])
            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)
            for (j,β) in enumerate(tri[1:2])
                aʲ, bʲ, cʲ = abc(tri..., j, p)
                α₁ = aⁱ*bʲ + bⁱ*aʲ
                α₂ = aⁱ*aʲ
                α₃ = bⁱ*bʲ
                α₄ = aⁱ*cʲ + cⁱ*aʲ
                α₅ = bⁱ*cʲ + cⁱ*bʲ
                I = α₁*(2x₁*y₁ + x₁*y₂ + x₂*y₁ + 2x₂*y₂) +
                   2α₂*(x₁^2 + x₁*x₂ + x₂^2) +
                   2α₃*(y₁^2 + y₁*y₂ + y₂^2) +
                   3α₄*(x₁ + x₂) +
                   3α₅*(y₁ + y₂) +
                   6*cⁱ*cʲ
                h[α] += hf(p[β].x,p[β].y) * l/6 * I
            end
        end
    end
    h
end

# в декартовых с нулем в центре масс треугольника
function H2(bSet::Matrix, p::Vector, hf::Function)
    N = size(p,1)
    h = zeros(N)
    for tri in eachcol(bSet)
        tri = collect(tri)
        x₁, y₁ = p[tri[1]].x, p[tri[1]].y
        x₂, y₂ = p[tri[2]].x, p[tri[2]].y
        len = L(x₁, y₁, x₂, y₂)
        x₀, y₀ = (x₁+x₂)/2, (y₁+y₂)/2
        x₁, x₂ = x₁-x₀, x₂-x₀
        y₁, y₂ = y₁-y₀, y₂-y₀
        for (i,α) in enumerate(tri[1:2])
            aⁱ, bⁱ, cⁱ = abc(tri..., i, p)
            for (j,β) in enumerate(tri[1:2])
                aʲ, bʲ, cʲ = abc(tri..., j, p)
                α₁ = aⁱ*aʲ
                α₂ = bⁱ*bʲ
                α₃ = aⁱ*bʲ + bⁱ*aʲ
                α₄ = (aⁱ*x₀ + bⁱ*y₀ + cⁱ)*(aʲ*x₀ + bʲ*y₀ + cʲ)
                I = α₁*x₁^2 + α₂*y₁^2 + α₃*x₁*y₁ + 3α₄
                h[α] += hf(p[β].x,p[β].y) * len/3 * I
            end
        end
    end
    h
end

# u(x,y) = x^2 + y^2
# q(x,y) = 1
# f(x,y) = x^2 + y^2 - 4
# h(x,y) = 2
# p(x,y) = 0

# u(x,y) = cos(x^2+y^2)
# q(x,y) = 1
# f(x,y) = 4*(x^2+y^2)*cos(x^2+y^2) + 4*sin(x^2+y^2) + cos(x^2+y^2)
# p(x,y) = 0
# h(x,y) = -2*sin(1)

# u(x,y) = cos(π*(x^2+y^2))
# q(x,y) = 1
# f(x,y) = 4*(x^2+y^2)*pi^2*cos(pi*(x^2 + y^2)) +
#          4*pi*sin(pi*(x^2 + y^2)) + cos(pi*(x^2 + y^2))
# p(x,y) = 10
# h(x,y) = -10

# u(x,y) = x^2+y^2
# f(x,y) = -4
# q(x,y) = 0
# p(x,y) = -10
# h(x,y) = 2-10

# u(x,y) = x^2
# f(x,y) = 10x^2-2
# q(x,y) = 10
# p(x,y) = -2
# h(x,y) = 0

# u(x,y) = y^2
# f(x,y) = 10y^2-2
# q(x,y) = 10
# p(x,y) = -2
# h(x,y) = 0

# u(x,y) = y^2
# f(x,y) = 2y^2-2
# q(x,y) = 2
# p(x,y) = 0
# h(x,y) = 2y^2

# u(x,y) = x^2
# f(x,y) = 2x^2-2
# q(x,y) = 2
# p(x,y) = 0
# h(x,y) = 2x^2

# u(x,y) = exp(-x^2)
# q(x,y) = 4x^2
# f(x,y) = 2*exp(-x^2)
# p(x,y) = 0
# h(x,y) = -2x^2*exp(-x^2)

# u(x,y) = exp(-x^2)
# q(x,y) = 4x^2
# f(x,y) = 2*exp(-x^2)
# p(x,y) = 3x^2
# h(x,y) = x^2*exp(-x^2)

u(x,y) = exp(-x^2 - y^2)
q(x,y) = -4
f(x,y) = -4*(x^2 + y^2)*exp(-x^2 - y^2)
p(x,y) = 2y^2
h(x,y) = -2x^2*exp(-x^2 - y^2)

using IterativeSolvers
using PyPlot

function test(;N=6, F=(x,y)->circle(x,y,1.0))
    D = Domain2D(xmin=-1.0, ymin=-1.0, xmax=1.0, ymax=1.0)
    points, Δh = triangularGrid(N, F, D)
    remove_close!(points, 0.5*Δh)
    points = sort(points, lt=pointIsLess)
    filter!(el->F(el...)<=0, points)
    Np = length(points)
    println("N: $N")

    v = [Point(points[j]..., j) for j=1:Np]

    dealunay(1, Np+1, v, true)

    # adjMat = adjacencyMatrix(v) # матрица смежности
    triSet = triangleSet(v) # набор треугольников
    bound = boundary(v) # индексы граничных узлов
    boundSet = boundarySet(triSet, bound) # набор граничных треугольников

    Kij = K(triSet, v)
    Fij = M3(triSet, v) * [f(V.x,V.y) for V in v]
    Tij = T3(triSet, v, q)
    Hij = H2(boundSet, v, h)
    Pij = P2(boundSet, v, p)

    A = Kij + Tij + Pij
    b = Fij + Hij
    # gm, min
    # C = gmres(A, b)
    U = A \ b


    X = getfield.(v,:x) #-1:0.05:1
    Y = getfield.(v,:y) #-1:0.05:1
    Uexact = [u(X[i],Y[i]) for i in axes(X,1)]

    clf()
    fig1, ax = plt.subplots(2, 2, figsize=(15, 15))
    ax[1,1].set_title("exact")
    ax[1,1].set_aspect(:equal)
    img1 = ax[1,1].tricontourf(X, Y, Uexact)
    colorbar(img1, ax=ax[1,1])

    ax[1,2].set_title("numerical")
    ax[1,2].set_aspect(:equal)
    img2 = ax[1,2].tricontourf(X, Y, U)
    colorbar(img2, ax=ax[1,2])

    ax[2,1].set_title("error")
    ax[2,1].set_aspect(:equal)
    img3 = ax[2,1].tricontourf(X, Y, abs.(Uexact-U))
    colorbar(img3, ax=ax[2,1])

    ax[2,2].set_title("mesh")
    ax[2,2].set_aspect(:equal)
    img4 = ax[2,2].triplot(X, Y, label="") # сетка
    gcf()

    # plt.savefig("img.png")
end

test(N=200) #, F = (x,y)->ellipse(x,y,1.0,0.5))
