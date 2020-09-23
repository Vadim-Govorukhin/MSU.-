module Grids

export Domain2D, squareGrid, triangularGrid, randomPoints, remove_close!, turn!

struct Domain2D{T<:Real}
    xmin::T
    ymin::T
    xmax::T
    ymax::T
    Lx::T
    Ly::T
    function Domain2D(xmin::T, ymin::T, xmax::T, ymax::T) where {T<:Real}
        new{T}(xmin, ymin, xmax, ymax, xmax-xmin, ymax-ymin)
    end
end

Domain2D(;xmin::T, ymin::T, xmax::T, ymax::T) where {T<:Real} = Domain2D(xmin, ymin, xmax, ymax)

using LinearAlgebra: normalize, dot

function derivative(x::Real, F, h=1e-2)
    (F(x+h/2)-F(x-h/2))/h
end

function derivative2(x::Real, F, h=1e-2)
    (F(x-h)-2*F(x)+F(x+h))/h^2
end

function grad(p, f)
    [ derivative(p[1], x->f(x,p[2])),
      derivative(p[2], y->f(p[1],y)) ]
end

function Argmin(x, δx, f)::Real
    F = α->f((x.+α.*δx)...)
    αₙ = αₙ₊₁ = 0.0
    # abs(F(αₙ₊₁) - F(αₙ)) < 1e-2*abs(F(αₙ)) && return αₙ₊₁
    for iter = 1:10
        dF  = derivative(αₙ, F)
        d2F = derivative2(αₙ, F)
        αₙ₊₁ = αₙ - dF/d2F
        abs((F(αₙ₊₁) - F(αₙ))/F(αₙ)) < 1e-2 && break #(println("iters: $iter");break)
        αₙ = αₙ₊₁
    end
    return αₙ₊₁
end

function conjgrad(x₀, f)
    δx₀ = tuple(normalize(-1 .* grad(x₀, f))...)
    α₀ = Argmin(x₀, δx₀, f)
    x₁ = x₀ .+ α₀.*δx₀
    isnan(α₀) && (x₁ = x₀)

    δx₁ = tuple(normalize(-1 .* grad(x₁, f))...)
    β₁ = max(0, dot(δx₁, δx₁.-δx₀)/dot(δx₀, δx₀))
    s₁ = δx₁ .+ β₁.*δx₀
    α₁ = Argmin(x₁, s₁, f)
    x₂ = x₁ .+ α₁.*s₁
    isnan(α₁) && (x₂ = x₁)
    x₂
end

function squareGrid(::Type{T}, N::Int, f, D::Domain2D{<:Real}) where {T<:Real}
    Lˣ, Lʸ = D.Lx, D.Ly
    # h = (Lˣ+Lʸ)/2/(N-1) * (1 + sqrt(1 + 4*(N-1)*Lˣ*Lʸ/(Lˣ+Lʸ)^2))
    h = sqrt(Lˣ*Lʸ/N)
    Nˣ, Nʸ = round.(Int, [Lˣ/h,Lʸ/h].+1)
    xₘᵢₙ = D.xmin + (Lˣ-h*(Nˣ-1))/2
    yₘᵢₙ = D.ymin + (Lʸ-h*(Nʸ-1))/2

    println("N: $(Nˣ*Nʸ)")
    println("h: $h")

    # points = Vector{Tuple{T,T}}(undef, Nˣ*Nʸ)# Matrix{T}(undef, 2, Nˣ*Nʸ)
    points = [(xₘᵢₙ+h*(i-1),yₘᵢₙ+h*(j-1)) for i = 1:Nˣ for j = 1:Nʸ]
    # for i = 1:Nˣ, j = 1:Nʸ
    #     points[:, (i-1)*Nʸ+j] .= [xₘᵢₙ, yₘᵢₙ] + [h*(i-1), h*(j-1)]
    # end

    for (i, p) in enumerate(points)
        if abs(f(p...)) <= h/2
            points[i] = conjgrad(p, (x,y)->f(x,y)^2)
            # println("p: $p")
        end
    end

    points, h
end
squareGrid(N::Int, f, D::Domain2D{<:Real}) = squareGrid(Float64, N, f, D)


function triangularGrid(::Type{T}, N::Int, f, D::Domain2D{<:Real}) where {T<:Real}
    Lˣ, Lʸ = D.Lx, D.Ly
    h₀ = sqrt(2*Lˣ*Lʸ/√3N)
    l = 1 + 2*round(Int, Lʸ/√3h₀)
    h = 2Lʸ/√3(l-1)
    Nˡ = 1 + round(Int, Lˣ/h)

    xₘᵢₙ = D.xmin + (Lˣ-h*(Nˡ-1))/2
    yₘᵢₙ = D.ymin

    println("N: $(l*Nˡ-l÷2)")
    println("l: $l, Nˡ: $Nˡ, h: $h")

    # points = Matrix{T}(undef, 2, l*Nˡ-l÷2)
    # points = zeros(2, l*Nˡ-l÷2)
    points = [(xₘᵢₙ+h*(j-1)+h/2*(1-i%2),yₘᵢₙ+√3/2*h*(i-1)) for i = 1:l for j = 1:(Nˡ+i%2-1)]

    # for i = 1:l, j = 1:(Nˡ+i%2-1)
    #     points[:,(i-1)*Nˡ-(i-1)÷2+j] .= [xₘᵢₙ, yₘᵢₙ] + [h*(j-1)+h/2*(1-i%2), √3/2*h*(i-1)]
    # end

    for (i, p) in enumerate(points)
        if abs(f(p...)) <= h/2
            points[i] = conjgrad(p, (x,y)->f(x,y)^2)
            # println("p: ", points[i])
        end
    end

    points, h
end
triangularGrid(N::Int, f, D::Domain2D{<:Real}) = triangularGrid(Float64, N, f, D)


using RandomExtensions: rand

function randomPoints(::Type{T}, N::Int, f, D::Domain2D{<:Real}) where {T<:Real}
    points = Vector{Tuple{T,T}}(undef, N)
    for j = 1:N
        while true
            p = (D.xmin,D.ymin) .+ (D.Lx,D.Ly).*rand(Tuple{T,T})
            f(p...) < 0 && (points[j] = p; break)
        end
    end
    points
end
randomPoints(N::Int, f, D::Domain2D{<:Real}) = randomPoints(Float64, N, f, D)


function remove_close!(p::Vector, h::Real, debug::Bool=false)
    delete_list = Int[]
    for i in 1:size(p,1), j in i+1:size(p,1)
        if i != j && sum((p[i].-p[j]).^2) < h^2
            p[i] = (p[i].+p[j])./2
            p[j] = p[i]
            push!(delete_list, j)
        end
    end
    sort!(delete_list)
    debug && println("delete list: ", delete_list)
    deleteat!(p, delete_list)
end

function turn!(p::Vector, θ::Real)
    iszero(θ) && return
    for i in axes(p,1)
        p[i] = (p[i][1]*cos(θ)-p[i][2]*sin(θ),
                p[i][1]*sin(θ)+p[i][2]*cos(θ))
    end
end

end # module
