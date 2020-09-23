include("../src/Grids.jl")
# include("../src/domain.jl")

using .Grids
using Plots; pyplot()


function main(N; f)
    D = Domain2D(xmin=-1.0, ymin=-1.0, xmax=1.0, ymax=1.0)
    p, h = triangularGrid(N, f, D) #squareGrid(N, f, D)
    remove_close!(p, 0.5*h)
    # filter!(el->f(el...)<=0, p)

    img = contour(D.xmin-0.01 : D.Lx/200 : D.xmax+0.01,
                  D.ymin-0.01 : D.Ly/200 : D.ymax+0.01,
                  f, levels=[-h/2,0,h/2,100], #border=:none,
                  xlim=(D.xmin-D.Lx/200,D.xmax).*1.1, lw=0.3, ratio=:equal, #fill=true,
                  ylim=(D.ymin-D.Ly/200,D.ymax).*1.1, lc=:black, legend=:none)

    scatter!(getfield.(p,1), getfield.(p,2), label="", ratio=:equal)
    for point in p
        if -h/2 <= f(point...) <= 0 #abs(f(point...)) <= h/2
            scatter!([point[1]], [point[2]], label="", c=:green)
        end
    end
    # for j=1:size(p)
    #     println(p[j])
    # end
    img
end



circle(x::T, y::T, r::Real) where {T<:Real} = sqrt(x^2+y^2)-r
ring(x::T, y::T) where {T<:Real} = max(circle(x,y,T(1)),-circle(x,y,T(0.5)))
ring2(x::T, y::T) where {T<:Real} = min(ring(x,y),ring(x-1.5,y))
ring3(x::T, y::T) where {T<:Real} = min(ring(x+0.75,y),ring(x-0.75,y),ring(x,y-0.75*sqrt(3)))
function rectangle(x::T, y::T, a::Real, b::Real) where {T<:Real}
    q = abs.([x,y]) .- [a/2, b/2]
    sqrt(sum(max.(q,0).^2)) + min(maximum(q),0)
end
ellipse(x, y, a, b) = sqrt((b*x)^2+(a*y)^2) - a*b


main(200, f = (x,y)->ellipse(x,y,1.0,0.5))
# main(200, f = (x,y)->circle(x,y,1.0))
# savefig("../img/trigrid2.png")
