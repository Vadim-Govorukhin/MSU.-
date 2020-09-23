include("../src/Grids.jl")
using .Grids: Domain2D, squareGrid, triangularGrid, remove_close!, turn!
include("../src/Dealunay.jl")
using .Dealunay: insert, triangle, hull, mergeHulls, pointIsLess, turn!

using Plots; pyplot()

function dealunay_test(L::Int, R::Int, p::Vector,
    film::Animation, debug::Bool=false) where {T<:Real}

    debug && println("L: $L, R: $R")

    L + 3 == R && (triangle(L, L+1, L+2, p); return)
    L + 2 == R && (insert(L, L+1, p); return)

    mid = (L + R) ÷ 2
    dealunay_test(L, mid, p, film, false)
    dealunay_test(mid, R, p, film, false)

    plot(xlim=(-1,1).*1.1, lw=0.3, ratio=:equal,
         ylim=(-1,1).*1.1, lc=:black, label="")
    for point in p
        for n in point.list
            plot!([point.x, p[n].x], [point.y, p[n].y], lc=:black, label="")
        end
    end
    scatter!([0], [0], label="before merge")
    frame(film)

    if L < mid < R
        BT = hull((L,mid-1), (mid,R-1), p)
        UT = reverse(hull((R-1,mid), (mid-1,L), p))
        mergeHulls(BT, UT, p)
    end

    plot(xlim=(-1,1).*1.1, lw=0.3, ratio=:equal,
         ylim=(-1,1).*1.1, lc=:black, label="")
    for point in p
        for n in point.list
            plot!([point.x, p[n].x], [point.y, p[n].y], lc=:black, label="")
        end
    end
    scatter!([0], [0], label="after merge")
    frame(film)
end


function dealunay_gif(;N=6, f=(x,y)->circle(x,y,1.0))
    D = Domain2D(xmin=-1.0, ymin=-1.0, xmax=1.0, ymax=1.0)
    p, h = triangularGrid(N, f, D)
    p = sort(p, lt=pointIsLess)
    filter!(el->f(el...)<=0, p)
    N = length(p);   println("N: $N")
    v = [Dealunay.Point(p[j]..., j) for j=1:length(p)]

    film = Animation()
    img = contour(D.xmin-0.01 : D.Lx/200 : D.xmax+0.01,
              D.ymin-0.01 : D.Ly/200 : D.ymax+0.01,
              f, levels=[-h/2,0,h/2,100], #border=:none,
              xlim=(D.xmin-D.Lx/200,D.xmax).*1.1, lw=0.3, ratio=:equal, #fill=true,
              ylim=(D.ymin-D.Ly/200,D.ymax).*1.1, lc=:black, legend=:none)
    for point in v
        for n in point.list
            plot!([point.x, v[n].x], [point.y, v[n].y], lc=:black)
        end
    end
    frame(film)

    dealunay_test(1, N+1, v, film, false)

    gif(film, "../img/test3.gif", fps=2)
end


function dealunay_test_2(L::Int, R::Int, p::Vector, debug::Bool=false)

    debug && println("L: $L, R: $R")
    L + 1 < R || println("shiiiiiiit")
# if L < R-1
    L + 3 == R && (triangle(L, L+1, L+2, p); return)
    L + 2 == R && (insert(L, L+1, p); return)

    # debug && println("L: $L, R: $R - out")

    mid = (L + R) ÷ 2
    dealunay_test_2(L, mid, p, debug)
    dealunay_test_2(mid, R, p, debug)

    if L < mid < R
        BT = hull((L,mid-1), (mid,R-1), p, true)
        UT = reverse(hull((R-1,mid), (mid-1,L), p, true))
        mergeHulls(BT, UT, p, true)
    end
# end
end

## картинка

function dealunay_img(;N=6, f=(x,y)->circle(x,y,1.0))
    D = Domain2D(xmin=-1.0, ymin=-1.0, xmax=1.0, ymax=1.0)
    p, h = triangularGrid(N, f, D) #squareGrid(N, f, D)
    remove_close!(p, 0.5*h)
    p = sort(p, lt=pointIsLess)
    filter!(el->f(el...)<=0, p)

    Np = length(p);   println("total N: $Np")
    θ = 0 #π/11 #-π/20*N/500 # π/19
    turn!(p, -θ)
    v = [Dealunay.Point(p[j]..., j) for j=1:length(p)]
    turn!(p, θ)

    println("Traingulation started")
    dealunay_test_2(1, Np+1, v, true)
    turn!(v, θ)
    println("Traingulation completed")

    img = contour(D.xmin-0.01 : D.Lx/200 : D.xmax+0.01,
              D.ymin-0.01 : D.Ly/200 : D.ymax+0.01,
              f, levels=[-h/2,0,h/2,100], #border=:none,
              xlim=(D.xmin-D.Lx/200,D.xmax).*1.1, lw=0.3, ratio=:equal, #fill=true,
              ylim=(D.ymin-D.Ly/200,D.ymax).*1.1, lc=:black, legend=:none)
    for point in v
      for n in point.list
          plot!([point.x, v[n].x], [point.y, v[n].y], lc=:black)
      end
    end
    scatter!(getfield.(p,1), getfield.(p,2), label="", ratio=:equal)
    for point in p
      if -h/2 <= f(point...) <= 0 #abs(f(point...)) <= h/2
          scatter!([point[1]], [point[2]], label="", c=:green)
      end
    end

    img
end

circle(x::T, y::T, r::Real) where {T<:Real} = sqrt(x^2+y^2)-r
ellipse(x, y, a, b) = sqrt((b*x)^2+(a*y)^2) - a*b

dealunay_img(N=400) #, f = (x,y)->ellipse(x,y,1.0,0.5))
