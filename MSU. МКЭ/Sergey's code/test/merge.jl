include("../src/Dealunay.jl")
using .Dealunay

using Plots; pyplot()

function main()
    v = [(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0), (1.0, 1.0)]
    v = sort(v, lt=pointIsLess)
    p = [Dealunay.Point(v[j]..., j) for j in 1:length(v)]

    L, R = 1, length(p)+1
    mid = (L+R)÷2
    insert(L, L+1, p)
    # insert(mid, mid+1, p)
    triangle(mid, mid+1, mid+2, p)


    BT = hull((L,mid-1), (mid,R-1), p)
    UT = reverse(hull((R-1,mid), (mid-1,L), p))
    mergeHulls(BT, UT, p)


    img = plot(legend=false)
    for point in p
      for n in point.list
          plot!([point.x, p[n].x], [point.y, p[n].y], lc=:black)
      end
    end
    for point in v
      scatter!([point[1]], [point[2]], label="", c=:green)
    end
    img
end

main()
