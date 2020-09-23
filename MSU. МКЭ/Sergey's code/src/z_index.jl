# sources:
#   http://asgerhoedt.dk/?p=276
#   http://www-graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN

B = [0x0000ffff0000ffff,
     0x00ff00ff00ff00ff,
     0x0f0f0f0f0f0f0f0f,
     0x3333333333333333,
     0x5555555555555555]
S = [16, 8, 4, 2, 1]

function separate(a::Int)::Int
    a &= 0xffffffff
    for i=1:5
        a = (a | (a << S[i])) & B[i]
    end
    return a
end

function mortonCode(x::Int, y::Int)
    separate(x) | (separate(y) << 1)
end

# function comb2(a::UInt, i=1)::Int
#     i > 5 && return a
#     comb2((a | (a << S[i])) & B[i], i+1)
# end
#
# function morton2(x::Int, y::Int)
#     comb2(UInt(x)&0xffffffff) | (comb2(UInt(y)&0xffffffff) << 1)
# end

function compact(x::Int)::Int
    x &= 0x5555555555555555;
    x = (x | (x >>  1)) & 0x3333333333333333
    x = (x | (x >>  2)) & 0x0f0f0f0f0f0f0f0f
    x = (x | (x >>  4)) & 0x00ff00ff00ff00ff
    x = (x | (x >>  8)) & 0x0000ffff0000ffff
    x = (x | (x >>  16))& 0x00000000ffffffff
    return x
end

function mortonDecode(c::Int)#::NTuple(2,Int)
    x = compact(c)
    y = compact(c >> 1)
    (x, y)
end

##
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

##
mutable struct Cell{T<:Real}
    active::Bool
    leaf::Bool
    h::T
end

struct Grid2d{T<:Real} <: AbstractArray{Cell{T}, 1}
    data::Vector{Cell{T}}
    activemap::BitVector
    leafmap::BitVector
    level::Int
end

function Grid2d{T}(level::Int, h::T) where T
    Grid2d{T}([Cell(false,false,h/2^level) for j=1:4^level],
              BitArray(zeros(4^level)), BitArray(zeros(4^level)), level)
end

# тип по-умолчанию Float64
Grid2d(level::Int, h::T) where T = Grid2d{Float64}(level, h)

Base.IndexStyle(::Type{<:Grid2d{T}}) where T = IndexLinear()

Base.parent(g::Grid2d{T}) where T = g.data

Base.size(g::Grid2d{T}) where T = size(parent(g))

Base.getindex(g::Grid2d{T}, i::Integer) where T = g.data[i]

Base.setindex!(g::Grid2d{T}, a, i::Integer) where T = (g.data[i] = a)

##
circle(x::T, y::T, r::T) where {T<:Real} = x^2+y^2-r^2
ring(x::T, y::T) where {T<:Real} = circle(x,y,T(1))*circle(x,y,T(0.5))
ring2(x::T, y::T) where {T<:Real} = min(ring(x,y),ring(x-1.5,y))
ring3(x::T, y::T) where {T<:Real} = min(ring(x,y),ring(x-1.5,y),ring(x-0.75,y-0.75*sqrt(3)))
rectangle(x::T, y::T) where {T<:Real} = abs(x-y) + abs(x+y) - 2.0

function vertices(ind::NTuple{2,Int}, h::T, xmin::T, ymin::T) where {T<:Real}
    (xmin.+h*(ind[1].+[0,0,1,1]),
     ymin.+h*(ind[2].+[0,1,1,0]))
end

function check(ind::NTuple{2,Int}, h::T, xmin::T, ymin::T, f)::Bool where {T<:Real}
    shifts = [(0,0),(0,1),(1,1),(1,0)]
    subsquares = [vertices(2 .* ind .+ shift, h/2, xmin, ymin) for shift in shifts]
    res = false
    for square in subsquares
        vals = map(f, square...)
        res |= minimum(vals.*circshift(vals,1)) <= 0
    end
    res
end

function top_down_traverse(grid::Grid2d{T}, D::Domain2D{<:Real}, f) where {T<:Real}
    for j = 1:length(grid)
        ind2d = mortonDecode(j-1)
        grid[j].active = check(ind2d, grid[j].h, D.xmin, D.ymin, f)
        grid.activemap[j] = grid[j].active
    end
end

# добавляем листья в сетку down_grid, обновляя информацию на сетке up_grid
function down_up_traverse(down_grid::Grid2d{T}, up_grid::Grid2d{T},
        D::Domain2D{<:Real}; last_lvl::Bool=false) where {T<:Real}

    if last_lvl
        for j = 1:4:length(down_grid)
            if true in down_grid.activemap[j:j+3]
                down_grid.leafmap[j:j+3] .= true
                up_grid[((j-1) >> 2) + 1].leaf = true
            end
        end
        return
    end

    neighbor_fltr = Bool[1 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 1]

    j = 1
    while j <= length(down_grid)
        if down_grid[j].leaf
            ind = ((j-1) >> 2) << 2 + 1  # индекс начала блока
            for i = 0:3
                down_grid.leafmap[ind+i] = ~down_grid[ind+i].leaf
            end
            # индекс клетки в up_grid, соответствующей блоку из down_grid
            up_grid[((j-1) >> 2) + 1].leaf = true

            ind2d = mortonDecode(ind-1)
            neighbor = 1 .+ map(mortonCode,ind2d[1].+2*[0,-1,1,0],ind2d[2].+2*[-1,0,0,1])
            neighbor = neighbor[vec((!).(down_grid.leafmap[ind:ind+3])' * neighbor_fltr) .> 0]
            for n in filter(z->(1<=z<=length(down_grid)),neighbor)
                m = ((n-1) >> 2) << 2 + 1
                if reduce(|,map(a->a.leaf,down_grid[m:m+3])) == false
                    down_grid.leafmap[m:m+3] .= true
                    up_grid[((n-1) >> 2) + 1].leaf = true
                end
            end

            j = ind+4
            continue
        end
        j += 1
    end
end

##
using Plots
pyplot()

function main(;f,levels::Int=4)
    # D = Domain2D(xmin=-1.0, ymin=-1.0, xmax=2.5, ymax=1+0.75*sqrt(3))
    D = Domain2D(xmin=-2.0, ymin=-1.0, xmax=2.0, ymax=1.0)
    h = max(D.Lx, D.Ly)    # сторона квадрата
    grids = Grid2d.(1:levels,h)

    for grid in grids
        top_down_traverse(grid, D, f)
    end

    down_up_traverse(grids[end], grids[end-1], D, last_lvl=true)
    for j = length(grids)-1:-1:2
        down_up_traverse(grids[j], grids[j-1], D)
    end

    p = contour(D.xmin : D.Lx/200 : D.xmax,
                D.ymin : D.Ly/200 : D.ymax, f, border=:none, levels=[0,100],
                xlim=(D.xmin,D.xmax).*1.01, lw=0.5, ratio=:equal, #fill=true,
                ylim=(D.ymin,D.ymax).*1.01, lc=:black, legend=:none)

    colour = [:blue, :red, :yellow, :green, :brown, :blue, :cyan]
    for i=length(grids):-1:1, j=1:length(grids[i])÷2
        if grids[i].leafmap[j]
            x, y = vertices(mortonDecode(j-1), grids[i][j].h, D.xmin, D.ymin)
            plot!([x; x[1]], [y; y[1]], label="", lc=colour[i])
        end
    end
    p
end

function thing(x::T, y::T) where {T<:Real}
    p, q = 2.0, 5.0
    min(circle(x+1,y,1.0+1e-6), abs((x-1.0)/p-y/q) + abs((x-1.0)/p+y/q) - 2.0)*circle(x+1,y,0.33)
end

main(f=thing, levels=6)

# savefig("thing.png")
