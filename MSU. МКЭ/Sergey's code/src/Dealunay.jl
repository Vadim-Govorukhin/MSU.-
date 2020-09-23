module Dealunay

export Point, pointIsLess, dealunay, turn!
# для тестов:
export triangle, hull, mergeHulls, insert


mutable struct Point{T<:Real}
    x::T
    y::T
    ind::Int
    # первый элемент - индекс точки
    # остальные - индексы соседних точек, лежащие в порядке против часовой стрелки, если
    # смотреть из самой точки
    list::Vector{Int}
    Point(x::T, y::T, ind::Int) where {T<:Real} = new{T}(x, y, ind, Vector{Int}())
end

function Base.isless(v::Point{T}, u::Point{T}) where {T<:Real}
    v.x > u.x  && (return false)
    v.x < u.x  && (return true)
    v.x == u.x && (return v.y < u.y)
end

function succ(v::Int, u::Int, p::Vector)
    v_list = p[v].list
    ind = findfirst(isequal(u), v_list)
    ind != nothing ? v_list[mod1(ind+1,length(v_list))] : 0
end
function pred(v::Int, u::Int, p::Vector)
    v_list = p[v].list
    ind = findfirst(isequal(u), v_list)
    ind != nothing ? v_list[mod1(ind-1,length(v_list))] : 0
end


using LinearAlgebra

mutable struct ConvexHull
    leftmost::Int
    rightmost::Int
end

function ccwTurn(Ax::T, Ay::T, Bx::T, By::T, Cx::T, Cy::T) where {T<:Real}
    M = [Ax Ay 1;
         Bx By 1;
         Cx Cy 1]
    # println("det: $(det(M))")
    det(M) > 0
end

function ccwTurn(A::Int, B::Int, C::Int, p::Vector)
    a, b, c = p[A], p[B], p[C]
    ccwTurn(a.x, a.y, b.x, b.y, c.x, c.y)
end

function isLeftOf(Z::Int, line::NTuple{2,Int}, p::Vector)
    X, Y = line
    ccwTurn(X, Y, Z, p)
end

function isRightOf(Z::Int, line::NTuple{2,Int}, p::Vector)
    X, Y = line
    ccwTurn(Y, X, Z, p)
end

function hull(Vl::NTuple{2,Int}, Vr::NTuple{2,Int}, p::Vector, debug::Bool=false)
    X = Vl[2]
    Y = Vr[1]
    Z = p[Y].list[1]
    Z2 = pred(X, p[X].list[1], p)

    while true
        debug && println("hull")
        if isRightOf(Z, (X,Y), p)
            Z,  Y = succ(Z,  Y, p), Z
        elseif isRightOf(Z2, (X,Y), p)
            Z2, X = pred(Z2, X, p), Z2
        else
            return (X,Y)
        end
        # println("yola: $(isRightOf(Z, (X,Y), p)), $(isRightOf(Z2, (X,Y), p))")
    end
end

function insertBtoA(A::Int, B::Int, p::Vector)
    A_list = p[A].list
    B == A || B in A_list && return

    isempty(A_list) && (push!(A_list, B); return)

    A_len = length(A_list)
    next_i = 1
    ccw_check = [ccwTurn(a, A, B, p) for a in A_list]
    for i = 1:A_len
        left  = ccw_check[i]
        next_i = mod1(i+1,A_len)
        right = ccw_check[next_i]
        left || right && break
    end
    if next_i == 1 && isLeftOf(B, (A,A_list[end]), p)
            next_i = A_len + 1
    end
    insert!(A_list, next_i, B)
end

function insert(A::Int, B::Int, p::Vector)
    insertBtoA(A, B, p)
    insertBtoA(B, A, p)
end

function deleteBfromA(A::Int, B::Int, p::Vector)
    A_list = p[A].list
    # B in A_list &&
    deleteat!(A_list, findfirst(isequal(B), A_list))
end

function delete(A::Int, B::Int, p::Vector)
    deleteBfromA(A, B, p)
    deleteBfromA(B, A, p)
end

function outOfCircle(A::Int, B::Int, C::Int, D::Int, p::Vector)
    # println("A: $A, B: $B, C: $C, D: $D")
    D == 0 && (return true)
    D in [A,B,C] && (return true)
    a, b, c, d = p[A], p[B], p[C], p[D]
    M = [a.x  a.y  a.x^2+a.y^2  1;
         b.x  b.y  b.x^2+b.y^2  1;
         c.x  c.y  c.x^2+c.y^2  1;
         d.x  d.y  d.x^2+d.y^2  1]
    det(M) <= 0
end

function mergeHulls(BT::NTuple{2,Int}, UT::NTuple{2,Int}, p::Vector, debug::Bool=false)
# BT - lower common tangent
# UT - upper common tangent
    L, R = BT
    count = 1
    while L != R && BT != UT
        A = B = false
        insert(L, R, p)
        debug && println(p[L].list)
        debug && println(p[R].list)

# Delete the edges in DT(Vr) which are not Delaunay edges in DT(V)
# by determining if L is within the circumcircle of △(R, R1, R2).
# If so, the edge (R, R1) is not a Delaunay edge and must be deleted.
        R1 = pred(R, L, p)
        if isLeftOf(R1, (L,R), p)
            R2 = pred(R, R1, p)
            debug && println("R1: ", R1, "; L: ", L, "; R: ", R, "; R2: ", R2)
            while !outOfCircle(R1, L, R, R2, p)
                delete(R, R1, p)
                R1 = R2
                R2 = pred(R, R1, p)
                debug && println("hey")
            end
        else
            A = true
        end

# Same operation on the edges in DT(Vl).
        L1 = succ(L, R, p)
        if isRightOf(L1, (R, L), p)
            L2 = succ(L, L1, p)
            debug && println("L: $L, R: $R, L1:$L1, L2:$L2")
            while !outOfCircle(L, R, L1, L2, p)
                delete(L, L1, p)
                L1 = L2
                L2 = succ(L, L1, p)
                debug && println("ho")
            end
        else
            B = true
        end

debug && println("L: $L, R: $R, R1:$R1, L1:$L1")
        if A
            L = L1
        elseif B
            R = R1
        elseif outOfCircle(L, R, R1, L1, p)
            R = R1
        else
            L = L1
        end

        BT = (L, R)
        count += 1
        debug && count > length(p)/2 && println("merge")
        count > length(p) && break
    end

    insert(UT..., p)
end

function triangle(A::Int, B::Int, C::Int, p::Vector)
    insert(A, B, p)
    insert(A, C, p)
    insert(B, C, p)
    # [A, p[A].list[1], p[p[A].list[1]].list[1]]
end

##

function pointIsLess(v::Tuple{T,T}, u::Tuple{T,T}) where {T<:Real}
    v[1] > u[1]  && (return false)
    v[1] < u[1]  && (return true)
    v[1] == u[1] && (return v[2] < u[2])
end

##
function dealunay(L::Int, R::Int, p::Vector, debug::Bool=false)

    debug && println("L: $L, R: $R")
    L + 1 < R || println("shiiiiiiit")

    L + 3 == R && (triangle(L, L+1, L+2, p); return)
    L + 2 == R && (insert(L, L+1, p); return)

    mid = (L + R) ÷ 2
    dealunay(L, mid, p, debug)
    dealunay(mid, R, p, debug)

    if L < mid < R
        BT = hull((L,mid-1), (mid,R-1), p, true)
        UT = reverse(hull((R-1,mid), (mid-1,L), p, true))
        mergeHulls(BT, UT, p, debug)
    end
end


function turn!(p::Vector{<:Point}, θ::Real)
    iszero(θ) && return
    for i in axes(p,1)
        p[i].x, p[i].y = p[i].x*cos(θ)-p[i].y*sin(θ),
                         p[i].x*sin(θ)+p[i].y*cos(θ)
    end
end

end # module
