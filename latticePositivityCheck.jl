const n = 1000000

function lattice_positive(Lf::ZZLatWithIsom, h::Union{Vector, nothing} = nothing)::(Bool, QQFieldElem)
    f = isometry(Lf)
    L = lattice(Lf)
    tau = get_tau(f)

    # step 1
    C0 = get_C0(Lf, tau)

    # step 2 - Check if C0 has obstructing roots => positive
    if !isone(C0)
        Cfancy = get_Cfancy(Lf, C0)
        if !isempty(Cfancy)
            return (false, QQ(Cfancy[0]))
        end
    end

    # step 3 - Prepare eigenvectors from tau and tau inverse 
    v = get_eigenvector(f, tau)
    w = get_eigenvector(f, tau^(-1))

    if dot(v,w)<0
        v = -v
    end

    # step 4 - Get the first h value if there is no in arguments
    if (h==nothing) 
        h = get_h(L,v,w)
    end
    # step 5 - Get the R set based on current h value
    Rh = get_R(L, h)
    # step 6 - Check all of the entries of R if there exists obstructing root => positive
    for r in Rh 
        result = check_R(r, v, w) 
        if !result[0] return result end
    end
    # step 7
    A = get_A(h, f)
    # step 8 - Calculate new set of h values, that can show us obstructing roots
    H = map((a,b)-> -b*h + a*(f*h), A)
    # step 9 - Check all h from H to check if there exists any obstructing root => positive
    for h in H
        Rh = get_R(h, L)
        for r in Rh
            result = check_R(r, v, w)
            if !result[0] return result end
        end
    end
    return (true, QQ(0))
end

function get_tau(f::QQMatrix)::QQFieldElem
    Qb = algebraic_closure(QQ);
    tau = QQ(0)
    for lambda in eigenvalues(Qb, f)
        if abs(lambda)>tau
            tau = abs(lambda)
        end
    end
    return tau
end

function get_eigenvector(f::QQMatrix, lambda::QQFieldElem)::Vector
    # should I use OSCAR solve functionality or Julia function to get eigenvectors?
    return solve(f-lambda*identity_matrix(f), zero(f),side == :right)
end

function get_C0(Lf::ZZLatWithIsom, tau::QQFieldElem)::PolyRingElem
    charPolyF = characteristic_polynomial(Lf)
    while divides(charPolyF, (x-1))
        div!(charPolyF, (x-1))
    end
    return div!(charPolyF, minpoly(QQ, tau)) #remove Salem polynomial of salem number tau
end

function get_Cfancy(Lf::ZZLatWithIsom, C0)::Array{Vector}
    # is it correct way to use short vectors?
    return map(short_vectors(lattice(kernel_lattice(Lf, C0)), 1.4 , 1.5)) do (v,n)
        if dot(v,v) == -2.0 return v
        else return nothing
        end
    end
end

function get_h(L::ZZLat, v::Vector, w::Vector)::Vector
    z = rand(short_vectors(L,0,2))
    return map(x->floor(x) , (z+n*(v+w)))
end

function get_R(L::ZZLat, h::Vector)::Array{Vector}
    return map(short_vectors(L, 1.4 , 1.5)) do (v,n)
        if dot(v,v) == -2.0 && dot(v,h) == 0.0 return v  # should I use ≈
        else return nothing
        end
    end
end

function get_A(h::Vector, f::QQMatrix)::Array{(Int, Int)}
    x = dot(h, h)
    y = dot(h, f*h)
    A = []

    bmin = trunc(-√(2(y^2-x^2)/x))
    for b = bmin:-1
        D = √((y^2-x^2)*(b^2+2x))
        amin = trunc(Int, (b*y-D)/x)
        amax = trunc(Int, (b*y+D)/x)
        for a = amin:amax
            push!(A,(a,b))
        end
    end
    return A
end

function check_R(r::Vector, v::Vector, w::Vector) :: (Bool, QQFieldElem)
    if dot(r, v)*dot(r, w) < 0 return (false, QQ(r))
    else return (true, QQ(0)) end
end