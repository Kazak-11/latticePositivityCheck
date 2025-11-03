const n = 1000000

function lattice_positive(Lf::ZZLatWithIsom, h::Vector) :: (Bool, QQFieldElem)
    f = isometry(Lf)
    L = lattice(Lf)
    tau = get_tau(f)

    # step 1
    C0 = get_C0(Lf)

    # step 2
    if C0 != 1
        Cfancy = get_Cfancy(Lf, C0)
        if Cfancy != empty
            return (false, Cfancy[0])
        end
    end

    # step 3
    v = get_eigenvector(f, tau)
    w = get_eigenvector(f, tau^(-1))

    if dot(v,w)<0
        v = -v
    end

    # step 4
    h = get_h(L,v,w)
    # step 5
    Rh = get_R(h, L)
    # step 6
    for r in Rh 
        result = check_R(r, v, w) 
        if !result[0] return result end
    end
    # step 7
    A = get_A(h, f)
    # step 8
    H = map((a,b)-> -b*h + a*(f*h), A)
    # step 9
    for h in H
        Rh = get_R(h, L)
        for r in Rh
            result = check_R(r, v, w)
            if !result[0] return result end
        end
    end
    return (true, 0)
end

function get_tau(f)
    Qb = algebraic_closure(QQ);
    tau = QQ(0)
    for lambda in eigenvalues(Qb, f)
        if abs(lambda)>tau
            tau = abs(lambda)
        end
        return nothing
    end
    return tau
end

function get_eigenvector(f, lambda)
    # should I use OSCAR functionality or Julia to get eigenvectors?
    return solve(f-lambda*identity_matrix(f), zero(f),side == :right)
end

function get_C0(Lf)
    charPolyF = characteristic_polynomial(Lf)

    # how to divide charPolyF with (x-1)? 
end

function get_Cfancy(Lf, C0)
    # is it correct way to use  short vectors?
    return map(short_vectors(lattice(kernel_lattice(Lf, C0)), 1.4 , 1.5)) do (v,n)
        if dot(v,v) == -2.0 return v
        else return nothing
        end
    end
end

function get_h(L, v, w)
    z = rand(short_vectors(L,0,2))
    return map(x->floor(x) , (z+n*(v+w)))
end

function get_R(h, L)
    return map(short_vectors(L, 1.4 , 1.5)) do (v,n)
        if dot(v,v) == -2.0 && dot(v,h) == 0.0 return v
        else return nothing
        end
    end
end

function get_A(h, f)
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

function check_R(r, v, w)
    if dot(r, v)*dot(r, w) < 0 return (false, r)
    else return (true, 0) end
end