function lattice_positive(Lf::ZZLatWithIsom, h::ZZFieldElem) :: (Bool, QQFieldElem)
    f = isometry(Lf)
    L = lattice(Lf)
    tau = get_tau(f)
    result = (true, 0) # initialize as positive without roots

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
    # is h integer as it is a part of L?
    # step 5
    Rh = get_R(h)
    # step 6
    for r in Rh 
        result = check_R(r) 
        if !result[0] return result end
    end
    # step 7
    A = get_A(h, f)
    # step 8
    H = map((a,b)-> -b*h + a*(f*h), A)
    # step 9
    for h in H
        Rh = get_R(h)
        for r in Rh
            result = check_R(r)
            if !result[0] return result end
        end
    end
    return result
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
    #charF

    charPolyF = characteristic_polynomial(Lf)

    # how to divide charPolyF with (x-1)? 
end

function get_Cfancy(Lf, C0)
    # is it correct way to use 
    return map(short_vectors(lattice(kernel_lattice(Lf, C0)), 1.4 , 1.5)) do (v,n)
        if dot(v,v) == -2.0 return v
        else return nothing
        end
    end
end

function get_R(h)
    return map(short_vectors(lattice(kernel_lattice(Lf, C0)), 1.4 , 1.5)) do (v,n)
        if dot(v,v) == -2.0 && dot(v,h) == 0.0 return v
        else return nothing
        end
    end
# each v in L, s.t. dot(h,v)=0 and dot(v,v)=-2
end

function check_R(r)
    if dot(r, v)*dot(r, w) < 0 return (false, r)
    else return (true, 0) end
end