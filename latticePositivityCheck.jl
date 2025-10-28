function latticePositive(Lf::ZZLatWithIsom, h) :: (Bool, QQFieldElem)
    f = isometry(Lf)
    L = lattice(Lf)
    tau = getTau(f)

    # step 1
    getC0(Lf)

    # step 2
    if C0 != 1
        Cfancy = getCfancy()
        if Cfancy != empty
            return (false, Cfancy[0])
        end
    end

    # step 3
    v = getEigenvector(f, tau)
    w = getEigenvector(f, tau^(-1))

    if (v,w)<0
        v = -v
    end

    # step 4

    # step 5
    Rh = getR(h)
    # step 6
    result = (0, false)
    result = foreach(checkR, Rh) 
    # step 7
    A = getA(h, f)
    # step 8
    H = getH(h, f, A)
    # step 9
    result = foreach(H) do h
        Rh = getR(h)
        result = foreach(checkR, Rh)
        return result
    end
    return result
end

# (true,0)

function getTau(f)
    Qb = algebraic_closure(QQ);
    tau = QQ(0)
    foreach(eigenvalues(Qb, f)) do lambda
        if abs(lambda)>tau
            tau = abs(lambda)
        end
        return nothing
    end
    return tau
end

function getEigenvector(f, lambda)
    # should I use OSCAR functionality or Julia to get eigenvectors?
    return solve(f-lambda*identity_matrix(f), zero(f),side == :right)
end

function getC0(Lf)
    #charF

    let charPolyF = characteristic_polynomial(Lf)
end

function getCfancy(Lf, C0)
    kernel_lattice(Lf, C0) 
    # how to get roots from the resulting lattice?
    #https://docs.oscar-system.org/dev/NumberTheory/QuadFormAndIsom/latwithisom/#kernel_lattice-Tuple%7BZZLatWithIsom,%20Integer%7D
end

function getR(h)
end

function checkR(r)
    if dot(r, v)*dot(r, w) < 0
        return r
    end
    return nothing
end
