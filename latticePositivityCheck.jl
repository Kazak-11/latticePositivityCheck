function latticePositive(Lf::ZZLatWithIsom, h) :: [Bool, QQFieldElem]
    f = isometry(Lf)
    L = lattice(Lf)
    tau = getTau(f)
    result = [true, 0] # initialize as positive without roots

    # step 1
    getC0(Lf)

    # step 2
    if C0 != 1
        Cfancy = getCfancy(Lf, C0)
        if Cfancy != empty
            return [false, Cfancy[0]]
        end
    end

    # step 3
    v = getEigenvector(f, tau)
    w = getEigenvector(f, tau^(-1))

    if dot(v,w)<0
        v = -v
    end

    # step 4

    # step 5
    Rh = getR(h)
    # step 6
    for r in Rh 
        result = checkR(r) 
        if !result[0] return result end
    end
    # step 7
    A = getA(h, f)
    # step 8
    #H = foreach
    # step 9
    for h in H
        Rh = getR(h)
        for r in Rh
            result = checkR(r)
            if !result[0] return result end
        end
    end
    return result
end

# (true,0)

function getTau(f)
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

function getEigenvector(f, lambda)
    # should I use OSCAR functionality or Julia to get eigenvectors?
    return solve(f-lambda*identity_matrix(f), zero(f),side == :right)
end

function getC0(Lf)
    #charF

    let charPolyF = characteristic_polynomial(Lf)

    # how to divide charPolyF with (x-1)? 
end

function getCfancy(Lf, C0)
    kernel_lattice(Lf, C0) 
    # how to get roots from the resulting lattice?
    #https://docs.oscar-system.org/dev/NumberTheory/QuadFormAndIsom/latwithisom/#kernel_lattice-Tuple%7BZZLatWithIsom,%20Integer%7D
end

function getR(h)
end

function checkR(r)
    if dot(r, v)*dot(r, w) < 0 return [false, r]
    else return [true, 0] end
end