function latticePositive(Lf::ZZLatWithIsom, h) 
    f = isometry(Lf)
    L = lattice(Lf)
    tau = getTau(f)

    # step 1
    getC0(f)

    # step 2
    if C0 != 1
        Cfancy = getCfancy()
        if Cfancy != empty
            return (Cfancy[0])
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
    result = foreach(checkR, Rh) 
    # step 7
    A = getA(h, f)
    # step 8
    H = getH(h, f, A)
    # step 9
    result = foreach(function (h)
        Rh = getR(h)
        result = foreach(checkR, Rh)
        return result
    end, H)
    return result
end

function getTau(f)
end

function getCfancy()
end


function getR(h)
end

function checkR(r)
    if (r, v)(r, w) < 0
        return r
    end
end
