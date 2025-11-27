const Qb = algebraic_closure(QQ);

function lattice_positive(Lf::ZZLatWithIsom, h::Union{Vector, Nothing} = nothing)::Tuple{Bool, Any}
    f = isometry(Lf)
    L = lattice(Lf)
    tau = get_tau(f)
    bi_form = get_bilinealform(L)

    # step 1
    C0 = get_C0(Lf, tau)

    # step 2 - Check if C0 has obstructing roots => positive
    if !isone(C0)
        Cfancy = get_Cfancy(Lf, C0)
        if !isempty(Cfancy)
            return (false, QQ(Cfancy[1][1]))
        end
    end

    # step 3 - Prepare eigenvectors from tau and tau inverse 
    v = get_eigenvector(f, tau)
    w = get_eigenvector(f, tau^(-1))

    if bi_form(v,w)<0
        v = -v
    end

    # step 4 - Get the first h value if there is no in arguments
    if (h==nothing) 
        h = get_h(L,v,w)
    end
    # step 5 - Get the R set based on current h value
    Rh = get_R(L, h, bi_form)
    # step 6 - Check all of the entries of R if there exists obstructing root => positive
    for r in Rh 
        result = check_R(r, v, w, bi_form) 
        if !result[1] return result end
    end
    # step 7
    A = get_A(h, f, bi_form)
    # step 8 - Calculate new set of h values, that can show us obstructing roots
    H = map((a,b)-> -b*h + a*(f*h), A)
    # step 9 - Check all h from H to check if there exists any obstructing root => positive
    for h in H
        Rh = get_R(h, L, bi_form)
        for r in Rh
            result = check_R(r, v, w, bi_form)
            if !result[1] return result end
        end
    end
    return (true, QQ(0))
end

function get_tau(f::QQMatrix) ::QQBarFieldElem
    tau = QQ(0)
    for lambda in eigenvalues(Qb, f)
        if abs(lambda)>tau
            tau = abs(lambda)
        end
    end
    return tau
end

function get_bilinealform(L::ZZLat)
    return (a,b)-> ((change_base_ring(Qb, a))*change_base_ring(Qb, gram_matrix(ambient_space(L)))*transpose(change_base_ring(Qb, b)))[1]
end

function get_eigenvector(f::QQMatrix, lambda::QQBarFieldElem)
    return eigenspace(change_base_ring(Qb, f), lambda; side=:left)
end

function get_C0(Lf::ZZLatWithIsom, tau::QQBarFieldElem)::PolyRingElem
    charPolyF = characteristic_polynomial(Lf)
    x = gen(parent(charPolyF))
    (n, remainder) = remove(charPolyF, x-1) # remove all x-1 factors
    return div(remainder, minpoly(parent(charPolyF), tau)) #remove Salem polynomial of salem number tau
end

function get_Cfancy(Lf::ZZLatWithIsom, C0)::Array{Vector}
    return short_vectors(lattice(kernel_lattice(Lf, C0)), 2 , 2)
end

function get_h(L::ZZLat, v, w)::QQMatrix
    RR = ArbField(64)
    l = number_of_rows(basis_matrix(L))
    h0 = (v+w)*change_base_ring(Qb, inv(basis_matrix(L)))
    n = 10
    z = matrix(Qb,1,l,rand(-20:20, l)) #vector rand(-10:10, l) is in lattice basis => z vector is in ambient space basis
    h = (z+h0*n) #in lattice basis
    h = map(x->QQ(ZZ(floor(RR(x)))) , h)*basis_matrix(L) # in ambient space basis and rounded
    while bi_form(h,h) <= 0 || n == 10000
        n+=1
        h = (z+h0*n) #in lattice basis
        h = map(x->QQ(ZZ(floor(RR(x)))) , h)*basis_matrix(L) # in ambient space basis and rounded
    end
    #z = ones(Int64,number_of_rows(basis_matrix(L)), 1) # for test purposes
    #z = matrix(Qb, transpose(z))
    print(n)
    return h
end

function get_R(L::ZZLat, h::QQMatrix)::Vector{QQMatrix}
    return short_vectors_affine(L,h,0,-2)
end

function process_finite_sets_of_h(h, f::QQMatrix, bi_form, L::ZZLat)
    x = bi_form(h, h)
    y = bi_form(h, h*f)
    z = y^2-x^2
    if (z<0) return [] end #then discrininant for a will be <0 
    RR = ArbField(64)
    b_min = trunc(convert(Float64,RR(-sqrt(2*z/x))))
    for b = b_min:-1
        a_max = trunc(convert(Float64,RR((ZZ(b)*y+sqrt(z*(ZZ(b)^2+2*x)))/x)))
        #print(a_max)
        for a = 1:a_max
            h_new = -ZZ(b)*h + ZZ(a)*(h*f)
            Rh = get_R(L, h_new)
            for r in Rh
                result = check_R(r, v, w, bi_form)
                if !result[1] 
                    return result
                end
            end
        end
    end
    return (true, zero(h))
end

function get_A(h, f::QQMatrix, bi_form) # need to combine it with new h calculations, so no set A or H should be created(memory optimization)
    #::Array{Tuple{Int, Int}}
    x = bi_form(h, h)
    y = bi_form(h, h*f)
    z = y^2-x^2
    if (z<0) return [] end #then discrininant for a will be <0 
    RR = ArbField(64)
    A = []
    b_min = trunc(convert(Float64,RR(-sqrt(2*z/x))))
    #print(bmin)
    for b = b_min:-1
        a_max = (ZZ(b)*y+sqrt(z*(ZZ(b)^2+2*x)))/x
        print(a_max)
        #for a = 1:a_max
        #    push!(A,(ZZ(a),ZZ(b)))
        #end
    end
    return bmin
end

function check_R(r, v, w, bi_form) :: Tuple{Bool, Any}
    if bi_form(r, v)*bi_form(r, w) < 0 return (false, r)
    else return (true, zero(r)) end
end