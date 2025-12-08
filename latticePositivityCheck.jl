function _get_tau(f::QQMatrix, Qb) ::QQBarFieldElem
  tau = QQ(0)
  for lambda in eigenvalues(Qb, f)
    if abs(lambda)>tau
      tau = abs(lambda)
    end
  end
  return tau
end

function _get_bilinealform(L::ZZLat, Qb)
  #return (a, b)->inner_product(ambient_space(L), a,b)[1,1] -> it doesn't go well as there is no convienient way to convert ZZLat in QQ field to the Qb field
  return (a,b)-> ((change_base_ring(Qb, a))*change_base_ring(Qb, gram_matrix(ambient_space(L)))*transpose(change_base_ring(Qb, b)))[1]
end

function _get_eigenvector(f::QQMatrix, lambda::QQBarFieldElem)
  return eigenspace(f, lambda; side=:left)
end

function _get_C0(Lf::ZZLatWithIsom, tau::QQBarFieldElem)::PolyRingElem
  charPolyF = characteristic_polynomial(Lf)
  x = gen(parent(charPolyF))
  (n, remainder) = remove(charPolyF, x-1) # remove all x-1 factors
  return div(remainder, minpoly(parent(charPolyF), tau)) #remove Salem polynomial of salem number tau
end

function _get_Cfancy(Lf::ZZLatWithIsom, C0)
  return short_vectors_iterator(lattice(kernel_lattice(Lf, C0)), 2 , 2)
end

function _get_h(L::ZZLat, v, w, Qb, bi_form)::QQMatrix
  RR = ArbField(64)
  l = number_of_rows(basis_matrix(L))
  h0 = (v+w)*change_base_ring(Qb, inv(basis_matrix(L)))
  n = 1
  i = 1
  z = matrix(Qb,1,l,rand(-10:10, l)) # random small vector in lattice basis
  h = (z+h0*n) # in lattice basis
  h = map(x->QQ(ZZ(floor(RR(x)))) , h)*basis_matrix(L) # in ambient space basis and rounded
  while bi_form(h,h) <= 0 || n == 10000  # block on 10000, because for a bigger n time consumption is way too big
    i+=1
    if i == 100
      n+=1
      i = 1
    end
    h = (z+h0*n) #in lattice basis
    h = map(x->QQ(ZZ(floor(RR(x)))) , h)*basis_matrix(L) # in ambient space basis and rounded
  end
  return h
end

function _get_R(L::ZZLat, h::QQMatrix)::Vector{QQMatrix}
  return short_vectors_affine(L,h,0,-2)
end

@doc raw"""
  _process_finite_sets_of_h(h, f::QQMatrix, v, w, bi_form, L::ZZLat) -> Tuple{Bool, QQMatrix}

The function calculates all $r$, that can be obstructing roots of 'L'.
Calculation is made one by one and then the $r$ is checked
$r$ is based on pairs of integer $(a,b)$ of $-2x^2+2y^2+2aby>=x(a^2+b^2)$ with $a>0$, $b<0$, $x>0$,$y>0$
See Steps 7,8,9 of Algorithm 5.8 of "MINIMUM POSITIVE ENTROPY OF COMPLEX
ENRIQUES SURFACE AUTOMORPHISMS" by KEIJI OGUISO and XUN YU in Duke Mathematical Journal, 2020,
volume 169, number 18, pages 3565 - 3606

Return a tuple of a boolean that represents if isometry f of the lattice is positive and QQMatrix,
that represents an obstructing root
"""
function _process_finite_sets_of_h(h, f::QQMatrix, v, w, bi_form, L::ZZLat):: Tuple{Bool, QQMatrix}
  x = bi_form(h, h)
  y = bi_form(h, h*f)
  z = y^2-x^2
  if (z<0) return (true, zero(h)) end #then discrininant for a will be <0 and there is no roots (a,b)=> no obstructing roots => positive
  #trunc(ZZRingElem,-sqrt(2*z/x)) #check issues on github as round works in the same context
  b_min = round(ZZRingElem, -sqrt(2*z/x), RoundUp)
  for b = b_min:-1
    a_max = round(ZZRingElem, (ZZ(b)*y+sqrt(z*(ZZ(b)^2+2*x)))/x, RoundDown)
    for a = 1:a_max
      h_new = -ZZ(b)*h + ZZ(a)*(h*f)
      Rh = _get_R(L, h_new)
      for r in Rh
        result = _check_R(r, v, w, bi_form)
        if !result[1] 
          return result
        end
      end
    end
  end
  return (true, zero(h))
end

function _check_R(r, v, w, bi_form) :: Tuple{Bool, QQMatrix}
  if bi_form(r, v)*bi_form(r, w) < 0 return (false, r)
  else return (true, zero(r)) end
end

@doc raw"""
  isometry_is_positive(Lf::ZZLatWithIsom, h::Union{QQMatrix, Nothing} = nothing) -> Tuple{Bool, QQMatrix}

Given lattice with isometry 'Lf' and given 'h' vector in this lattice
according to the ambient space basis, such that $h^2>0$. 
If no such vector is given, 'h' will be calculated by function inside based on a random vector in 'Lf'.

Return a tuple of a boolean that represents if isometry f of the lattice is positive and QQMatrix,
that represents an obstructing root (see Algorithm 5.8 of "MINIMUM POSITIVE ENTROPY OF COMPLEX
ENRIQUES SURFACE AUTOMORPHISMS" by KEIJI OGUISO and XUN YU in Duke Mathematical Journal, 2020,
volume 169, number 18, pages 3565 - 3606)

For positive isometries it returns a vector of zeros with the same dimension as 'h'.
"""

function isometry_is_positive(Lf::ZZLatWithIsom, h::Union{QQMatrix, Nothing} = nothing):: Tuple{Bool, QQMatrix}
  Qb = algebraic_closure(QQ);
  f = ambient_isometry(Lf)
  L = lattice(Lf)
  tau = _get_tau(f, Qb)
  bi_form = _get_bilinealform(L, Qb)

  # step 1
  C0 = _get_C0(Lf, tau)
  # step 2 - Check if C0 has obstructing roots => positive
  if !isone(C0)
    Cfancy = _get_Cfancy(Lf, C0)
    first_element = iterate(Cfancy)
    if first_element !== nothing
      return (false, first_element[1][1])
    end
  end
  # step 3 - Prepare eigenvectors from tau and tau inverse 
  v = _get_eigenvector(f, tau)
  w = _get_eigenvector(f, tau^(-1))

  if bi_form(v,w)<0
    v = -v
  end
  # step 4 - Get the first h value if there is no in arguments
  if (h===nothing)
    h = _get_h(L,v,w, Qb, bi_form)
  end
  # step 5 - Get the R set based on current h value
  Rh = _get_R(L, h)
  # step 6 - Check all of the entries of R if there exists obstructing root => positive
  for r in Rh
    result = _check_R(r, v, w, bi_form) 
    if !result[1] return result end
  end
  # step 6,7,8 are combined to process them iteratively
  return _process_finite_sets_of_h(h, f, v, w, bi_form, L)
end