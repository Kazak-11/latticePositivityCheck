function test()
  #Qb = algebraic_closure(QQ);
  R, x = QQ[:x]
  # Gram matrix
  G = QQ[0 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 1 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 -2 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1; 0 0 0 -2 0 0 0 0 0 0 1 0 1 0 1 1 1 0; 0 0 0 0 -2 0 0 0 0 0 0 0 1 1 1 0 0 1; 0 0 0 0 0 -2 0 0 0 0 1 1 0 1 0 1 0 0; 0 0 0 0 0 0 -2 0 0 0 1 0 0 1 1 0 1 1; 0 0 0 0 0 0 0 -2 0 0 1 1 1 0 0 0 1 0; 0 0 0 0 0 0 0 0 -2 0 0 1 0 0 0 0 1 1; 0 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 1 0 0; 1 0 0 1 0 1 1 1 0 0 -2 2 2 2 2 0 0 2; 1 0 1 0 0 1 0 1 1 0 2 -2 1 0 2 2 1 1; 1 0 1 1 1 0 0 1 0 0 2 1 -2 0 0 1 2 1; 1 0 1 0 1 1 1 0 0 0 2 0 0 -2 1 2 2 1; 1 0 0 1 1 0 1 0 0 1 2 2 0 1 -2 1 1 1; 1 0 1 1 0 1 0 0 0 1 0 2 1 2 1 -2 2 2; 1 0 0 1 0 0 1 1 1 0 0 1 2 2 1 2 -2 1; 1 0 1 0 1 0 1 0 1 0 2 1 1 1 1 2 1 -2]

  # isometry
  f = QQ[3 1 0 -1 0 1 1 -1 -1 -1 1 0 -1 1 0 0 0 0; 2 1 0 -1 0 0 0 -1 0 0 0 0 -1 1 0 0 0 0; 2 1 0 -1 0 1 1 -1 0 -1 1 0 -1 1 0 0 0 0; 2 1 0 -1 0 1 1 -1 -1 0 1 0 -1 1 0 0 0 0; 3 1 0 -1 0 1 0 -1 -1 -1 1 0 -1 1 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 1 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 3 1 0 -1 0 0 1 -1 -1 -1 1 0 -1 1 0 0 0 0; 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 11 5 -3//2 -5//2 -1//2 1//2 1//2 -7//2 -3//2 -3//2 0 -1 -3 2 0 0 0 0; 6 2 0 -3 0 3 2 -2 -2 -2 2 1 -3 3 0 0 0 0; 0 -1 1//2 -1//2 1//2 3//2 3//2 -1//2 -1//2 -1//2 1 0 -1 2 0 0 0 0; 2 0 0 -1 0 2 2 -1 -1 -1 2 1 -2 2 0 0 0 0; 3 1 0 -1 0 1 2 -1 -1 -2 2 0 -1 2 0 -1 0 0; 6 2 -1 -2 0 1 2 -3 -2 -1 1 -1 -3 3 1 0 0 0; 14 7 -3//2 -5//2 -3//2 -1//2 -3//2 -3//2 -3//2 -5//2 0 0 -1 0 -1 -1 0 -1; 2 1 1 -1 0 1 1 -1 -1 -1 0 0 -1 2 -1 1 1 0]

  # check that f is an isometry
  @assert f*G*transpose(f)  == G

  # define the lattice
  B = QQ[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1//2 1//2 1//2 1//2 1//2 1//2 1//2 1//2 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]

  V = quadratic_space(QQ, G)
  Vf = quadratic_space_with_isometry(V, f)
  L = lattice(Vf, B)
  @assert lattice(Vf,basis_matrix(L)*f) == L

  tau = get_tau(f)

  C0 = x^14 + 2*x^13 + 5*x^12 + 8*x^11 + 11*x^10 + 14*x^9 + 15*x^8 + 16*x^7 + 15*x^6 + 14*x^5 + 11*x^4 + 8*x^3 + 5*x^2 + 2*x + 1

  bi_form = get_bilinealform(lattice(L))

  v = get_eigenvector(f, tau)
  w = get_eigenvector(f, tau^(-1))
  h = get_h(lattice(L), v, w)
  @testset verbose = true "Lattice Positivity Checker" begin
    @testset "Vector Sanity Check" begin
      @test bi_form(v,v) == 0
      @test bi_form(w,w) == 0
    end
    @testset "C Sets Check" begin
      #@test get_C0(L, tau) == C0
      @test get_Cfancy(L, C0) == []
    end
    @testset "h Check" begin
      @test bi_form(h,h)>0
    end
    @testset "r from h Check" begin
      Rh = get_R(lattice(L), h)
      if Rh != []
        @test check_R(Rh[1], v, w, bi_form)[1] == true
      end
    end
    @testset "Whole Algorithm Check" begin
      #@test lattice_positive(L, h)[1] == true
    end
  end
    nothing
end
