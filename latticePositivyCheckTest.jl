function test()
    B = matrix(QQ, 3, 5, [  1 0 0 0 0;
                            0 0 1 0 1;
                            0 0 0 1 0])

    G = matrix(QQ, 5, 5, [  2 -1  0  0  0;
                            -1  2 -1  0  0;
                            0 -1  2 -1  0;
                            0  0 -1  2 -1;
                            0  0  0 -1  2])

    f = matrix(QQ, 5, 5, [ 1  0  0  0  0;
                             -1 -1 -1 -1 -1;
                              0  0  0  0  1;
                              0  0  0  1  0;
                              0  0  1  0  0])
    L =  integer_lattice(B; gram = G)
    Lf = integer_lattice_with_isometry(L, f)

    tau = get_tau(f)
    v = get_eigenvector(f, tau)
    
end