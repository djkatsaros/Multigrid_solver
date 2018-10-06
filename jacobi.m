
function  u = jacobi(A,f, u0, nit)


    u = u0;


        D = diag(diag(A));
        B = A - D;
        Di = 1 ./ diag(A);

    for i = 1:nit
        u = 2/3* ( Di .* ( f - B *u) ) + 1/3 *u;
    end



end
