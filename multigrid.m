% Alternative numerical solver  - mulitscale in  nature.
%
% Require the sparse matrix A
% Jacobi Iteration relaxes high freq. modes in the error.
%       e^(0)(mix of high and low energies) -- Jacobi -->e^(n)
%  Smooth multigrid uses coarse and fine meshes
%
%   \nu cycle for A_n u_n = f
%   compute r_n = f -A_n u_n
%   Transfer to coarse grid
%       r_2n = I_nr_n
%       solve A _2n u_corr = r_2n
%   add u <-- u + I_n u_corr
%
%
%
%  A = D +B where D are diagonal, and B are off-diagonal terms.
% Compare Jacobi and multigrid performance.
%


% Require the sparse matrix A

%  A = D +B where D are diagonal, and B are off-diagonal terms.
function  u = multigrid(A,f, u0, nit)

    %  Dirichlet problem paramters.
    M = size(f,1);
    N = sqrt(M)+1;


    if N < 10
        u = A\f; %umfpack(A, "\", f);
    else
        %some steps of weighted Jacobi
        u = u0;

        D = diag(diag(A));
        B = A - D;
        %The inverse.
        Di = 1 ./ diag(A);

        for i = 1:nit/2
            u = 2/3* ( Di .* ( f - B*u) ) + 1/3 *u;
        end
        HM = (N/2-1)^2;
        ind1 = (1:HM)';
        ind2 = 2*ind1 - 1 + N* ceil( ind1/ (N/2-1) );
        Ih = sparse([ind1;ind1],[ ind2-1; ind2+1], [1/8*ones(2*HM,1)], HM, M) ...
        +  sparse([ind1;ind1], [ind2+N-1; ind2-N+1], [1/8*ones(2*HM,1)], HM, M) ...
        +  sparse([ind1;ind1], [ind2+N-2; ind2-N], [1/16*ones(2*HM,1)], HM, M) ...
        +  sparse([ind1;ind1], [ind2+N; ind2-N+2], [1/16*ones(2*HM,1)], HM, M) ...
        +  sparse(ind1, ind2, [1/4*ones(HM,1)], HM, M);
       
        Iht = 4*Ih';

        r = Ih*(f-A*u);


        ucorr = multigrid(Ih*A*Iht, r, 0*r, nit);
        u = u + Iht*ucorr;

        for i = 1:nit/2
            u = 2/3* (Di .* ( f - B *u)) + 1/3 *u;
        end
    end

end
