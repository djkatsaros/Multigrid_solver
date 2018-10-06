% Finite difference scheme for elliptic PDES
% PDE  is  \grad \cdot ( a(X) \grad u) = f
% discretize second derivative:
%       u'' = u(x_h) - 2u(x) + u(x-h)/ h^2 + O(h^2)  + ...
%       = (D^+ * D^-) u
%       Ex:, 1D, this is  -(a(x)u'(x))' = D^+(a(x-h/2)D^- u) = D^+(a(x-h/2) (u(x) -
%       u(x-h)/h) =
%           = a(x+h/2)((u(x+h) -u(x))/h^2 - a(x-h/2)((u(x)-u(x-h))/h^2 +
%           O(h^2) [//second order accurate , dropping the O(h^2) Terms//]
%
%       Higher D,becomes  - \grad \cdot (a(x) \grad u) = -\partial x_1 (a(x)
%       \partial x_1 u) - \partial x_2(a(x) \partial x_2 u) + O(h^2)
%
%       Stability etc.:
%       if PDE is Lu  = f, second order means that (L - L_h) u = O(h^2),
%       and have ||L_h^{-1}|| < p
%
%       THIS PROB: -(a(x/eps)u'(x))' = x^2
%       Test 1: let a(x/eps) = 2, then soln u = x^4/24
%
%       || u_{eps,h} - u_{hom}||_{L^2}  = ( sum of (squares*dx) )^{1/2}
%
%       2D MultiGrid  Num. Solver for AU = F (assume discretized, so linear
%       alg. problem)A \in \IR^n \times\IR^n
%           LU decomp ~ O(n^3) iterative schemes m(#iterations)* n^2 (m*n for
%           sparse sys.)
%           multigrid for sparse matrices makes it ~O(n) for sure.
%
%           Process
%           Take Jacobi Iteration: U^{n+1} = D^{-1}(f + (L+U)U^n) where
%           A = D - L - U, e^{n+1} = M e^n where M  = D^{-1}(L+U). Take
%           weighted sum of jacobi step;
%               M_\omega = (1- \omega)I + \omega(D^{-1}(L+U))
%               ~~Converges if and only if spectrum M is in the unit circle ~~
%
%               Model Prob - u'' = f on [0,1]
%           fintite diff. n intervals, kth eigenvector
%                   (v_k)j= sin(pi k j /  n) 0 \leq j \leq n
%                   lambda_k(A) = 4 sin^2 ( pi k /2n)
%                    lambda_k (M_\omega) = 1-2\omega sin^2(k pi /2n)
%                           max'd at k = 1, lambda_k(M_\omega)
%                                    \approx 1 - \omega/2 * pi^2/n^2
%                           take n=4, \omega = 1
%                           so \lambda_k(M_\oemga) \leq \approx 2/3, so
%                                      spec(M_\omega) \approx 2/3
%                    errors:
%                      e_0 = \sum_k a_k v_k
%                       M^je_0 = \sum_j a_k^jlambda_k^j v_k ~error goes to
%                           zero quickly for small e.val e.vect.s,
%                 iDeas
%                   1. Use Coarse Grid to get better initial guess - interpolate
%                 to a finer mesh, giving hgihly oscillatory errors, which
%                 can be treated well.
%                   2. compute residual, pose on a coarse grid to get a
%                   corrector
%
%                   TODO: solve -u'' = f with weighted Jacobi \omega = 2/3
%                           use u_0 = randn(n-1,1)



clf;

for j = 0:2
    eps = 10^(-2-j);

    y = @(x) (1.1 + sin(2 * pi * x ./ eps));%(2 +0*x); %
    %  There is a huge oscillation in the stiffness.
    s = @(x) x.^4/24 + x;
    %homogonized equation, bar{a} = 0.45
    h = @(x) (1/(12*0.45))*(x - x.^4);
    % y= 2 + 0*x;
    usv = ones(9,8);
    l2_err = ones(1,8);
    for i = 1:8

        N = 10*2^(i-1); % (0:3)

        %N = 100;
        dx = 1/N;
        x = (1:N-1)*dx;  % grid for a single side

        A =  N^2 * diag( ( y([x+dx/2]) + y([x-dx/2]) ));
        A = A - N^2 * diag(y([x(2:end)-dx/2]),-1);
        A = A - N^2 * diag(y([x(1:end-1)+dx/2]),1);


        f= x' .* x'; % ones(N-1,1);

        u = A\f;
        %errs = u' - s(x);
        errs = u' - h(x);
        l2_err(i)  = sqrt(sum(errs.^2));

        figure(1)
        subplot(2,3,j+1)
        plot([0,x,1], [0,u',0])
        title(['U_{h,eps}, eps = ' ,num2str(eps)]);
        hold on
        subplot(2,3,4+j)
        plot([0,x,1],[0,h(x),0])
        title('U_{hom}');

%         hold on
%         subplot(2,1,2)
%         plot([0,x,1],[0,errs,0] )
%         %hold off
%         title('l2 norm of u_{h,eps} - u_{hom}')
          J = N/10;
        usv(:,i) = u(J:J:N-1);

    end
    %suptitle(' Graphs of U_{h,eps} and U_{hom}')
    figure(2)
    subplot(3,1,j+1)
    plot(l2_err)
    title(['l2 error of U_{h,eps} - U_{hom} for eps = ',num2str(eps)])
    ylabel('l2 errors')
end
    xlabel('Mesh size h = 1/(10*2^{i-1}) with i = 1:8')

suptitle('Graphs of l2 Errors vs mesh size (in terms of i) for 3 Different Epsilon')
hold off
%legend([{'1'},{'2'},{'3'},{'4'},{'5'},{'6'}])
