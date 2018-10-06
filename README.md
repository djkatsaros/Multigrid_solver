#  Multigrid_solver

Preamble to Code:

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

Code is split into several functions. 
Project is motivated by theory that can be found in  "Multiscale Methods", G.A. Pavliotis and A.M. Stuart.

Compare the solutions to "homogenized equation" in the spirit of the book to the equation with coefficient
given by direct integration. 

Compare performance when using multigrid numerical solver (https://en.wikipedia.org/wiki/Multigrid_method) 
vs direct or Jacobi  solvers. In particular, it was obeserved that the multigrid solver was faster than both 
a direct solver via LU decomposition (O(n^3) in this problem) and the Jacboi iteration - see post project report. 

Much of the study was repeated twice: once each for the 1d and 2d versions of the problem respectively. 


