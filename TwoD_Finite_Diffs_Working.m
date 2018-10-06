%
%     2D versoin of OneD_Finite_Diffs.
%
%     2D MultiGrid  Num. Solver for AU = F (assume discretized, so linear
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

function rrr

% function y = a(x)
%     y=  1+ (1.1 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)/eps)).*x(:,1);
% end

function y = a(x)
    y= 1+ 0.9192.*x(:, 1); %0.9192 calculated homoegnized coefficien from homogd_coeff.
    % compare with integrated coefficient using integ.m
end

for j = 2
    eps = 10^-j;

    irng = 6;
    for i = irng

        N = 4 * 2^(i+2); %(0:3)
        dx = 1/N;
        s = (1:N-1)*dx;  % grid for a single side
        [X,Y] = meshgrid(s,s);

        x=X(:);
        y=Y(:);
        M = (N-1)^2;

        fprintf('N %d M %d\n', N, M)
        inds = reshape(1:M, N-1, N-1);
        tic();

                A =  N^2 * sparse((1:M)', (1:M)', ( a([x+dx/2,y]) + a([x-dx/2,y]) + a([x, y-dx/2]) + a([x,y+dx/2])));
        indm = inds(:,1:end-1);indm = indm(:);
        indp = inds(:,2:end);indp = indp(:);
        A = A - N^2 * sparse(indm, indp, a([x(indm)+dx/2, y(indm)]), M, M);
        A = A - N^2 * sparse(indp, indm, a([x(indp)-dx/2, y(indp)]), M, M);

        indp = inds(1:end-1,:);indp = indp(:);
        indm = inds(2:end,:);indm = indm(:);
        A = A - N^2 * sparse(indm, indp, a([x(indm), y(indm)+dx/2]), M, M);
        A = A - N^2 * sparse(indp, indm, a([x(indp), y(indp)-dx/2]), M, M);



        fprintf('Build %f\n',toc())



        %f=  2* x .* (1-x) + 2 * y .* (1 - y);  %ones(N-1,1);
        %u_ex = x.*(1-x).*y.*(1-y);

        %f=  1000* (2* (exp(x)-1) .* (1-x) + (1+x) .* exp(x) .* y .* (1 - y));  %ones(N-1,1);
        u_ex = 1000*(exp(x)-1).*(1-x).*y.*(1-y); 

        f =x.^2;

        tic();


         u =A\f;% umfpack(A,"\",f)
        u0=ones(size(f,1),1);% I'M THE INITIAL CONDITION
       %u= jacobi(A,f, u0, 10000);%  CHANGE ME FOR JACOBI!!!

        fprintf('umfpack   Solve %f\n', toc())

        fprintf('L2 error direct %f   Int %f \n', norm((u-u_ex)*dx), norm((u_ex)*dx))
   shg;
 % subplot(2,ceil(size(irng,2)),i)
U=reshape(u, N-1, N-1);
 s = surf(X,Y,U);
% alpha 0.7 ;
% colormap winter ;
 s.EdgeColor = 'None';

%Titles for comparison plots
%title(['dx = ',num2str(dx)]);
%    if i == 4
%        zlabel({'Equation with Numerically',' Computed Averaged coefficient'});
%
%    end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Repeat of the Above Code for Computing the Errors %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function y = b(x)
%     y=  1+ (1.1 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)/eps)).*x(:,1);
% end
%
% for j = 2
%     eps = 10^-j;
%
%     irng = 6;
%     for i = irng
%
%         N = 4 * 2^(i+2); %(0:3)
%         dx = 1/N;
%         s = (1:N-1)*dx;  % grid for a single side
%         [X,Y] = meshgrid(s,s);
%
%         x=X(:);
%         y=Y(:);
%         M = (N-1)^2;
%
%         fprintf('N %d M %d\n', N, M)
%         inds = reshape(1:M, N-1, N-1);
%         tic();
%
%                 A =  N^2 * sparse((1:M)', (1:M)', ( b([x+dx/2,y]) + b([x-dx/2,y]) + b([x, y-dx/2]) + b([x,y+dx/2])));
%         indm = inds(:,1:end-1);indm = indm(:);
%         indp = inds(:,2:end);indp = indp(:);
%         A = A - N^2 * sparse(indm, indp, b([x(indm)+dx/2, y(indm)]), M, M);
%         A = A - N^2 * sparse(indp, indm, b([x(indp)-dx/2, y(indp)]), M, M);
%
%         indp = inds(1:end-1,:);indp = indp(:);
%         indm = inds(2:end,:);indm = indm(:);
%         A = A - N^2 * sparse(indm, indp, b([x(indm), y(indm)+dx/2]), M, M);
%         A = A - N^2 * sparse(indp, indm, b([x(indp), y(indp)-dx/2]), M, M);
%
%
%
%         fprintf('Build %f\n',toc())
%
%
%
%         %f=  2* x .* (1-x) + 2 * y .* (1 - y);  %ones(N-1,1);
%         %u_ex = x.*(1-x).*y.*(1-y);
%
%         %f=  1000* (2* (exp(x)-1) .* (1-x) + (1+x) .* exp(x) .* y .* (1 - y));  %ones(N-1,1);
%         u_ex = 1000*(exp(x)-1).*(1-x).*y.*(1-y);
%
%         f =x.^2;
%
%         tic();
%
%
%          u =A\f;% umfpack(A,"\",f)
%         u0=ones(size(f,1),1);% I'M THE INITIAL CONDITION
%        %u= jacobi(A,f, u0, 10000);%  CHANGE ME FOR JACOBI!!!
%
%         fprintf('umfpack   Solve %f\n', toc())
%
%         fprintf('L2 error direct %f   Int %f \n', norm((u-u_ex)*dx), norm((u_ex)*dx))
%    shg;
% %   subplot(2,ceil(size(irng,2)),i)
% U2=reshape(u, N-1, N-1);
% % s = surf(X,Y,U);
% % alpha 0.7 ;
% % colormap winter ;
% % s.EdgeColor = 'None';
%
% % function y = a(x)
% %     y=  1+ (1.1 + sin(2*pi*x(:,1)).*sin(2*pi*x(:,1)/eps)).*x(:,1);
% % end
%
% %Titles for comparison plots
%
% %title(['dx = ',num2str(dx)]);
% %    if i == 4
% %        zlabel({'Equation with Analytically Found',' Averaged coefficient'});
% %
% %    end
%     end
%
% end
%
% s = surf(X,Y, abs(U-U2));
% s.EdgeColor = 'None';
%





 end
% hs = [32,64,128,256,512,1024]
% J = [0.0041,0.0145851,5.72633,25.849968, 37.77,41.704210]
% MG  = [0.0041,0.001033,0.001117,0.001315,0.001269, 0.00121]
% D = [0.0041,0.001025,0.000256,0.000064,0.000016,0.000004]
