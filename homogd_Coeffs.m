    %
    %       - \grad_y \cdot ( a(x,y) \grad_y\chi) = \grad_y(a(x,y))
    %               for all x of interest solve for \chi(x,y) for y \in \IT^d
    %
    %           Compute Abar = \int_{\IT^d} (a(x,y) + grad_y\chi(a(x,y))dy)
    %           Used in the Finite Diff. Solvers.
 function homog

    function z = ay(x,y)
            z = 1.1 + sin(2 * pi * x(:,1)) .* sin(2*pi * y(:,1));
    end

    function z = d2ay(x,y)
            z = 2* pi * sin(2 * pi * x(:,1)) .* cos(2* pi * y(:,1));
    end

    x = [0.1, 0.1];
    Ahom(x)

    function z = Ahom(x)
            eps = 10^(-2);
            i=2;
            N = 4 * 2^i; %(0:3)
            dx = 1/N;
            s = (1:N)*dx;  % grid for a single side
            [y1,y2] = meshgrid(s,s);

            y1=y1(:);
            y2=y2(:);


            M = N^2;

            %printf('N %d M %d\n', N, M)
            inds = reshape(1:M, [N, N]);
            tic();

            A =  N^2 * sparse((1:M)', (1:M)', ( ay(x, [y1+dx/2,y2]) + ...
            ay(x,[y1-dx/2,y2]) + ay(x,[y1, y2-dx/2]) + ay(x,[y1,y2+dx/2])));

            indm = [inds(:,end),inds(:,1:end-1)];indm = indm(:);
            indp = [inds(:,2:end),inds(:,1)];indp = indp(:);
            ind = inds(:);
            A = A - N^2 * sparse(ind,indm, ay(x, [y1(ind)-dx/2,y2(ind)]), M,M);
            A = A - N^2 * sparse(ind,indp, ay(x, [y1(ind)+dx/2,y2(ind)]), M,M);

            indp = [inds(end,:);inds(1:end-1,:)];indp = indp(:);
            indm = [inds(2:end,:);inds(1,:)];  indm = indm(:);
            A = A - N^2 * sparse(ind, indm, ay(x,[y1(ind), y2(ind)-dx/2]), M, M);
            A = A - N^2 * sparse(ind, indp, ay(x,[y1(ind), y2(ind)+dx/2]), M, M);

            chi  = A\d2ay(x,[y1,y2]); %   \ \grad a(x,y) solve for chi, solution to elliptic PDE.
            Abar = diag([(trapz((1./ay(x,[y1,y2]/eps)))/size(chi,1) +...
                trapz((1./ay(x,[y1,y2]/eps)).*[diff(chi);0])/size(chi,1))^(-1) ...
                trapz(ay(x,[y1,y2]/eps))/size(chi,1) +...
                trapz(ay(x,[y1,y2]/eps).*[diff(chi);0])/size(chi,1)]) ;

            ff=@(x1) double((10*1i)/((10*sin(2*pi*x1) - 11)^(1/2)*(10*sin(2*pi*x1) + 11)^(1/2)));
            av = (ff(x(1))); %  x = [0.1,0.1]

            AbarAn = 1/av;
%             Abar = diag([trapz(ay(x,[y1,y2]/eps))/size(chi,1) +...
%                 trapz(ay(x,[y1,y2]/eps).*[diff(chi);0])/size(chi,1) ...
%                 (trapz((1./ay(x,[y1,y2]/eps)))/size(chi,1) +...
%                 trapz((1./ay(x,[y1,y2]/eps)).*[diff(chi);0])/size(chi,1))^(-1)]);
            z = Abar(1:2,1:2);  %diag([Abar(1,1) Abar(3,3)]) %
    end
 end
