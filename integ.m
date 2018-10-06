

syms eps pi  x1 x2
%  a=@(x1) 1.1 + sin(2*pi*x1)*sin(2*pi*x2/eps) 
 eps=vpa(sym('eps','real'))
 x1=vpa(sym('x1','real'))
  a=@(x2) 1/(1.1 + sin(2 * pi * .1) * sin(2* pi * x2));
 
% int(int(diff(a(x1,x2),x2),x2),x2)

(int(1/(1.1 + sin(2 * pi * x1) * sin(2* pi * x2)),x2,0,1)) % int ( 1/a(x,y))


ff=@(x1) double((10*1i)/((10*sin(2*pi*x1) - 11)^(1/2)*(10*sin(2*pi*x1) + 11)^(1/2)));
av = (ff(.1)); %  x = [0.1,0.1]

harAvg = 1/av;

Abar = diag([harAvg 1.1]) % diag([1.1 harAvg])


