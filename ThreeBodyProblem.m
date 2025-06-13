function dydt = ThreeBodyProblem(t,y)
% Three-body problem
%   System of 1st order, non-linear differential equations representing the
%   three-body problem

mu = 0.012277471; 
nu = 1-mu;

D1 = ((y(1) + mu)^2 + y(2)^2)^(3/2);
D2 = ((y(1) - nu)^2 + y(2)^2)^(3/2);

dydt = zeros(4,1);
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = y(1) + 2*y(4) - nu*(y(1) + mu)/D1 - mu*(y(1) - nu)/D2;
dydt(4) = y(2) - 2*y(3) - nu*y(2)/D1 - mu*y(2)/D2;

end