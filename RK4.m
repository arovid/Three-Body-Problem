function out = RK4(fun,tInit,tEnd,yInit,h)
% RK4 solver
%   Inputs: 
%       - fun:          System of ODEs to solve (function defined elsewhere)
%       - tInit:        Starting time
%       - tEnd:         Ending time
%       - yInit:        Initial conditions of the system
%       - h:            Step size
%   Outputs:
%       - out.t:        Time steps
%       - out.y:        Solutions at out.t times
%       - out.fevals:   Number of function evaluations

N = floor((tEnd-tInit)/h);      % Step number
t = linspace(tInit,N*h,N+1);    % Grid
y = zeros(4,N+1);               % Vector of the numerical approximation
y(:,1) = yInit;                 % Initial conditions
nfeval = 0;                     % Function evaluations counter

for j=1:N

    y_current = y(:,j);
    t_current = t(j);

     k1 = fun(t_current,y_current);
    k2 = fun(t_current,y_current + 0.5*h*k1);
    k3 = fun(t_current,y_current + 0.5*h*k2);
    k4 = fun(t_current,y_current + h*k3);
    nfeval = nfeval + 4;

    y(:,j+1) = (y_current + h/6 * (k1 + 2*k2 + 2*k3 + k4));
end
    out.y = y;
    out.t = t;
    out.nsteps = nfeval;
end