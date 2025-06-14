function out = RKF45(fun,tInit,tEnd,yInit,RelTol,AbsTol) 
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

T = tInit;  % Grid
Y = yInit;  % Vector of the numerical approximation (4th order)
Z = yInit;  % Vector of the numerical approximation (5th order)

t = T(1);
y = Y(:,1);


h = (tEnd-tInit)/1000; % Guess for initial step size
 
nsteps = 0;
nfeval = 0;
nfailed = 0;

% Computing the numerical approximation using RK4:
while t < tEnd
    k1 = fun(t, y);
    k2 = fun(t + 1/4*h, y + h*(1/4*k1));
    k3 = fun(t + 3/8*h, y + h*(3/32*k1 + 9/32*k2));
    k4 = fun(t + 12/13*h, y + h*(1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3));
    k5 = fun(t + 1*h, y + h*(439/216*k1 - 8*k2 +  3680/513*k3 - 845/4104*k4));
    k6 = fun(t + 1/2*h, y + h*(-8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5));

    yTemp = y + h*(25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5);
    zTemp = y + h*(16/135*k1 + 6656/12825*k3 + 28651/56430*k4 - 9/50*k5 + 2/55*k6);

    e = abs(yTemp-zTemp);                   % Error in the step
    etol = RelTol*abs(yTemp) + AbsTol;      % Error tolerance defined by user

    if all(e<etol)
        y = yTemp;
        z = zTemp;
        t = t + h;
        Y(:,end+1) = y;
        Z(:,end+1) = z;
        T(end+1) = t;
        nsteps = nsteps + 1;
    else
        nfailed = nfailed + 1;
    end

    frac = 0.9;
    h = min(frac*h*(etol./e).^(1/5));    
    
    if t+h > tEnd
        h = tEnd - t;
    end
    nfeval = nfeval + 6;
    
end

% out.t = cell2mat(T);
% out.y = cell2mat(Y);
% out.z = cell2mat(Z);
out.t = T;
out.y = Y;
out.z = Z;
out.stats.nfeval = nfeval;
out.stats.nsteps = nsteps;
out.stats.nfailed = nfailed;
end
