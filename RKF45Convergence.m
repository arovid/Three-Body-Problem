clear
clc
format long
clf

% Set the initial conditions and parameters for the simulations
tInit = [0; 0; 0];
tEnd = [6.2; 11.2; 17.1];
u1Init = [1.2; 0.994; 0.994];
u2dInit = [-1.049357510; -2.0317326295573368357302057924; -2.00158510637908252240537862224];
ICs = @(i)[u1Init(i); 0; 0; u2dInit(i)];

Tol =  @(i) 1e-2./10.^((i-1)/2);
Ni = 9;
Nj = 17;

RKF45Conv.TwoLoop.endpoints = zeros(4,Ni,Nj);
RKF45Conv.ThreeLoop.endpoints = zeros(4,Ni,Nj);
RKF45Conv.FourLoop.endpoints = zeros(4,Ni,Nj);

for i = 1:Ni
    for j = 1:Nj
        % Two loops
        TwoLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(1),tEnd(1),ICs(1),Tol(i),Tol(j));

        % Three loop
        ThreeLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(2),tEnd(2),ICs(2),Tol(i),Tol(j));

        % Three loop
        FourLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(3),tEnd(3),ICs(3),Tol(i),Tol(j));

        RKF45Conv.TwoLoop.endpoints(:,i,j) = TwoLoop.RKF45.y(:,end);
        RKF45Conv.ThreeLoop.endpoints(:,i,j) = ThreeLoop.RKF45.y(:,end);
        RKF45Conv.FourLoop.endpoints(:,i,j) = FourLoop.RKF45.y(:,end);
    end
end

surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(1,:,:))')