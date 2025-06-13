clear
clc
format long

% Set the initial conditions and parameters for the simulations
tInit = [0; 0; 0];
tEnd = [6.2; 11.2; 17.1];
u1Init = [1.2; 0.994; 0.994];
u2dInit = [-1.049357510; -2.0317326295573368357302057924; -2.00158510637908252240537862224];
ICs = @(i)[u1Init(i); 0; 0; u2dInit(i)];
AbsTol = 1e-4;
RelTol = 1e-4;

h_RK4 = 1e-5;

opts = odeset('AbsTol',1e-4,'RelTol',1e-12,'Stats','on');

% Two loops
t1 = tic;
TwoLoop.RK4 = RK4(@ThreeBodyProblem,tInit(1),tEnd(1),ICs(1),h_RK4);
TwoLoop.RK4.CompTime = toc(t1);

t2 = tic;
TwoLoop.ode45 = ode45(@ThreeBodyProblem,[tInit(1) tEnd(1)], ICs(1),opts);
TwoLoop.ode45.CompTime = toc(t2);

t3 = tic;
TwoLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(1),tEnd(1),ICs(1),AbsTol,RelTol);
TwoLoop.RKF45.CompTime = toc(t3);

% Three loop
t4 = tic;
ThreeLoop.RK4 = RK4(@ThreeBodyProblem,tInit(2),tEnd(2),ICs(2),h_RK4);
ThreeLoop.RK4.CompTime = toc(t1);

t5 = tic;
ThreeLoop.ode45 = ode45(@ThreeBodyProblem,[tInit(2) tEnd(2)], ICs(2),opts);
ThreeLoop.ode45.CompTime = toc(t5);

t6 = tic;
ThreeLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(2),tEnd(2),ICs(2),AbsTol,RelTol);
ThreeLoop.RKF45.CompTime = toc(t6);

% Four loop
t7 = tic;
FourLoop.RK4 = RK4(@ThreeBodyProblem,tInit(3),tEnd(3),ICs(3),h_RK4); 
FourLoop.RK4.CompTime = toc(t7);

t8 = tic;
FourLoop.ode45 = ode45(@ThreeBodyProblem,[tInit(3) tEnd(3)], ICs(3),opts);
FourLoop.ode45.CompTime = toc(t8);

t9 = tic;
FourLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(3),tEnd(3),ICs(3),AbsTol,RelTol);
FourLoop.RKF45.CompTime = toc(t9);