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

h =  @(i) 1e-2./10.^((i-1)/20);

N = 61;

RK4Conv.TwoLoop.endpoints = zeros(4,N);
RK4Conv.ThreeLoop.endpoints = zeros(4,N);
RK4Conv.FourLoop.endpoints = zeros(4,N);

for i = 1:N

% Two loops
TwoLoop.RK4 = RK4(@ThreeBodyProblem,tInit(1),tEnd(1),ICs(1),h(i));

% Three loop
ThreeLoop.RK4 = RK4(@ThreeBodyProblem,tInit(2),tEnd(2),ICs(2),h(i));

% Three loop
FourLoop.RK4 = RK4(@ThreeBodyProblem,tInit(3),tEnd(3),ICs(3),h(i));

RK4Conv.TwoLoop.endpoints(:,i) = TwoLoop.RK4.y(:,end);
RK4Conv.ThreeLoop.endpoints(:,i) = ThreeLoop.RK4.y(:,end);
RK4Conv.FourLoop.endpoints(:,i) = FourLoop.RK4.y(:,end);

end


% Plotting
tiledlayout(2,2)
nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(1,21:N))
hold on 
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(1,21:N))
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(1,21:N))
set(gca, 'XDir', 'reverse')

nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(2,21:N))
hold on 
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(2,21:N))
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(2,21:N))
set(gca, 'XDir', 'reverse')

nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(3,21:N))
hold on 
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(3,21:N))
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(3,21:N))
set(gca, 'XDir', 'reverse')

nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(4,21:N))
hold on  
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(4,21:N))
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(4,21:N))
set(gca, 'XDir', 'reverse')
