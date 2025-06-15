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
fig1 = figure(1);
fig1.Position = [100 100 1600 1200];
tiledlayout(2,2,TileSpacing="compact")
nexttile
plot(log10(h(1:N)),RK4Conv.TwoLoop.endpoints(1,:),"LineWidth",2)
hold on 
plot(log10(h(1:N)),RK4Conv.ThreeLoop.endpoints(1,:),"LineWidth",2)
plot(log10(h(1:N)),RK4Conv.FourLoop.endpoints(1,:),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_1 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'FourLoop'
set(gca, 'TickLabelInterpreter','latex')

nexttile
plot(log10(h(1:N)),RK4Conv.TwoLoop.endpoints(2,:),"LineWidth",2)
hold on 
plot(log10(h(1:N)),RK4Conv.ThreeLoop.endpoints(2,:),"LineWidth",2)
plot(log10(h(1:N)),RK4Conv.FourLoop.endpoints(2,:),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_2 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'FourLoop'
set(gca, 'TickLabelInterpreter','latex')

nexttile
plot(log10(h(1:N)),RK4Conv.TwoLoop.endpoints(3,:),"LineWidth",2)
hold on 
plot(log10(h(1:N)),RK4Conv.ThreeLoop.endpoints(3,:),"LineWidth",2)
plot(log10(h(1:N)),RK4Conv.FourLoop.endpoints(3,:),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_3 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'FourLoop'
set(gca, 'TickLabelInterpreter','latex')

nexttile
plot(log10(h(1:N)),RK4Conv.TwoLoop.endpoints(4,:),"LineWidth",2)
hold on  
plot(log10(h(1:N)),RK4Conv.ThreeLoop.endpoints(4,:),"LineWidth",2)
plot(log10(h(1:N)),RK4Conv.FourLoop.endpoints(4,:),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_4 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'Four Loop'
set(gca, 'TickLabelInterpreter','latex')


fig2 = figure(2);
fig2.Position = [100 100 1600 1200];
tiledlayout(2,2,TileSpacing="tight")
nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(1,21:N),"LineWidth",2)
hold on 
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(1,21:N),"LineWidth",2)
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(1,21:N),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_1 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'Four Loop'
set(gca, 'TickLabelInterpreter','latex')

nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(2,21:N),"LineWidth",2)
hold on 
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(2,21:N),"LineWidth",2)
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(2,21:N),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_2 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'Four Loop'
set(gca, 'TickLabelInterpreter','latex')

nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(3,21:N),"LineWidth",2)
hold on 
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(3,21:N),"LineWidth",2)
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(3,21:N),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_3 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'FourLoop'
set(gca, 'TickLabelInterpreter','latex')

nexttile
plot(log10(h(21:N)),RK4Conv.TwoLoop.endpoints(4,21:N),"LineWidth",2)
hold on  
plot(log10(h(21:N)),RK4Conv.ThreeLoop.endpoints(4,21:N),"LineWidth",2)
plot(log10(h(21:N)),RK4Conv.FourLoop.endpoints(4,21:N),"LineWidth",2)
set(gca, 'XDir', 'reverse')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
ylabel('$$u_4 [-]$$')
xlabel('$$ \lg(h) [-]$$')
legend 'Two loop' 'Three loop' 'Four Loop'
set(gca, 'TickLabelInterpreter','latex')


exportgraphics(fig1, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RK4Convergence_full.pdf', 'ContentType', 'vector');
exportgraphics(fig2, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RK4Convergence_zoom.pdf', 'ContentType', 'vector');