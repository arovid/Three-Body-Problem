clear
clc
format long

% Set the initial conditions and parameters for the simulations
tInit = [0; 0; 0];
tEnd = [6.2; 11.2; 17.1];
u1Init = [1.2; 0.994; 0.994];
u2dInit = [-1.049357510; -2.0317326295573368357302057924; -2.00158510637908252240537862224];
ICs = @(i)[u1Init(i); 0; 0; u2dInit(i)];
h = 10^(-3.9);

opts = odeset('AbsTol',1e-8,'RelTol',1e-8,'Stats','on');

% RK4 - Two loops
TwoLoop.RK4 = RK4(@ThreeBodyProblem,tInit(1),tEnd(1),ICs(1),h);

% RK4 - Three loop
ThreeLoop.RK4 = RK4(@ThreeBodyProblem,tInit(2),tEnd(2),ICs(2),h);

% RK4 - Four loop10e
FourLoop.RK4 = RK4(@ThreeBodyProblem,tInit(3),tEnd(3),ICs(3),h);

% ode45 - Two loops
TwoLoop.ode45 = ode45(@ThreeBodyProblem,[tInit(1) tEnd(1)], ICs(1),opts);

% ode45 - Three loop
ThreeLoop.ode45 = ode45(@ThreeBodyProblem,[tInit(2) tEnd(2)], ICs(2),opts);

% ode45 - Four loop
FourLoop.ode45 = ode45(@ThreeBodyProblem,[tInit(3) tEnd(3)], ICs(3),opts);

% RKF45 - Two loops
TwoLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(1),tEnd(1),ICs(1),1e-5,10^(-4.8));

% RKF45 - Three loop
ThreeLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(2),tEnd(2),ICs(2),1e-5,10^(-4.8));

% RKF45 - Four loop
FourLoop.RKF45 = RKF45(@ThreeBodyProblem,tInit(3),tEnd(3),ICs(3),1e-5,10^(-4.8));




timeVectors = {TwoLoop.RK4.t, ThreeLoop.RK4.t, FourLoop.RK4.t;
               TwoLoop.ode45.x, ThreeLoop.ode45.x, FourLoop.ode45.x;
               TwoLoop.RKF45.t, ThreeLoop.RKF45.t, FourLoop.RKF45.t};
StateVectors = {TwoLoop.RK4.y(:,:), ThreeLoop.RK4.y(:,:), FourLoop.RK4.y(:,:);
                TwoLoop.ode45.y(:,:), ThreeLoop.ode45.y(:,:), FourLoop.ode45.y(:,:);
                TwoLoop.RKF45.y(:,:), ThreeLoop.RKF45.y(:,:), FourLoop.RKF45.y(:,:)};
colors = ['b', 'g', 'r']; 

for i = 1:3
    fig(i) = figure(i);
    fig(i).Position = [40 40 1200 800];
    tiledlayout(4,1,TileSpacing="compact")
    for j = 1:4
        nexttile
        for k = 1:3
            t = cell2mat(timeVectors(k,i));
            x = cell2mat(StateVectors(k,i));
            plot(t,x(j,:),colors(k),'LineWidth',2);
            hold on
            grid on
            xlabel('t [s]')
            if j < 3
                ylabel(strcat('$$u_',string(j),' [-]$$'))
            else
                ylabel(strcat('$$\dot{u}_',string(j-2),' [-]$$'))
            end
        end
        plotTfinal = [8 14 20];
        xlim([tInit(i) plotTfinal(i)]);
        legend 'RK4' 'ode45' 'RKF45'
    end

    set(findall(fig(i),'-property','FontSize'),'FontSize', 16)
    set(findall(fig(i),'-property','Box'),'Box', 'on')
    set(findall(fig(i),'-property','Interpreter'),'Interpreter', 'latex')
    set(findall(fig(i),'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex')
end

exportgraphics(fig(1), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\States_TwoLoop.png', 'ContentType', 'vector');
exportgraphics(fig(2), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\States_ThreeLoop.png', 'ContentType', 'vector');
exportgraphics(fig(3), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\States_FourLoop.png', 'ContentType', 'vector');

close all

mu = 0.012277471;
for i = 1:3
    fig(i) = figure(i);
    fig(i).Position = [40 40 600 600];
    for j =1:3
        x = cell2mat(StateVectors(j,i));
        plot(x(1,:),x(2,:),colors(j),'LineWidth',2);
        hold on
    end
    scatter(-mu,0,40,'k','filled')
    scatter(1-mu,0,40,'k','filled')
    theta = 0:0.01:2*pi;
    plot((1-mu)*cos(theta),(1-mu)*sin(theta),'k--','LineWidth',1.5)
    plot(mu*cos(theta),mu*sin(theta),'k--','LineWidth',1.5)
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    grid on
    axis equal
    xlabel('$$u_1 [-]$$')
    ylabel('$$u_2 [-]$$')
    legend 'RK4' 'ode45' 'RKF45'
    set(findall(fig(i),'-property','FontSize'),'FontSize', 16)
    set(findall(fig(i),'-property','Box'),'Box', 'on')
    set(findall(fig(i),'-property','Interpreter'),'Interpreter', 'latex')
    set(findall(fig(i),'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex')
end

exportgraphics(fig(1), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\Orbit_TwoLoop.png', 'ContentType', 'vector');
exportgraphics(fig(2), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\Orbit_ThreeLoop.png', 'ContentType', 'vector');
exportgraphics(fig(3), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\Orbit_FourLoop.png', 'ContentType', 'vector');
close all

% Comparison 
Reference.Loop(1).x  = cell2mat(timeVectors(2,1));
Reference.Loop(1).u  = cell2mat(StateVectors(2,1));
Reference.Loop(2).x  = cell2mat(timeVectors(2,2));
Reference.Loop(2).u  = cell2mat(StateVectors(2,2));
Reference.Loop(3).x  = cell2mat(timeVectors(2,3));
Reference.Loop(3).u  = cell2mat(StateVectors(2,3));

% RK4
x = cell2mat(timeVectors(1,1));
v = cell2mat(StateVectors(1,1));
Comp.Loop(1).RK4.u = zeros(4,length(Reference.Loop(1).x));
Comp.Loop(1).RK4.u(1,:) = interp1(x,v(1,:),Reference.Loop(1).x,"spline");
Comp.Loop(1).RK4.u(2,:) = interp1(x,v(2,:),Reference.Loop(1).x,"spline");
Comp.Loop(1).RK4.u(3,:) = interp1(x,v(3,:),Reference.Loop(1).x,"spline");
Comp.Loop(1).RK4.u(4,:) = interp1(x,v(4,:),Reference.Loop(1).x,"spline");
Loop(1).RK4.u = zeros(4,length(TwoLoop.RK4.t));
Loop(1).RK4.t = x;
Loop(1).RK4.u(1,:) = v(1,:);
Loop(1).RK4.u(2,:) = v(2,:);
Loop(1).RK4.u(3,:) = v(3,:);
Loop(1).RK4.u(4,:) = v(4,:);

x = cell2mat(timeVectors(1,2));
v = cell2mat(StateVectors(1,2));
Comp.Loop(2).RK4.u = zeros(4,length(Reference.Loop(2).x));
Comp.Loop(2).RK4.u(1,:) = interp1(x,v(1,:),Reference.Loop(2).x,"spline");
Comp.Loop(2).RK4.u(2,:) = interp1(x,v(2,:),Reference.Loop(2).x,"spline");
Comp.Loop(2).RK4.u(3,:) = interp1(x,v(3,:),Reference.Loop(2).x,"spline");
Comp.Loop(2).RK4.u(4,:) = interp1(x,v(4,:),Reference.Loop(2).x,"spline");
Loop(2).RK4.u = zeros(4,length(ThreeLoop.RK4.t));
Loop(2).RK4.t = x;
Loop(2).RK4.u(1,:) = v(1,:);
Loop(2).RK4.u(2,:) = v(2,:);
Loop(2).RK4.u(3,:) = v(3,:);
Loop(2).RK4.u(4,:) = v(4,:);

x = cell2mat(timeVectors(1,3));
v = cell2mat(StateVectors(1,3));
Comp.Loop(3).RK4.u = zeros(4,length(Reference.Loop(3).x));
Comp.Loop(3).RK4.u(1,:) = interp1(x,v(1,:),Reference.Loop(3).x,"spline");
Comp.Loop(3).RK4.u(2,:) = interp1(x,v(2,:),Reference.Loop(3).x,"spline");
Comp.Loop(3).RK4.u(3,:) = interp1(x,v(3,:),Reference.Loop(3).x,"spline");
Comp.Loop(3).RK4.u(4,:) = interp1(x,v(4,:),Reference.Loop(3).x,"spline");
Loop(3).RK4.u = zeros(4,length(FourLoop.RK4.t));
Loop(3).RK4.t = x;
Loop(3).RK4.u(1,:) = v(1,:);
Loop(3).RK4.u(2,:) = v(2,:);
Loop(3).RK4.u(3,:) = v(3,:);
Loop(3).RK4.u(4,:) = v(4,:);

% RKF45
x = cell2mat(timeVectors(3,1));
v = cell2mat(StateVectors(3,1));
Comp.Loop(1).RKF45.u = zeros(4,length(Reference.Loop(1).x));
Comp.Loop(1).RKF45.u(1,:) = interp1(x,v(1,:),Reference.Loop(1).x,"spline");
Comp.Loop(1).RKF45.u(2,:) = interp1(x,v(2,:),Reference.Loop(1).x,"spline");
Comp.Loop(1).RKF45.u(3,:) = interp1(x,v(3,:),Reference.Loop(1).x,"spline");
Comp.Loop(1).RKF45.u(4,:) = interp1(x,v(4,:),Reference.Loop(1).x,"spline");
Loop(1).RKF45.u = zeros(4,length(TwoLoop.RKF45.t));
Loop(1).RKF45.t = x;
Loop(1).RKF45.u(1,:) = v(1,:);
Loop(1).RKF45.u(2,:) = v(2,:);
Loop(1).RKF45.u(3,:) = v(3,:);
Loop(1).RKF45.u(4,:) = v(4,:);

x = cell2mat(timeVectors(3,2));
v = cell2mat(StateVectors(3,2));
Comp.Loop(2).RKF45.u = zeros(4,length(Reference.Loop(2).x));
Comp.Loop(2).RKF45.u(1,:) = interp1(x,v(1,:),Reference.Loop(2).x,"spline");
Comp.Loop(2).RKF45.u(2,:) = interp1(x,v(2,:),Reference.Loop(2).x,"spline");
Comp.Loop(2).RKF45.u(3,:) = interp1(x,v(3,:),Reference.Loop(2).x,"spline");
Comp.Loop(2).RKF45.u(4,:) = interp1(x,v(4,:),Reference.Loop(2).x,"spline");
Loop(2).RKF45.u = zeros(4,length(ThreeLoop.RKF45.t));
Loop(2).RKF45.t = x;
Loop(2).RKF45.u(1,:) = v(1,:);
Loop(2).RKF45.u(2,:) = v(2,:);
Loop(2).RKF45.u(3,:) = v(3,:);
Loop(2).RKF45.u(4,:) = v(4,:);


x = cell2mat(timeVectors(3,3));
v = cell2mat(StateVectors(3,3));
Comp.Loop(3).RKF45.u = zeros(4,length(Reference.Loop(3).x));
Comp.Loop(3).RKF45.u(1,:) = interp1(x,v(1,:),Reference.Loop(3).x,"spline");
Comp.Loop(3).RKF45.u(2,:) = interp1(x,v(2,:),Reference.Loop(3).x,"spline");
Comp.Loop(3).RKF45.u(3,:) = interp1(x,v(3,:),Reference.Loop(3).x,"spline");
Comp.Loop(3).RKF45.u(4,:) = interp1(x,v(4,:),Reference.Loop(3).x,"spline");
Loop(3).RKF45.u = zeros(4,length(FourLoop.RKF45.t));
Loop(3).RKF45.t = x;
Loop(3).RKF45.u(1,:) = v(1,:);
Loop(3).RKF45.u(2,:) = v(2,:);
Loop(3).RKF45.u(3,:) = v(3,:);
Loop(3).RKF45.u(4,:) = v(4,:);


% plot abs(ode45-otherMethod) as a function of time
for i = 1:3
    fig(i+3) = figure(i+3);
    fig(i+3).Position = [40 40 1200 800];
    tiledlayout(4,1,TileSpacing="compact")
    for j = 1:4
        nexttile
        plot(Reference.Loop(i).x, abs(Comp.Loop(i).RK4.u(j,:) - Reference.Loop(i).u(j,:)),'r','LineWidth',2)
        hold on
        plot(Reference.Loop(i).x, abs(Comp.Loop(i).RKF45.u(j,:) - Reference.Loop(i).u(j,:)),'k','LineWidth',2)
    end
    set(findall(fig(i+3),'-property','FontSize'),'FontSize', 16)
    set(findall(fig(i+3),'-property','Box'),'Box', 'on')
    set(findall(fig(i+3),'-property','Interpreter'),'Interpreter', 'latex')
    set(findall(fig(i+3),'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex')
end

% Jacobi constant - conserved
mu = 0.012277471;
nu = 1 - mu;
Theta = @(x,y) 1/2*(x.^2 + y.^2) + nu./sqrt((x+mu).^2 + y.^2) + mu./sqrt((x-nu).^2 + y.^2);
C = @(x,y,xd,yd) 2*Theta(x,y) - (xd.^2 + yd.^2);


fig(7) = figure(7);
fig(7).Position = [40 40 1200 800];
tiledlayout(3,1,TileSpacing="compact")
for i = 1:3
nexttile
x = Loop(i).RK4.u(1,:);
y = Loop(i).RK4.u(2,:);
xd = Loop(i).RK4.u(3,:);
yd = Loop(i).RK4.u(4,:);
C_values = C(x, y, xd, yd);
plot(Loop(i).RK4.t, abs(C_values-C_values(1)),'b','LineWidth',2)
hold on

x = Reference.Loop(i).u(1,:);
y = Reference.Loop(i).u(2,:);
xd = Reference.Loop(i).u(3,:);
yd = Reference.Loop(i).u(4,:);
C_values = C(x, y, xd, yd);
plot(Reference.Loop(i).x, abs(C_values-C_values(1)),'g','LineWidth',2)

x = Loop(i).RKF45.u(1,:);
y = Loop(i).RKF45.u(2,:);
xd = Loop(i).RKF45.u(3,:);
yd = Loop(i).RKF45.u(4,:);
C_values = C(x, y, xd, yd);
plot(Loop(i).RKF45.t, abs(C_values-C_values(1)),'r','LineWidth',2)
%ylim([0,1e-7])
xlabel('t [s]')
ylabel('$$ \mathrm{C_{diff}} [-]$$')
legend 'RK4' 'ode45' 'RKF45'
end
set(findall(fig(7),'-property','FontSize'),'FontSize', 16)
set(findall(fig(7),'-property','Box'),'Box', 'on')
set(findall(fig(7),'-property','Interpreter'),'Interpreter', 'latex')
set(findall(fig(7),'-property','TickLabelInterpreter'),'TickLabelInterpreter', 'latex')

% exportgraphics(fig(4), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\Orbit_TwoLoop.png', 'ContentType', 'vector');
% exportgraphics(fig(5), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\Orbit_ThreeLoop.png', 'ContentType', 'vector');
% exportgraphics(fig(6), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\Orbit_FourLoop.png', 'ContentType', 'vector');
exportgraphics(fig(7), 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\JacobiConst.png', 'ContentType', 'vector');