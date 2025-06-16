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

Tol =  @(i) 1e-4./10.^((i-1)/4);
Ni = 17;
Nj = 17;

RKF45Conv.TwoLoop.endpoints = zeros(4,Ni,Nj);
RKF45Conv.ThreeLoop.endpoints = zeros(4,Ni,Nj);
RKF45Conv.FourLoop.endpoints = zeros(4,Ni,Nj);

RKF45Conv.TwoLoop.nsteps = zeros(Ni,Nj);
RKF45Conv.ThreeLoop.nsteps = zeros(Ni,Nj);
RKF45Conv.FourLoop.nsteps = zeros(Ni,Nj);

RKF45Conv.TwoLoop.nfailed = zeros(Ni,Nj);
RKF45Conv.ThreeLoop.nfailed = zeros(Ni,Nj);
RKF45Conv.FourLoop.nfailed = zeros(Ni,Nj);

RKF45Conv.TwoLoop.nfevals = zeros(Ni,Nj);
RKF45Conv.ThreeLoop.nfevals = zeros(Ni,Nj);
RKF45Conv.FourLoop.nfevals = zeros(Ni,Nj);
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

        RKF45Conv.TwoLoop.nsteps(i,j) = TwoLoop.RKF45.stats.nsteps;
        RKF45Conv.ThreeLoop.nsteps(i,j) = ThreeLoop.RKF45.stats.nsteps;
        RKF45Conv.FourLoop.nsteps(i,j) = FourLoop.RKF45.stats.nsteps;

        RKF45Conv.TwoLoop.nfailed(i,j) = TwoLoop.RKF45.stats.nfailed;
        RKF45Conv.ThreeLoop.nfailed(i,j) = ThreeLoop.RKF45.stats.nfailed;
        RKF45Conv.FourLoop.nfailed(i,j) = FourLoop.RKF45.stats.nfailed;

        RKF45Conv.TwoLoop.nfevals(i,j) = TwoLoop.RKF45.stats.nfevals;
        RKF45Conv.ThreeLoop.nfevals(i,j) = ThreeLoop.RKF45.stats.nfevals;
        RKF45Conv.FourLoop.nfevals(i,j) = FourLoop.RKF45.stats.nfevals;
        clear TwoLoop ThreeLoop FourLoop
    end
end

fig1 = figure(1);
fig1.Position = [100 100 1600 800];
tiledlayout(2,2,TileSpacing="tight")
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(1,:,:))')
fontsize(12,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_1 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(2,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_2 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(3,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_3 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(4,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_4 [-]$$')


fig2 = figure(2);
fig2.Position = [100 100 1600 800];
tiledlayout(2,2,TileSpacing="tight")
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(1,:,:))')
fontsize(12,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_1 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(2,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_2 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(3,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_3 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(4,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_4 [-]$$')

fig3 = figure(3);
fig3.Position = [100 100 1600 800];
tiledlayout(2,2,TileSpacing="tight")
nexttile
surf(squeeze(log10(Tol(1:Ni))),squeeze(log10(Tol(1:Nj))),squeeze(RKF45Conv.FourLoop.endpoints(1,:,:))')
fontsize(12,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_1 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.FourLoop.endpoints(2,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_2 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.FourLoop.endpoints(3,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_3 [-]$$')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),squeeze(RKF45Conv.FourLoop.endpoints(4,:,:))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_4 [-]$$')

% zoom
fig4 = figure(4);
fig4.Position = [100 100 1600 800];
tiledlayout(2,2,TileSpacing="tight")
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(1,6:end,6:end))')
fontsize(12,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_1 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(2,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_2 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(3,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_3 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.TwoLoop.endpoints(4,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_4 [-]$$')


fig5 = figure(5);
fig5.Position = [100 100 1600 800];
tiledlayout(2,2,TileSpacing="tight")
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(1,6:end,6:end))')
fontsize(12,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_1 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(2,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_2 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(3,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_3 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.ThreeLoop.endpoints(4,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_4 [-]$$')


fig6 = figure(6);
fig6.Position = [100 100 1600 800];
tiledlayout(2,2,TileSpacing="tight")
nexttile
surf(squeeze(log10(Tol(6:Ni))),squeeze(log10(Tol(6:Nj))),squeeze(RKF45Conv.FourLoop.endpoints(1,6:end,6:end))')
fontsize(12,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_1 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.FourLoop.endpoints(2,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_2 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.FourLoop.endpoints(3,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_3 [-]$$')
nexttile
surf(log10(Tol(6:Ni)),log10(Tol(6:Nj)),squeeze(RKF45Conv.FourLoop.endpoints(4,6:end,6:end))')
fontsize(12,"points")
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
zlabel('$$u_4 [-]$$')

fig7 = figure(7);
fig7.Position = [100 100 1600 600];
tiledlayout(1,3,TileSpacing="tight")
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.TwoLoop.nsteps(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Steps')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.TwoLoop.nfailed(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel(' lg(ATOL) [-]')
ylabel(' lg(RTOL) [-]')
title('Failed steps')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.TwoLoop.nfevals(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Function evaluations')

fig8 = figure(8);
fig8.Position = [100 100 1600 600];
tiledlayout(1,3,TileSpacing="tight")
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.ThreeLoop.nsteps(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Steps')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.ThreeLoop.nfailed(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Failed steps')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.ThreeLoop.nfevals(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel(' lg(ATOL) [-]')
ylabel(' lg(RTOL) [-]')
title('Function evaluations')

fig9 = figure(9);
fig9.Position = [100 100 1600 600];
tiledlayout(1,3,TileSpacing="tight")
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.FourLoop.nsteps(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Steps')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.FourLoop.nfailed(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Failed steps')
nexttile
surf(log10(Tol(1:Ni)),log10(Tol(1:Nj)),RKF45Conv.FourLoop.nfevals(:,:)')
fontsize(15,"points")
set(0,'defaulttextinterpreter','latex')
set(gca, 'TickLabelInterpreter','latex')
xlabel('lg(ATOL) [-]')
ylabel('lg(RTOL) [-]')
title('Function evaluations')


exportgraphics(fig1, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Convergence_TwoLoop.png', 'ContentType', 'vector');
exportgraphics(fig2, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Convergence_ThreeLoop.png', 'ContentType', 'vector');
exportgraphics(fig3, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Convergence_FourLoop.png', 'ContentType', 'vector');

exportgraphics(fig4, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Convergence_TwoLoop_zoom.png', 'ContentType', 'vector');
exportgraphics(fig5, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Convergence_ThreeLoop_zoom.png', 'ContentType', 'vector');
exportgraphics(fig6, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Convergence_FourLoop_zoom.png', 'ContentType', 'vector');

exportgraphics(fig7, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Stats_TwoLoop.png', 'ContentType', 'vector');
exportgraphics(fig8, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Stats_ThreeLoop.png', 'ContentType', 'vector');
exportgraphics(fig9, 'Y:\Egyetem\MSc\1Semester\Math\project\Three-Body-Problem\figures\RKF45_Stats_FourLoop.png', 'ContentType', 'vector');

