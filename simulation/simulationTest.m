function simulationTest()
%
%   simulationTest()
%
%   This function runs a simulation testing the newly written dynamics 
%   model written to describe the dynamics of a two lobed foam robot
%
%   By: Martin Majkut
%   Date: May 24, 2018
%

% Clear workspace and add dynmics to path
close all; clear; clc;
addpath('dynamicModelsForMatlabSimulation/');

% Set initial conditions and parameters
tStart = 0;     % [s]
tEnd   = 50;    % [s]
nGrid  = 150;   
tGrid  = linspace(tStart, tEnd, nGrid);

phaseShift  = 0;        
% Note: phaseShift must be the same for both to keep total mass constant
param.massA = @(t) ( cos(t + pi/4 + phaseShift).^2 + 1);
param.massB = @(t) ( sin(t + pi/4 + phaseShift).^2 + 1);
param.motor = 10;       % [N]
param.mus   = 0.5;
param.muk   = 0.1;
param.k     = 10;       % [N/m]
param.lc    = 0.1;      % [m]
param.g     = 10;       % [m/s^2]

zInit    = zeros(6,1);    
zInit(2) = param.lc;

dynFun = @(t,z) (twoLobeSimpleDynamics(t, z, param));

% Run Simulation
[zGrid, ~] = runSimulation(dynFun, tGrid, zInit, 'rk4');

% Plot Results
figure(52418); clf;

subplot(2,2,1);
hold on;
plot(tGrid, zGrid(1,:), 'r', tGrid, zGrid(2,:), 'b');
title('Postion of Lobes vs. Time');
xlabel('Time [s]');
ylabel('Position [m]');
legend('A', 'B');

subplot(2,2,2);
hold on;
plot(tGrid,zGrid(3,:), 'k');
title('Position of Robot COM vs. Time');
xlabel('Time [s]');
ylabel('Position [m]');
legend('Robot COM');

subplot(2,2,3);
hold on;
plot(tGrid, zGrid(4,:), 'r', tGrid, zGrid(5,:), 'b');
title('Velocity of Lobes vs. Time');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
legend('A', 'B');

subplot(2,2,4);
hold on;
plot(tGrid,zGrid(6,:), 'k');
title('Velocity of Robot COM vs. Time');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
legend('Robot COM');

end
