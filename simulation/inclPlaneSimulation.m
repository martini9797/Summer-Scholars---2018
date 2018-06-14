function inclPlaneSimulation(isFrictionVariable)
%
%   This function runs a simulation of a point mass / box on an inclined
%   plane. Two dynamic models are available: constant friction inclined
%   plane or variable friction plane.
%
%   INPUTS:
%       (optional) isFrictionVariable = bool = if true, will run simulation
%           with a variable friction plane, where there are "patches" of 
%           friction with different friction coefficients. Otherwise, will 
%           run simulation with constant friction plane.
%
%   NOTES:
%       Default is to run constant friction model.
%
%   By: Martin Majkut
%   Date: Jun. 14, 2018
%

% Clear workspace and add dynaics to path (a clear workspace is nice)
close all; clc;
addpath('dynamicModelsForMatlabSimulation/');

% Check for which dynamics model to run
if nargin < 1
    isFrictionVariable = false;
end

% Set initial conditions
tStart = 0;    % [s]
tEnd   = 1;    % [s]
nGrid  = 100;
tGrid  = linspace(tStart, tEnd, nGrid);

param.theta = pi/6;
param.mus   = 0.5;
param.muk   = 0.2;
param.g     = 9.81;

zInit = zeros(2,1);

if isFrictionVariable
    error('TODO: Implement the variable friction dynamics model :)');
else
    dynFun = @(t,z)( inclinedPlaneDynamics(z, param) );
end

% Run simulation
[zGrid, ~] = runSimulation(dynFun, tGrid, zInit, 'rk4');

% Plot Results
figure(61418); clf;
hold on;
title('Simulated Displacment of Point Mass on an Inclined Plane');

subplot(2,1,1);
hold on;
plot(tGrid, zGrid(1,:), 'k');
title('Position vs. Time');
xlabel('Time [s]');
ylabel('Position [m]');

subplot(2,1,2)
hold on;
plot(tGrid, zGrid(2,:), 'k');
title('Velocity vs. Time');
xlabel('Time [s]');
ylabel('Velocity [m/s]');

end