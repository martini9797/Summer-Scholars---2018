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
%           run simulation with constant friction plane. Default = false
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

% Set initial conditions, parameters for dynamics and simulation
tStart = 0;     % [s]
tEnd   = 10;    % [s]
nGrid  = 1000;
tGrid  = linspace(tStart, tEnd, nGrid);

param.theta = pi/6;
param.g     = 9.81;

zInit = zeros(2,1);

% Set friction model
mus      = 0.5;
muk      = 0.2;
patchLen = 10;     %[m]
if isFrictionVariable
    % Uncomment for on/off friction model 
    param.mus = @(x) mus * (mod(round(x/patchLen),2) - 1);
    param.muk = @(x) muk * (mod(round(x/patchLen),2) - 1);
else
    param.mus = @(x) mus;
    param.muk = @(x) muk;
end

dynFun = @(t,z)( inclinedPlaneDynamics(z, param) );

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

% Animate and Record the Simulation
figure(62818); clf;
hold on;
title('Animation of Sliding Box');

% Draw Inclined Plane
planeThick = .01; %[m]
xFlatStart = zGrid(1,1);
r          = [cos(param.theta) -sin(param.theta); 
              sin(param.theta)  cos(param.theta)];

if isFrictionVariable
    xFlatEnd = round(zGrid(1,end)/patchLen) * patchLen;
else
    xFlatEnd = zGrid(1,end);
    xFlat    = [xFlatStart xFlatEnd xFlatEnd xFlatStart];
    yFlat    = planeThick * [0 0 -1 -1];
end

planeCoor = [xFlat; yFlat]' * r;
planeCoor = planeCoor';
axis([planeCoor(1,1) inf -inf planeCoor(2,end)]);
daspect([1 1 1]);
fill(planeCoor(1,:),planeCoor(2,:),'k');




end