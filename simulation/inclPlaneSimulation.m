function inclPlaneSimulation(isFrictionVariable)
%
%   This function runs a simulation of a point mass / box on an inclined
%   plane. Two dynamic models are available: constant friction inclined
%   plane or variable friction plane.
%
%   INPUTS:
%       isFrictionVariable (optional) = bool = if true, will run simulation
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

% Add dynaics to path
addpath('dynamicModelsForMatlabSimulation/');

% Check for which dynamics model to run
if nargin < 1
    isFrictionVariable = false;
end

% Set initial conditions, parameters for dynamics and simulation
tStart = 0;     % [s]
tEnd   = 5;     % [s]
nGrid  = 1000;
tGrid  = linspace(tStart, tEnd, nGrid);

param.theta = pi/6;
param.g     = 9.81;

zInit = zeros(2,1);

% Set friction model
mus      = 0.4;
muk      = 0.1;
patchLen = 10;  %[m]
if isFrictionVariable
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

% Animation settings
boxH   = 0.05 * zGrid(1,end);   %[m]
boxW   = 0.05 * zGrid(1,end);   %[m]
boxCol = [48 15 0] ./ 255;  % rgb, Tufts brown for a brown box (:
%dt     = tGrid(2) - tGrid(1);   %[s]

% Animate Simulation
hFig = figure(62818); clf;
hold on;
title('Animation of Sliding Box');
daspect([1 1 1]);
set(gca,'xtick',[],'ytick',[]);

% Create rotation matrix to convert from local to global coordinate system
r = [cos(param.theta) -sin(param.theta);
     sin(param.theta)  cos(param.theta)];
 
% Draw Inclined Plane
xFlatStart = zGrid(1,1);

if isFrictionVariable
    
    % Determine the vertices of each patch along the plane
    numPatches = round(zGrid(1,end)/patchLen);
    xFlatEnd   = numPatches * patchLen;
    xFlat      = xFlatStart:patchLen:xFlatEnd;
    xMat       = [ xFlat(1:end-1); xFlat(2:end); xFlat(2:end); ...
                   xFlat(1:end-1); xFlat(1:end-1)];
    planeThick = 0.01 * xFlatEnd;      %[m]
    yFlat      = planeThick * [0 0 -1 -1 0];
    yMat       = repmat(yFlat,1,numPatches);
    
    % Reshape into coumn vectors to do math to rotate plane
    xMat       = reshape(xMat,[],1);
    yMat       = reshape(yMat,[],1);
    planeCoor  = [xMat yMat] * r;
    
    % Reshape back into matrices that patch() accepts
    xMat       = reshape(planeCoor(:,1),[],numPatches);
    yMat       = reshape(planeCoor(:,2),[],numPatches);
    patch(xMat,yMat,rand(numPatches,1));
    axis( [ xMat(1,1) xMat(3,end) yMat(2,end) yMat(1,1)+boxH ] );

else
    
    %Determine the vertices of the plane, rotate and draw it
    xFlatEnd   = zGrid(1,end);
    xFlat      = [xFlatStart xFlatEnd xFlatEnd xFlatStart]';
    planeThick = 0.0025 * xFlatEnd;      %[m]
    yFlat      = planeThick * [0 0 -1 -1]';
    planeCoor  = [xFlat yFlat] * r;
    patch(planeCoor(:,1),planeCoor(:,2),'k');
    axis( [ planeCoor(1,1) planeCoor(2,1) planeCoor(end - 1,2) ...
                 planeCoor(1,2)+boxH ] );

end

% Determine vertices of the box as it slides down the ramp
xBox    = [ zGrid(1,:); zGrid(1,:); zGrid(1,:) + boxW; zGrid(1,:) + boxW ];
yBox    = repmat([0; boxH; boxH; 0],1,size(zGrid,2));
boxCoor = [reshape(xBox,[],1) reshape(yBox,[],1)] * r;
xBox    = reshape(boxCoor(:,1),[],size(zGrid,2));
yBox    = reshape(boxCoor(:,2),[],size(zGrid,2));

% Animate
startT = tic;
box    = patch(xBox(:,1),yBox(:,1),boxCol);
currT  = toc(startT);
while currT < tGrid(end)
    
    delete(box);
    
    xNow = interp1(tGrid,xBox',currT,'linear');
    yNow = interp1(tGrid,yBox',currT,'linear');
    
    box  = patch(xNow',yNow',boxCol);
    drawnow;
    
    currT = toc(startT);
end

end