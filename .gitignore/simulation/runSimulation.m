function [zGrid, nEval] = runSimulation(dynFun, tGrid, zInit, method)
%  [zGrid, nEval] = runSimulation(dynFun, tGrid, zInit, method)
%
% Simulate a dynamical system, given an integration routine, time grid,
% and initial state. Select from several explicit Runge-Kutta methods.
%
% INPUTS:
%    dynFun = a function handle:  dz = dynFun(t, z)
%        IN:  t = [1, nTime] = row vector of time
%        IN:  z = [nState, nTime] = matrix of states corresponding to each time
%        OUT: dz = [nState, nTime] = time-derivative of the state at each point
%    tGrid = [1, nGrid] = time grid to evaluate the simulation
%    zInit = [nDim, 1] = initial state
%    method = string = name of the desired method
%      'euler' = Euler's method (first-order)
%      'heun' = Heun's method (second-order)
%      'midpoint' = the midpoint method (second-order)
%      'ralston' = Ralston's method (second-order)
%      'rk4' = "The" Runge--Kutta method (forth-order)
%
% OUTPUTS:
%   zGrid = [nDim, nGrid] = state at each point in tGrid. zGrid(:,1) = zInit
%   nEval = scalar = total number of calls to the dynamics function
%
% Modified By: Martin Majkut
% Date:        Feb. 8, 2018

% Simulate based on three different methods
switch lower(method)
    
    case 'euler' 
        [zGrid, nEval] = EulerMethod(dynFun, tGrid, zInit);
        
    case 'midpoint'
        [zGrid, nEval] = MidpointMethod(dynFun, tGrid, zInit);
        
    case 'rk4'
        [zGrid, nEval] = RK4Method(dynFun, tGrid, zInit);
    
    % throw error if none of above method selected    
    otherwise
        error('Inputted method not implemented!');
        
end

end %%%%% End runSimulation.m %%%%%

%~~~~~~~~~~~~~~~~~~~~~ Function Definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [zGrid, nEval] = EulerMethod(dynFun, tGrid, zInit)
% function [zGrid, nEval] = EulerMethod(dynFun, tGrid, zInit)
%
% Simulate a dynamical system, given an integration routine, time grid,
% and initial state using the Forward Euler Method
% 
% INPUTS:
%    dynFun = a function handle:  dz = dynFun(t, z)
%        IN:  t = [1, nTime] = row vector of time
%        IN:  z = [nState, nTime] = matrix of states corresponding to each time
%        OUT: dz = [nState, nTime] = time-derivative of the state at each point
%    tGrid = [1, nGrid] = time grid to evaluate the simulation
%    zInit = [nDim, 1] = initial state
%
% OUTPUTS:
%   zGrid = [nDim, nGrid] = state at each point in tGrid. zGrid(:,1) = zInit
%   nEval = scalar = total number of calls to the dynamics function
%

% Set up for Euler's method
nTime  = length(tGrid);
nState = size(zInit,1);
zGrid  = zeros(nState, nTime);

% Set initial conditions
zGrid(:,1) = zInit;

% Calculate the next step
for kCurr = 1:(nTime - 1)
   
    % Unpack to make lines shorter
    kNext = kCurr + 1;
    zCurr = zGrid(:,kCurr);
    tCurr = tGrid(kCurr);
    h     = tGrid(kNext) - tGrid(kCurr);
    
    zGrid(:,kNext) = zCurr + h * dynFun(tCurr,zCurr);

end

% Calculate number of evaluations based on how many times dynFunc is called
nEval = nTime - 1;

end %%%%% End EulerMethod %%%%%

function [zGrid, nEval] = MidpointMethod(dynFun, tGrid, zInit)
% function [zGrid, nEval] = MidpointMethod(dynFun, tGrid, zInit)
%
% Simulate a dynamical system, given an integration routine, time grid,
% and initial state using the Midpoint Method
% 
% INPUTS:
%    dynFun = a function handle:  dz = dynFun(t, z)
%        IN:  t = [1, nTime] = row vector of time
%        IN:  z = [nState, nTime] = matrix of states corresponding to each time
%        OUT: dz = [nState, nTime] = time-derivative of the state at each point
%    tGrid = [1, nGrid] = time grid to evaluate the simulation
%    zInit = [nDim, 1] = initial state
%
% OUTPUTS:
%   zGrid = [nDim, nGrid] = state at each point in tGrid. zGrid(:,1) = zInit
%   nEval = scalar = total number of calls to the dynamics function
%

% Set up for Euler's method
nTime  = length(tGrid);
nState = size(zInit,1);
zGrid  = zeros(nState, nTime);

% Set initial conditions
zGrid(:,1) = zInit;

% Calculate the next step
for kCurr = 1:(nTime - 1)
   
    % Unpack to make lines shorter
    kNext = kCurr + 1;
    zCurr = zGrid(:,kCurr);
    tCurr = tGrid(kCurr);
    h     = tGrid(kNext) - tGrid(kCurr);
    
    zGrid(:,kNext) = zCurr + h * dynFun(tCurr + h/2, zCurr + (h/2) * ...
                                        dynFun(tCurr, zCurr));

end

% Calculate number of evaluations based on how many times dynFunc is called
nEval = 2 * (nTime - 1);

end %%%%% End MidpointMethod %%%%%%

function [zGrid, nEval] = RK4Method(dynFun, tGrid, zInit)
% function [zGrid, nEval] = RK4Method(dynFun, tGrid, zInit)
%
% Simulate a dynamical system, given an integration routine, time grid,
% and initial state using "The" Runge-Kutta Method
% 
% INPUTS:
%    dynFun = a function handle:  dz = dynFun(t, z)
%        IN:  t = [1, nTime] = row vector of time
%        IN:  z = [nState, nTime] = matrix of states corresponding to each time
%        OUT: dz = [nState, nTime] = time-derivative of the state at each point
%    tGrid = [1, nGrid] = time grid to evaluate the simulation
%    zInit = [nDim, 1] = initial state
%
% OUTPUTS:
%   zGrid = [nDim, nGrid] = state at each point in tGrid. zGrid(:,1) = zInit
%   nEval = scalar = total number of calls to the dynamics function
%

% Set up for Euler's method
nTime  = length(tGrid);
nState = size(zInit,1);
zGrid  = zeros(nState, nTime);

% Set initial conditions
zGrid(:,1) = zInit;

% Calculate the next step
for kCurr = 1:(nTime - 1)
   
    % Unpack to make lines shorter
    kNext = kCurr + 1;
    zCurr = zGrid(:,kCurr);
    tCurr = tGrid(kCurr);
    h     = tGrid(kNext) - tGrid(kCurr);
    
    % Calculate the intermediate steps
    k1 = dynFun(tCurr, zCurr);
    k2 = dynFun(tCurr + h/2, zCurr + (h/2)*k1);
    k3 = dynFun(tCurr + h/2, zCurr + (h/2)*k2);
    k4 = dynFun(tCurr + h  , zCurr + h * k3);
    
    % Calculate the full step
    zGrid(:,kNext) = zCurr + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end

% Calculate number of evaluations based on how many times dynFunc is called
nEval = 4 * (nTime - 1);

end %%%%% End RK4Method %%%%%%
