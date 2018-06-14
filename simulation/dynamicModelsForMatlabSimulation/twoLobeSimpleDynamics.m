function dz = twoLobeSimpleDynamics(t, z, param)
%
%   This function computes the dynamics of a two lobed, foam robot that has
%   been simplified as two point masses with oscilating mass values between
%   1 and 2 kgs. 
%
%   INPUTS:
%       t = [1, 1] = time scalar
%       z = [6; 1] = [position; velocity] = current state of robot
%       param = struct = parameters for system
%           .massA = function handle = mA = massA(t)
%           .massB = function handle = mB = massB(t)
%           .motor = scalar = max force of oscilating motor
%           .mus = scalar = coefficient of static friction
%           .muk = scalar = coefficient of kinetic friction
%           .k  = scalar = spring constant
%           .lc = scalar = uncompressed spring length (i.e. the distance
%                         the "points" of contact on robot)
%           .g  = scalar = gravity constant
%
%   OUTPUTS:
%       dz = [6, nTime] = [velocity; acceleration] = time derivative of
%       current state
%
%   NOTES:
%       Dynamics model was simplified from a two lobed robot to a simple
%       two point mass system with a spring modeling internal forces of
%       compression and expansion during locomotion.
%       
%       5/24/18 - will first implement this to only take scalars only to
%       work out the dynamics. Eventually will (fully) vectorize code.
%
%       6/14/18 - updated comments - output of function is a [6, nTime]
%       vector of the time derivative of currrent state.
%
%   By: Martin Majkut
%   Date: May 24, 2018
%

% Pull out velocities
iVel = 4:6;                 % indecies of velocities
vel  = z(iVel, :);          % velocities

% Calculate accelerations
delX = abs( z(1,:) - z(2,:) );  % distances between point masses
mA   = param.massA(t);          % vector of masses for point A
mB   = param.massB(t);          % vector of masses for point B
mCOM = mA + mB;

if (param.k*(param.lc - delX)) > (param.mus * mA * param.g)
    aA = param.muk * param.g + (param.k * (param.lc - delX)) ./ mA;
else
    aA = 0;
end

fMotor = param.motor * cos(t);

if (param.k*(param.lc - delX) + fMotor) > (param.mus * mB * param.g)
    aB = param.muk * param.g - (fMotor + param.k * (param.lc - delX)) ./ mB;
else
    aB = 0;
end

aCOM   = ( mA.*vel(1,:) + mB.*vel(2,:) ) ./ mCOM;
accel  = [ aA; aB; aCOM ];

% Form output
dz = [ vel; accel ];

end 
