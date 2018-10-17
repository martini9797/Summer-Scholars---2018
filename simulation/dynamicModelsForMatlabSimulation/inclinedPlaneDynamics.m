function dz = inclinedPlaneDynamics(z, param)
%   dz = inclinedPlaneDynamics(z,param)
%
%   This function coputes the dynamics of a point mass / box on an inclined
%   plane. The plane can have variable coefficient of friction or a
%   constant coefficient of friction.
%
%   INPUTS:
%       z = [2, 1] = [position; velocity] = current state of mass
%       param = struct = parameters for system and it's discretization 
%           .theta = scalar = angle between the plane and horizontal axis
%           .mus   = function handle = coefficient of static friction 
%                    mus = mus(x)
%           .muk   = function handle = coefficient of kinetic friction
%                    muk = muk(x)
%           .g     = scalar = gravity constant
%           .tol   = (optional) = scalar = tolerance for determining when 
%                    the problem is a statics vs. dynamics problem 
%                    default = 1e-08
%
%   OUTPUTS:
%       dz = [2, 1] = [velocity; acceleration] = time derivative of current
%       state
%
%   NOTES:
%       Since friction is velocity dependent and that is not predefined but
%       rather calculatated, dynamics model cannot be vectorized
%
%       Mass parameter not needed for dynamics since analysis simplifies to
%       not include mass when calculating acceleration. This model is also
%       not time dependant, so t parameter not required.
%
%   By: Martin Majkut
%   Date: Jun. 14, 2018
%

% Pull out states and initizlize acceleration
pos = z(1);
vel = z(2);
acc = 0;

% Pull out parameters to make equations a little more readable
theta = param.theta;
mus   = param.mus;
muk   = param.muk;
g     = param.g;

if ~isfield(param,'tol')
    tol = 1e-08;
else
    tol = param.tol;
end

% Determine acceleration 
if vel > tol            % Dynamics problem
    acc = g * ( sin(theta) - muk(pos)*cos(theta)*sign(vel) );
elseif vel < tol        % If crossed zero (tolerance within 0), stopped moving and turns into statics problem
    vel = 0;
    if sin(theta) > mus(pos)*cos(theta)
        acc = g * (sin(theta) - mus(pos)*cos(theta));
    end
end

% Form output
dz = [vel; acc];

end
