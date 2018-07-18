function dz = inclinedPlaneDynamics(z, param)
%   dz = inclinedPlaneDynamics(z,param)
%
%   This function coputes the dynamics of a point mass / box on an inclined
%   plane. The plane can have variable coefficient of friction or a
%   constant coefficient of friction.
%
%   INPUTS:
%       z = [2, 1] = [position; velocity] = current state of mass
%       param = struct = parameters for system
%           .theta = scalar = angle between the plane and horizontal axis
%           .mus = function handle = coefficient of static friction 
%               mus = mus(x)
%           .muk = function handle = coefficient of kinetic friction
%               muk = muk(x)
%           .g = scalar = gravity constant
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

% Pull out velocity and initialize acceleration
pos = z(1);
vel = z(2);
acc = 0;

% Determine acceleration 
if vel ~= 0
    acc = param.g*(sin(param.theta) - param.muk(pos)*cos(param.theta))*sign(vel);
else
    if sin(param.theta) > param.mus(pos)*cos(param.theta)
        acc = param.g*(sin(param.theta) - param.muk(pos)*cos(param.theta)*sign(vel));
    else
        acc=0;
        vel=0;
    end
end

% Form output
dz = [vel; acc];

end
