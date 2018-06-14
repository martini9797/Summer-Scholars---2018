function dz = constFrictionInclPlaneSimpleDynamics(z, param)
%
%   This function coputes the dynamics of a point mass / box on an inclined
%   plane with a constant friction.
%
%   INPUTS:
%       z = [2, 1] = [position; velocity] = current state of mass
%       param = struct = parameters for system
%           .theta = scalar = angle between the plane and horizontal axis
%           .mus = scalar = coefficient of static friction
%           .muk = scalar = coefficient of kinetic friction
%           .g = scalar = gravity constant
%
%   OUTPUTS:
%       dz = [2, 1] = [velocity; acceleration] = time derivative of current
%       state
%
%   NOTES:
%       Mass parameter not needed for dynamics since analysis simplifies to
%       not include mass when calculating acceleration. This model is also
%       not time dependant, so t parameter not required.
%
%   By: Martin Majkut
%   Date: Jun. 14, 2018
%

% Pull out velocity
vel = z(2);
acc = 0;

% Determine acceleration 
if vel ~= 0
    acc = param.g*(sin(param.theta) - param.muk*cos(param.theta));
else
    if sin(param.theta) > param.mus*cos(param.theta)
        acc = param.g*(sin(param.theta) - param.muk*cos(param.theta));
    end
end

% Form output
dz = [vel; acc];

end
    
    
    
    