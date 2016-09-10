function uf = circlegaussianrelation(u0,ds)
% This will take an initial condition and evolve the ODE for the scattering
% relation. It will also make sure to get when it exits, and return the
% exit position and velocity. 
% Maybe use a fixed step so that way it's more modulated? 

options = odeset('Events',@cgEventsFcn); 

[~,u,~,ue] = ode45(@gaussianmetric, [0,ds] , u0, options); 
% [~,u] = ode45(@gaussianmetric, [0,1], u0);  % for example
% it should be kept adapative for the interval of length

% uf = u0;   % initialize it, but I don't think we need to?

% When the ray has not left the domain yet
if isempty(ue)==1
    u0 = u(end,:);
    uf = squaregaussianrelation(u0,ds);
    return;
end

% When the ray does leave/ODE evolution stops
if isempty(ue)==0
   uf = ue;
   % Renormalize the velocity (since we are interested in direction):
   uf(3) = uf(3)/sqrt(uf(3)^2 + uf(4)^2);
   uf(4) = uf(4)/sqrt(uf(3)^2 + uf(4)^2);
end

end
