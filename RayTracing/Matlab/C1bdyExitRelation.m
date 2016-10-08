function uf = C1bdyExitRelation(u0,ds,metric,options)
% This will take an initial condition and evolve the ODE for the scattering
% relation. It will also make sure to get when it exits, and return the
% exit position and velocity. 
% Designed so that it can take any boundary, given the options in events.

[~,u,~,ue] = ode45(metric, [0,ds] , u0, options);
uf = u(end,:);
% [~,u] = ode45(@gaussianmetric, [0,1], u0);  % for example
% it should be kept adapative for the interval of length

% uf = u0;   % initialize it, but I don't think we need to?

% When the ray has not left the domain yet:
% Need to be careful! If it is trapped? Maybe instead, use a while loop and
% have a stopping # of iterations (and return like all NaNs)? 
iter = 1;
scale = max(ds,1/ds);
while isempty(ue)==1 && iter < 100*scale
    u0 = u(end,:);
    [~,u,~,ue] = ode45(metric, [0,ds] , u0, options);
    iter = iter + 1;
    % uf = C1bdyExitRelation(u0,ds,metric,options);
end

if isempty(ue)==0
    uf=ue(1,:);   % Saw a bug where ue was matrix and not vector...??
    % uf(3) = mod(uf(3),2*pi);  % Do we want to mod by 2*pi?
    return;
end

if isempty(ue) == 1  % Suspect trapped ray then.
   uf = u(end,:);
   i = 1:length(u0);
   uf(i) = Inf;
   return;
end

% % When the ray does leave/ODE evolution stops
% if isempty(ue)==0
%    uf = ue;
%    % Renormalize the velocity (since we are interested in direction):
%    uf(3) = uf(3)/sqrt(uf(3)^2 + uf(4)^2);
%    uf(4) = uf(4)/sqrt(uf(3)^2 + uf(4)^2);
% end

end
