function uf = squaregaussianrelation(u0,ds)
% This will take an initial condition and evolve the ODE for the scattering
% relation. It will also make sure to get when it exits, and return the
% exit position and velocity. 
% Maybe use a fixed step so that way it's more modulated? 

options = odeset('Events',@sgEventsFcn); 

[~,u,~,ue] = ode45(@gaussianmetric, [0,ds] , u0, options); 
% [~,u] = ode45(@gaussianmetric, [0,1], u0);  % for example
% it should be kept adapative for the interval of length

% uf = u0;   % initialize it, but I don't think we need to?
if isempty(ue) ==1
    u0 = u(end,:);
    uf = squaregaussianrelation(u0,ds);
    return;
end

if isempty(ue)==0
   uf = ue;
   % Renormalize the velocity (since we are interested in direction):
   uf(3) = uf(3)/sqrt(uf(3)^2 + uf(4)^2);
   uf(4) = uf(4)/sqrt(uf(3)^2 + uf(4)^2);
end

end







%%%%% Old one with non-adaptive stepping: 

% % Goes 2:end in case for general boundary, roundoff error?
% kx = find(abs(u(2:end,1)) > 1, 1, 'first');
% ky = find(abs(u(2:end,2)) > 1, 1, 'first');
% 
% if isempty(kx) && isempty(ky)
%     u0 = u(end,:);
%     uf = squaregaussianrelation(u0,ds);
%     return;
% end
% 
% if isempty(kx)
%    kx = ky+1; 
% end
% if isempty(ky)
%    ky = kx+1;
% end

% % See if the ray has left the boundary. We start from the first index 
% % because the ray could possibly re-enter the domain. 
% % Trying to improve prior interpolation method:
% k = min(kx,ky) + 1;
% uf = u(k,:);
% dy = u(k,2)-u(k-1,2);
% dx = u(k,1)-u(k-1,1);
% m = dy/dx;
% % Can we just use slope m (or, dy and dx) for the directions/velocities?
% 
% if abs(uf(1)) > 1
%     uf(1) = uf(1)/abs(uf(1));   % So it preserves the sign.
%     
%     % Interpolated y component
%     uf(2) = u(k-1,2) + m*(uf(1)-u(k-1,1));
%     
%     % Interpolated velocities
%     uf(3) = dx/(dy^2+dx^2)^0.5;
%     uf(4) = dy/(dy^2+dx^2)^0.5;
%     
%     % To check difference between two interpolation approaches:
%     % Use r to interpolate exit velocities with weighted averaging:
%     % r = (uf(1)-u(k-1,1))/(u(k,1)-u(k-1,1));
%     % uf(3) = uf(3) - (u(k-1,3) + r*(u(k,3) - u(k-1,3)));
%     % uf(4) = uf(4) - (u(k-1,4) + r*(u(k,4) - u(k-1,4)));
% 
%     return;
% end
% 
% if abs(uf(2)) > 1
%     uf(2) = uf(2)/abs(uf(2));   % So it preserves the sign.
%     
%     % Interpolated x component
%     uf(1) = u(k-1,1) + (uf(2)-u(k-1,2))/m;
%     
%     % Interpolated velocities
%     uf(3) = dx/(dy^2+dx^2)^0.5;
%     uf(4) = dy/(dy^2+dx^2)^0.5;
%     
%     % To check difference between two interpolation approaches:
%     % Use r to interpolate exit velocities with weighted averaging:
%     % r = (uf(1)-u(k-1,1))/(u(k,1)-u(k-1,1));
%     % uf(3) = uf(3) - (u(k-1,3) + r*(u(k,3) - u(k-1,3)));
%     % uf(4) = uf(4) - (u(k-1,4) + r*(u(k,4) - u(k-1,4)));
% 
%     return;
% end


%%%%%% Even Older one: %%%%%%

% See if the ray has left the boundary. We start from s=0 because the ray
% could possibly re-enter the domain. 

% i=2;
% while i <= length(u(:,1))
%    if abs(u(i,1)) > 1  || abs(u(i,2)) > 1
%        % To linearly interpolate exit positions: 
%        m = (u(i,2)-u(i-1,2))/(u(i,1)-u(i-1,1));  
%        uf = u(i,:);
%        if abs(uf(1)) > 1
%           uf(1) = u(i,1)/abs(u(i,1));   % So it preserves the sign.
%           % To interpolate exit velocities.
%           r = (uf(1)-u(i-1,1))/(u(i,1)-u(i-1,1));  
%           % Interpolated y component 
%           uf(2) = u(i-1,2) + m*(uf(1)-u(i-1,1)); 
%           % Interpolated velocities
%           uf(3) = u(i-1,3) + r*(u(i,3) - u(i-1,3)); 
%           uf(4) = u(i-1,4) + r*(u(i,4) - u(i-1,4));
%           %  return;
%        end
%        
%        % Also will take care of the case if both |x|,|y| > 1
%        % If |y| > 1 still, means |y| hits 1 before |x| does...?
%        if abs(uf(2)) > 1
%           uf(2)=u(i,2)/abs(u(i,2));  % To preserve sign.
%           r = (uf(2)-u(i-1,2))/(u(i,2)-u(i-1,2));
%           % Interpolated x component
%           uf(1) = u(i-1,1) + (uf(2)-u(i-1,2))/m; 
%           % Interpolated velocities
%           uf(3) = u(i-1,3) + r*(u(i,3) - u(i-1,3));
%           uf(4) = u(i-1,4) + r*(u(i,4) - u(i-1,4));
%           %  return;
%        end
%        return;
%    end
%    i = i + 1;
% end
% 
% % If the method hasn't returned, we need to evolve the ODE more.
% u0 = u(end,:);
% uf = squaregaussianrelation(u0);

% end


%%%%%%%%%%%%
% Any More Scratch Work:


