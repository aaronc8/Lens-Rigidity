%% Scripts to run ode45 for ray tracing
% This is just for the square case and the metric g is a Gaussian. 
% First define the function:

% function dH = gaussianmetric(s,u)
% % u(1) = x coordinate, u(2) = y coordinate, 
% % u(3) = x velocity, u(4) = y velocity.
%     dH = zeros(4,1);
%     % The positions:
%     dH(1) = u(3).*exp(u(1).^2 + u(2).^2);
%     dH(2) = u(4).*exp(u(1).^2 + u(2).^2);
%     % The momenta: 
%     dH(3) = (u(3).^2 + u(4).^2).*u(1).*exp(u(1).^2 + u(2).^2);
%     dH(4) = (u(3).^2 + u(4).^2).*u(2).*exp(u(1).^2 + u(2).^2);
% end

%% Checking the ODE evolution is fine and iterating.
tic;
Nedge = 2^5;   % number of partitions of sides.
Nangle = 2^5;   % number of partitions of angle.
dphi = pi/Nangle;   % This needs to be taken better care of
dl = 2/Nedge; 
% Should double for loop over n,N!
% Can switch the order of for loops to see the change in ray traces with 
% entry position or entry angle graphically.
%%
% To terminate when either |x|,|y| >= 1
options = odeset('Events',@sgEventsFcn); 
%Here,
% function [value,isterminal,direction] = sgEventsFcn(s,u)
% % Terminate ODE evolution when the position is outside of square.
% value = max(0,1 - max(abs(u(1)), abs(u(2))));
% isterminal = 1;
% direction = -1;
% end
%%
for i = 1:Nedge-2
    clf
    for j = 1:Nangle-2
        % u0 = [-1;1-i*dl;cos(j*dphi-pi/2); sin(j*dphi-pi/2)];   % Only for left edge!
        % u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   % This is for bottom edge!
         u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   % Only for right edge!
        % u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   % This is for top edge!
        [s,u,~,ue] = ode45(@gaussianmetric, [0,1], u0, options);
        % Noting how some rays don't make it out yet
        plot(u(:,1), u(:,2)); axis([-2 2 -2 2]);
        xlabel('x position'); ylabel('y position');
        x = -1:0.1:1;
        y = ones(length(x),1);
        hold on
        plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);
        pause(0.01);
    end
end
toc;
% The last spread of angles is displayed

%% Scattering Relation, Exit Data. 
% Then search for when it FIRST leaves the box somehow. 
% Uses squaregaussianrelation.m --> if working correctly, then use to
% get all the data with SGscatteringrelation.m 

% Two velocity interpolations used in squaregaussianrelation.m , 
% Could subtract two methods to see how their accuracy is w.r.t. each other

tic;
ds = 1;  % Or make it smaller if desired?
uNSEW = cell((Nedge-2),(Nangle-2)); 
% The cell of all exit-pos/vel for each given incidence pos/angle. 
% To access: Use for ex. cellfun(@(v) v(4), uNSEW(:,1)) ?
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. So,
% The 1st entry ~ location on the edge. 2nd entry ~ angle of incidence.

for i = 1:Nedge-2
    for j = 1:Nangle-2
        % u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   % Only for left (West) edge!
        % u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   % This is for bottom (South) edge!
         u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   % Only for right (East) edge!
        % u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   % This is for top (North) edge!
        uNSEW{i,j} = squaregaussianrelation(u0,ds);
    end
end
toc;
%%
% To check qualitatively with the plot above:
celldisp(uNSEW(end,:))
% Be careful if we wanted to switch the order of for loops
% It looked consistent with all 4 sides. 
%% Entire scattering relation
% Use: (Takes a while to compute)
tic;
[uW, uN, uE, uS] = SGscatteringrelation(Nedge, Nangle);
toc;
% To check with prior notes I guess:
celldisp(uE(end,:))





