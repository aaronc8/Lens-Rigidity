%% Scripts to run ode45 for ray tracing
% This is just for the square case and the metric g is a Gaussian. 
% These are just notes / scratch work primarily. 

% function dH = gaussianmetrictheta(s,u)
% % u(1) = x coordinate, u(2) = y coordinate, 
% % u(3) = x velocity, u(4) = y velocity.
% % This has no bearing on the domain anyways!
% 
%     dH = zeros(3,1);
%     % The positions:
%     dH(1) = cos(u(3)).*exp(u(1).^2/2 + u(2).^2/2);
%     dH(2) = sin(u(3)).*exp(u(1).^2/2 + u(2).^2/2);
%     % The momenta: 
%     dH(3) = (-u(1).*sin(u(3))+u(2).*cos(u(3))).*exp(u(1).^2/2 + u(2).^2/2);
%     
% end

%% Checking the ODE evolution is fine and iterating.
tic;
Nedge = 2^5;   % number of partitions of sides.
Nangle = 2^5;   % number of partitions of angle.
dphi = pi/Nangle;   % This needs to be taken better care of
dl = 2/Nedge; 
% Two loops over n,N
% Can switch the order of for loops to see the change in ray traces with 
% entry position or entry angle graphically.
%%
% To terminate when either |x|,|y| >= 1
options = odeset('Events',@sgEventsFcn); 
% Here, the events function was defined as:
% function [value,isterminal,direction] = sgEventsFcn(s,u)
% % Terminate ODE evolution when the position is outside of square.
% value = max(0,1 - max(abs(u(1)), abs(u(2))));
% isterminal = 1;
% direction = -1;
% end
%%
ds=2;
for i = 1:Nedge-1
    clf
    % INCLUDE THIS: So that the event will trigger when the starting point
    % is close to the boundary: 
    % ds = dl*2*i*(Nedge-i)/Nedge;    % leave off to see rays fully exit.
    
    for j = 1:Nangle-1
        % u0 = [-1;1-i*dl;cos(j*dphi-pi/2); sin(j*dphi-pi/2)];   % Only for left edge!
        % u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   % This is for bottom edge!
        % u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   % Only for right edge!
         u0 = [1-i*dl;1;-j*dphi];   % This is for top edge!
        [s,u,~,ue] = ode45(@gaussianmetrictheta, [0,ds], u0, options);
        % The arclength interval of ODE evolution has to be small enough
        % for the event to trigger sometimes, perhaps relate with dl? 
        
        % Noting how some rays don't make it out yet
        plot(u(:,1), u(:,2)); 
        if j==1
            axis([-2 2 -2 2]);
            xlabel('x position'); ylabel('y position');
            x = -1:0.1:1;
            y = ones(length(x),1);
            hold on
            plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);
        end
    end
    pause(0.001);
end
toc;
%% Scattering Relation, Exit Data. 
% Then search for when it 1st leaves the box somehow. 
% Uses squaregaussianrelation.m --> if working correctly, then use to
% get all the data with SGscatteringrelation.m 

% Two velocity interpolations used in squaregaussianrelation.m , 
% Could subtract two methods to see how their accuracy is w.r.t. each other

tic;
ds = 0;  % Or can make it smaller/larger?
uNSEW = zeros((Nedge-1),(Nangle-1),3); 
% The cell of all exit-pos/vel for each given incidence pos/angle. 
% To access: Use for ex. cellfun(@(v) v(4), uNSEW(:,1)) ?
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. So,
% The 1st entry ~ location on the edge. 2nd entry ~ angle of incidence.

for i = 1:Nedge-1
    ds = dl*2*i*(Nedge-i)/Nedge;
    for j = 1:Nangle-1
        % u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   % Only for left (West) edge!
        % u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   % This is for bottom (South) edge!
        % u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   % Only for right (East) edge!
         u0 = [1-i*dl;1;-j*dphi];   % This is for top (North) edge!
        uNSEW(i,j,:) = SGthetaexitrelation(u0,ds);
    end
end
%%
% To check qualitatively with the plot above:
disp(uNSEW(end,:,:))
toc;
% Be careful if we wanted to switch the order of for loops
% It looked consistent with all 4 sides. 
%% Entire scattering relation
% Use: (Takes a while to compute)
tic;
[uW, uN, uE, uS] = SGthetascatteringrelation(Nedge, Nangle);

% To check with prior notes I guess:
disp(uN(end,:,:))
%%
% To get the components, I guess just do a for loop, e.g. for x:
xexit = zeros(Nedge-1,Nangle-1);
i = 1 : Nedge-1;
j = 1 : Nangle-1;
xexit(i,j) = uNSEW(i,j,1);
yexit(i,j) = uNSEW(i,j,2);
disp(xexit(end,:) - uNSEW(end,:,1))
%%
% Check exit:
M = max(abs(xexit),abs(yexit)) - 1;
disp(norm(M(:),Inf))
toc;

%% Scripts to run ode45 for ray tracing
% This is just for the square case and the metric g is a Gaussian. 
% These are just notes / scratch work primarily. 

% function dH = gaussianmetrictheta(s,u)
% % u(1) = x coordinate, u(2) = y coordinate, 
% % u(3) = x velocity, u(4) = y velocity.
% % This has no bearing on the domain anyways!
% 
%     dH = zeros(3,1);
%     % The positions:
%     dH(1) = cos(u(3)).*exp(u(1).^2/2 + u(2).^2/2);
%     dH(2) = sin(u(3)).*exp(u(1).^2/2 + u(2).^2/2);
%     % The momenta: 
%     dH(3) = (-u(1).*sin(u(3))+u(2).*cos(u(3))).*exp(u(1).^2/2 + u(2).^2/2);
%     
% end

%% Checking the ODE evolution is fine and iterating.
tic;
Nedge = 2^5;   % number of partitions of sides.
Nangle = 2^5;   % number of partitions of angle.
dphi = pi/Nangle;   % This needs to be taken better care of
dl = 2/Nedge; 
% Two loops over n,N
% Can switch the order of for loops to see the change in ray traces with 
% entry position or entry angle graphically.
%%
% To terminate when either |x|,|y| >= 1
options = odeset('Events',@sgEventsFcn); 
% Here, the events function was defined as:
% function [value,isterminal,direction] = sgEventsFcn(s,u)
% % Terminate ODE evolution when the position is outside of square.
% value = max(0,1 - max(abs(u(1)), abs(u(2))));
% isterminal = 1;
% direction = -1;
% end
%%
ds=2;
for i = 1:Nedge-1
    clf
    % INCLUDE THIS: So that the event will trigger when the starting point
    % is close to the boundary: 
    % ds = dl*2*i*(Nedge-i)/Nedge;    % leave off to see rays fully exit.
    
    for j = 1:Nangle-1
        % u0 = [-1;1-i*dl;cos(j*dphi-pi/2); sin(j*dphi-pi/2)];   % Only for left edge!
        % u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   % This is for bottom edge!
        % u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   % Only for right edge!
         u0 = [1-i*dl;1;-j*dphi];   % This is for top edge!
        [s,u,~,ue] = ode45(@gaussianmetrictheta, [0,ds], u0, options);
        % The arclength interval of ODE evolution has to be small enough
        % for the event to trigger sometimes, perhaps relate with dl? 
        
        % Noting how some rays don't make it out yet
        plot(u(:,1), u(:,2)); 
        if j==1
            axis([-2 2 -2 2]);
            xlabel('x position'); ylabel('y position');
            x = -1:0.1:1;
            y = ones(length(x),1);
            hold on
            plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);
        end
    end
    pause(0.001);
end
toc;
%% Scattering Relation, Exit Data. 
% Then search for when it 1st leaves the box somehow. 
% Uses squaregaussianrelation.m --> if working correctly, then use to
% get all the data with SGscatteringrelation.m 

% Two velocity interpolations used in squaregaussianrelation.m , 
% Could subtract two methods to see how their accuracy is w.r.t. each other

tic;
ds = 0;  % Or can make it smaller/larger?
uNSEW = zeros((Nedge-1),(Nangle-1),3); 
% The cell of all exit-pos/vel for each given incidence pos/angle. 
% To access: Use for ex. cellfun(@(v) v(4), uNSEW(:,1)) ?
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. So,
% The 1st entry ~ location on the edge. 2nd entry ~ angle of incidence.

for i = 1:Nedge-1
    ds = dl*2*i*(Nedge-i)/Nedge;
    for j = 1:Nangle-1
        % u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   % Only for left (West) edge!
        % u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   % This is for bottom (South) edge!
        % u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   % Only for right (East) edge!
         u0 = [1-i*dl;1;-j*dphi];   % This is for top (North) edge!
        uNSEW(i,j,:) = SGthetaexitrelation(u0,ds);
    end
end
%%
% To check qualitatively with the plot above:
disp(uNSEW(end,:,:))
toc;
% Be careful if we wanted to switch the order of for loops
% It looked consistent with all 4 sides. 
%% Entire scattering relation
% Use: (Takes a while to compute)
tic;
[uW, uN, uE, uS] = SGthetascatteringrelation(Nedge, Nangle);

% To check with prior notes I guess:
disp(uN(end,:,:))
%%
% To get the components, I guess just do a for loop, e.g. for x:
xexit = zeros(Nedge-1,Nangle-1);
i = 1 : Nedge-1;
j = 1 : Nangle-1;
xexit(i,j) = uNSEW(i,j,1);
yexit(i,j) = uNSEW(i,j,2);
disp(xexit(end,:) - uNSEW(end,:,1))
%%
% Check exit:
M = max(abs(xexit),abs(yexit)) - 1;
disp(norm(M(:),Inf))
toc;

