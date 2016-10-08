%% Scripts to run ode45 for ray tracing on more general cases.
% This is just for the square case and the metric g is a Gaussian. 
% These are just notes / scratch work primarily. 
% The boundary has to be C1 though, in particular from polar curves... 
%%
% Now, the metric and the options have to be inputted. 
% E.g. we already have CGeventsfcn and also the metric with 4 and 3 vars.
Nrotate = 2^5;   % number of partitions of rotation angle 0 to 2*pi.
Nangle = 2^5;   % number of partitions of angle.
dphi = pi/Nangle;   % This needs to be taken better care of
dtheta = 2*pi/Nrotate; 
options = odeset('Events',@cgEventsFcn); 
ds=2;
tic;
for i = 1:Nrotate-1
    clf
    for j = 1:Nangle-1
        % Include this to ensure events works near the boundary:
        % ds = 4*j*dphi*(Nangle-j)/Nangle;   % Leave off to see rays exit
        u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)]; 
        [s,u,~,ue] = ode45(@gaussianmetric, [0,ds], u0, options);
        % u0 = [cos(i*dtheta), sin(i*dtheta), i*dtheta + pi/2 + j*dphi]; 
        % [s,u,~,ue] = ode45(@gaussianmetrictheta, [0,ds], u0, options);
        % The arclength interval of ODE evolution has to be small enough
        % for the event to trigger sometimes, perhaps relate with dtheta? 
        
        % Noting how some rays don't make it out yet
        plot(u(:,1), u(:,2)); 
        if j==1
            axis([-2 2 -2 2]);
            xlabel('x position'); ylabel('y position');
            theta = 0:0.1:2*pi+0.1;
            x = cos(theta); y = sin(theta);
            hold on
            plot(x,y);
        end
    end
    pause(0.001);
end
toc;

%%
% Now mainly to see it takes either metric:
tic;
ds = 1;  % Or can make it smaller/larger?

% The cell of all exit-pos/vel for each given incidence pos/angle. 
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. So,
% The 1st entry ~ location on the edge. 2nd entry ~ angle of incidence.

% uExit = zeros((Nrotate-1),(Nangle-1),3); 
uExit = zeros(Nrotate-1,Nangle-1,4);
for i = 1:Nrotate-1
    for j = 1:Nangle-1
        ds = 4*j*dphi*(Nangle-j)/Nangle;
        u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)]; 
        uExit(i,j,:) = C1bdyExitRelation(u0,ds,@gaussianmetric,options);
        %u0 = [cos(i*dtheta), sin(i*dtheta), (i*dtheta + pi/2 + j*dphi)];
        %uExit(i,j,:) = C1bdyExitRelation(u0,ds,@gaussianmetrictheta,options);
    end
end

%%
% To check qualitatively with the plot above:
disp(uExit(end,:,:))
%%
% To get the components, I guess just do a for loop, e.g. for x:
xexit = uExit(:,:,1); yexit = uExit(:,:,2);
disp(xexit(end,:) - uExit(end,:,1))
%%
% Check boundary: 
M = xexit.^2 + yexit.^2 - 1;
disp(norm(M(:),Inf))
toc;


%% From a less "nice" Polar Curve: 
% Consider the equation $r = 2+cos(3t)$ ?
G = @(x,y) (x.^2+y.^2).^2 - 2*(x.^2+y.^2).^(3/2) + 3*x.*(x.^2+y.^2) - 4*x.^3;
dG = @(x,y) [ 4*x.*(x.^2+y.^2) - 6*x.*(x.^2+y.^2)^0.5 - 3*x.^2 + 3*y.^2, ...
            4*y.*(x.^2+y.^2) - 6*y.*((x.^2+y.^2))^0.5 + 6*x.*y ];
xpara = @(t) (2+cos(3.*t)).*cos(t);
ypara = @(t) (2+cos(3.*t)).*sin(t);

polaroptions = odeset('Events',@polarEventsFcn); 

Nrotate = 2^5;
Nangle = 2^5;
dphi = pi/Nangle;
dtheta = 2*pi/Nrotate;
ds = 5.0;

%%
% Just the evolution again:

for i = 1:Nrotate-1
    pause(.01)
    clf()
    for j = 1:Nangle-1
        % ds = j*dphi*(Nangle-j)/Nangle;  % May not want to, so that we can see the whole evolution.

        x = xpara(i*dtheta); y = ypara(i*dtheta);
        v = dG(x,y);
        % Make sure we get the right angle of the normal:
        angle = atan(v(2)/v(1));
        if sign(v(1)) > 0
          angle = angle+pi;
        end
%         u0 = [x,y, cos(angle + pi/2 - j*dphi), sin(angle + pi/2 - j*dphi)];
%         [~,u] = ode45(@gaussianmetric, [0.0,ds], u0);
        u0 = [x,y, angle + pi/2 - j*dphi];
        [~,u] = ode45(@gaussianmetrictheta, [0,ds], u0, polaroptions);

        % Noting how some rays don't make it out yet
        plot(u(:,1), u(:,2));
        if j == 1
          axis([-2 4 -3 3]);
          xlabel('x position'); ylabel('y position');
          theta = 0:0.1:2*pi+0.1;
          x = xpara(theta);
          y = ypara(theta);
          plot(x,y);
          hold on
        end
    end
end

%% 
% Now for the exit data gathering: 

tic;
uExit = zeros(Nrotate-1,Nangle-1,3); 
ds = 3.0;

for i = 1:Nrotate-1
    for j = 1:Nangle-1
        % ds = 8*j*dphi*(Nangle-j)/Nangle;
        x = xpara(i*dtheta); y = ypara(i*dtheta);
        v = dG(x,y);
        % Make sure we get the right angle of the normal:
        angle = atan(v(2)/v(1));
        if sign(v(1)) > 0
          angle = angle+pi;
        end
%         u0 = [x,y, cos(angle + pi/2 - j*dphi), sin(angle + pi/2 - j*dphi)];
%         uExit(i,j,:) = C1bdyExitRelation(u0,ds,@gaussianmetric,polaroptions);
        u0 = [x,y, angle + pi/2 - j*dphi];
        uExit(i,j,:) = C1bdyExitRelation(u0,ds,@gaussianmetrictheta,polaroptions);
    end
end

disp(uExit(end,:,:))

%%
% To get the components, I guess just do a for loop, e.g. for x:
xexit = uExit(:,:,1); yexit = uExit(:,:,2);
disp(xexit(end,:) - uExit(end,:,1))
%%
% Check boundary: 
M = G(xexit,yexit) - 1;
disp(norm(M(:),Inf))
toc;



