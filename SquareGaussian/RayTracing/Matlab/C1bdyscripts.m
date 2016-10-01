%% Scripts to run ode45 for ray tracing on more general cases.
% This is just for the square case and the metric g is a Gaussian. 
% These are just notes / scratch work primarily. 
% The boundary has to be C1 so we only can really just test the circle. 
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
        [s,u,~,ue] = ode45(@gaussianmetric, [0,1], u0, options);
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
toc;
%%
% To check qualitatively with the plot above:
disp(uExit(end,:,:))