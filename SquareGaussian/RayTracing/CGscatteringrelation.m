function [uTotExit] = CGscatteringrelation(Nrotate, Nangle)
% Gather the exit data for the incidence data along the 
% West(left),North(top),East(right),South(bot) edges. 

dtheta = 2*pi/Nrotate;
dphi = pi/Nangle;
uTotExit = cell((Nrotate-1),(Nangle-1));
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. 

ds = 1;

for i = 1:Nrotate-1
    for j = 1:Nangle-1
        u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)]; 
        uTotExit{i,j} = circlegaussianrelation(u0,ds);
    end
end

end