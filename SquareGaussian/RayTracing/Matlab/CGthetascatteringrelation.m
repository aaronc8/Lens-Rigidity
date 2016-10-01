function [uTotExit] = CGthetascatteringrelation(Nrotate, Nangle)
% Gather the exit data for the incidence data along the 
% West(left),North(top),East(right),South(bot) edges. 

dtheta = 2*pi/Nrotate;
dphi = pi/Nangle;
uTotExit = zeros((Nrotate-1),(Nangle-1),3);
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. 

for i = 1:Nrotate-1
    for j = 1:Nangle-1
        ds = 4*j*dphi*(Nangle-j)/Nangle;
        u0 = [cos(i*dtheta), sin(i*dtheta), i*dtheta + pi/2 + j*dphi]; 
        uTotExit(i,j,:) = CGthetaexitrelation(u0,ds);
    end
end

end