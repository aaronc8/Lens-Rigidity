function [uW, uN, uE, uS] = SGthetascatteringrelation(Nedge, Nangle)
% Gather the exit data for the incidence data along the 
% West(left),North(top),East(right),South(bot) edges. 

dl = 2/Nedge;
dphi = pi/Nangle;
uW = zeros((Nedge-1),(Nangle-1),3);
uN = zeros((Nedge-1),(Nangle-1),3);
uS = zeros((Nedge-1),(Nangle-1),3);
uE = zeros((Nedge-1),(Nangle-1),3);
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. 

for i = 1:Nedge-1
    ds = dl*2*i*(Nedge-i)/Nedge;
    for j = 1:Nangle-1
        u0 = [-1; 1-i*dl; j*dphi - pi/2];   % Only for left edge!
        uW(i,j,:) = SGthetaexitrelation(u0,ds);
        u0 = [-1 + i*dl; -1; j*dphi];   % Only for bottom edge!
        uS(i,j,:) = SGthetaexitrelation(u0,ds);
        u0 = [1; -1+i*dl; j*dphi + pi/2];   % Only for right edge!
        uE(i,j,:) = SGthetaexitrelation(u0,ds);
        u0 = [1 - i*dl; 1; -j*dphi];   % Only for top edge!
        uN(i,j,:) = SGthetaexitrelation(u0,ds);
    end
end

end

% Prior stuff just for checking:

% % Now for the South (Bottom) Edge:
% for i = 1:Nedge-2
%     for j = 1:Nangle-2
%         u0 = [-1 + i*dl; -1; cos(j*dphi); sin(j*dphi)];   % Only for bottom edge!
%         uS{i,j} = squaregaussianrelation(u0,ds);
%     end
% end
% 
% % Now the East (Right) Edge:
% for i = 1:Nedge-2
%     for j = 1:Nangle-2
%         u0 = [1; -1+i*dl; cos(j*dphi + pi/2); sin(j*dphi + pi/2)];   % Only for right edge!
%         uE{i,j} = squaregaussianrelation(u0,ds);
%     end
% end
% 
% % The North (Top) Edge:
% for i = 1:Nedge-2
%     for j = 1:Nangle-2
%         u0 = [1 - i*dl; 1; cos(-j*dphi); sin(-j*dphi)];   % Only for top edge!
%         uN{i,j} = squaregaussianrelation(u0,ds);
%     end
% end
