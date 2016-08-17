function [uW, uN, uE, uS] = SGscatteringrelation(Nedge, Nangle)
% Gather the exit data for the incidence data along the 
% West(left),North(top),East(right),South(bot) edges. 

dl = 2/Nedge;
dphi = pi/Nangle;
uW = cell((Nedge-2),(Nangle-2));
uN = cell((Nedge-2),(Nangle-2));
uS = cell((Nedge-2),(Nangle-2));
uE = cell((Nedge-2),(Nangle-2));
% For the cells, each row is a point on the boundary and each collumn is
% the angle of incidence. 

ds = 1;
% Only for the West (Left) Edge:
for i = 1:Nedge-2
    for j = 1:Nangle-2
        u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   % Only for left edge!
        uW{i,j} = squaregaussianrelation(u0,ds);
    end
end

% Now for the South (Bottom) Edge:
for i = 1:Nedge-2
    for j = 1:Nangle-2
        u0 = [-1 + i*dl; -1; cos(j*dphi); sin(j*dphi)];   % Only for bottom edge!
        uS{i,j} = squaregaussianrelation(u0,ds);
    end
end

% Now the East (Right) Edge:
for i = 1:Nedge-2
    for j = 1:Nangle-2
        u0 = [1; -1+i*dl; cos(j*dphi + pi/2); sin(j*dphi + pi/2)];   % Only for left edge!
        uE{i,j} = squaregaussianrelation(u0,ds);
    end
end

% The North (Top) Edge:
for i = 1:Nedge-2
    for j = 1:Nangle-2
        u0 = [1 - i*dl; 1; cos(-j*dphi); sin(-j*dphi)];   % Only for bottom edge!
        uN{i,j} = squaregaussianrelation(u0,ds);
    end
end


end