function dH = gaussianmetric(s,u)
% u(1) = x coordinate, u(2) = y coordinate, 
% u(3) = x velocity, u(4) = y velocity.
% This has no bearing on the domain anyways!

    dH = zeros(4,1);
    % The positions:
    dH(1) = u(3).*exp(u(1).^2 + u(2).^2);
    dH(2) = u(4).*exp(u(1).^2 + u(2).^2);
    % The momenta: 
    dH(3) = -(u(3).^2 + u(4).^2).*u(1).*exp(u(1).^2 + u(2).^2);
    dH(4) = -(u(3).^2 + u(4).^2).*u(2).*exp(u(1).^2 + u(2).^2);
    
end
