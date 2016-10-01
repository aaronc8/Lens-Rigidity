function dH = gaussianmetrictheta(s,u)
% u(1) = x coordinate, u(2) = y coordinate, 
% u(3) = x velocity, u(4) = y velocity.
% This has no bearing on the domain anyways!

    dH = zeros(3,1);
    % The positions:
    dH(1) = cos(u(3)).*exp(u(1).^2/2 + u(2).^2/2);
    dH(2) = sin(u(3)).*exp(u(1).^2/2 + u(2).^2/2);
    % The momenta: 
    dH(3) = (u(1).*sin(u(3)) - u(2).*cos(u(3))).*exp(u(1).^2/2 + u(2).^2/2);
    
end