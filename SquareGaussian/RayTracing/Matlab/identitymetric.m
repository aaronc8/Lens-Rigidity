function dH = identitymetric(s,u)
% This is just to make sure that all lines are straight in this case. 
% u(1) = x coordinate, u(2) = y coordinate, 
% u(3) = x velocity, u(4) = y velocity.

    dH = zeros(4,1);
    % The positions:
    dH(1) = u(3);
    dH(2) = u(4);
    % The momenta: 
    dH(3) = 0;
    dH(4) = 0;
    
end
