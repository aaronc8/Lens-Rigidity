function [value,isterminal,direction] = cgEventsFcn(s,u)
% Terminate ODE evolution when the position is outside of square.
value = max(0,1 - sqrt(u(1)^2 + u(2)^2));
isterminal = 1;
direction = -1;   
end

