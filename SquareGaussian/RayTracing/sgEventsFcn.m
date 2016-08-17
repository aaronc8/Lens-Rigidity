function [value,isterminal,direction] = sgEventsFcn(s,u)
% Terminate ODE evolution when the position is outside of square.
value = max(0,1 - max(abs(u(1)), abs(u(2))));
isterminal = 1;
direction = -1;
end

