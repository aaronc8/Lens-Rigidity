function [value,isterminal,direction] = polarEventsFcn(s,u)
% Terminate ODE evolution when the position is outside of square.
x = u(1); y = u(2);
value = max(0,(x.^2+y.^2).^2 - 2*(x.^2+y.^2).^(3/2) + 3*x.*(x.^2+y.^2) - 4*x.^3);
isterminal = 1;
direction = 0;  % This time direction needs to be 0 or 1... 
end