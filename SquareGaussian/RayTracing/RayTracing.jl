function gaussianmetric(s,u)
# u(1) = x coordinate, u(2) = y coordinate,
# u(3) = x velocity, u(4) = y velocity.
# This has no bearing on the domain anyways!

    dH = zeros(4);
    # The positions:
    dH[1] = u[3].*exp(u[1].^2 + u[2].^2);
    dH[2] = u[4].*exp(u[1].^2 + u[2].^2);
    # The momenta:
    dH[3] = -(u[3].^2 + u[4].^2).*u[1].*exp(u[1].^2 + u[2].^2);
    dH[4] = -(u[3].^2 + u[4].^2).*u[2].*exp(u[1].^2 + u[2].^2);
    return dH

end



#####################################################################


# function sgEventsFcn(s,u)
#   # Terminate ODE evolution when the position is outside of square.
#   value = max(0,1 - max(abs(u(1)), abs(u(2))));
#   isterminal = 1;
#   direction = -1;
#   return value,isterminal,direction;
# end


#####################################################################


function squaregaussianrelation(u0,ds)
# This will take an initial condition and evolve the ODE for the scattering
# relation. It will also make sure to get when it exits, and return the
# exit position and velocity.
# Maybe use a fixed step so that way it's more modulated?

# options = odeset('Events',sgEventsFcn);   #### This is where we have problems??

~,u = ode45(gaussianmetric, u0, [0,ds]) #options);
# [~,u] = ode45(@gaussianmetric, [0,1], u0);  % for example
# it should be kept adapative for the interval of length

# uf = u0;   % initialize it, but I don't think we need to?

kx = find(x -> ( abs(x[1]) >= 1 ), u[2:end]);
ky = find(x -> ( abs(x[2]) >= 1 ), u[2:end]);

if isempty(kx) && isempty(ky)
    u0 = u[end];
    uf = squaregaussianrelation(u0,ds);
    return uf;
end

if isempty(kx)
   k = ky[1]+1;
end
if isempty(ky)
   k = kx[1]+1;
end
if !isempty(kx) && !isempty(ky)
   k = min(kx[1],ky[1]) + 1;
end

uf = u[k];
dy = u[k][2]-u[k-1][2];
dx = u[k][1]-u[k-1][1];
m = dy/dx;
# Can we just use slope m (or, dy and dx) for the directions/velocities?

if abs(uf[1]) > 1

    uf[1] = uf[1]/abs(uf[1]);   # So it preserves the sign.

    # Interpolated y component
    uf[2] = u[k-1][2] + m*(uf[1]-u[k-1][1]);

    # Interpolated velocities
    uf[3] = dx;
    uf[4] = dy;

    # To check difference between two interpolation approaches:
    # Use r to interpolate exit velocities with weighted averaging:
    # r = (uf(1)-u(k-1,1))/(u(k,1)-u(k-1,1));
    # uf(3) = uf(3) - (u(k-1,3) + r*(u(k,3) - u(k-1,3)));
    # uf(4) = uf(4) - (u(k-1,4) + r*(u(k,4) - u(k-1,4)));

end

if abs(uf[2]) > 1

    uf[2] = uf[2]/abs(uf[2]);   # So it preserves the sign.

    # Interpolated x component
    uf[1] = u[k-1][1] + (uf[2]-u[k-1][2])/m;

    # Interpolated velocities
    uf[3] = dx;
    uf[4] = dy;

    # To check difference between two interpolation approaches:
    # Use r to interpolate exit velocities with weighted averaging:
    # r = (uf(1)-u(k-1,1))/(u(k,1)-u(k-1,1));
    # uf(3) = uf(3) - (u(k-1,3) + r*(u(k,3) - u(k-1,3)));
    # uf(4) = uf(4) - (u(k-1,4) + r*(u(k,4) - u(k-1,4)));

end

uf[3] = uf[3]/sqrt(uf[3]^2 + uf[4]^2);
uf[4] = uf[4]/sqrt(uf[3]^2 + uf[4]^2);

return uf;

end


###############################################################


function SGscatteringrelation(Nedge, Nangle)
# Gather the exit data for the incidence data along the
# West(left),North(top),East(right),South(bottom) edges.

dl = 2/Nedge;
dphi = pi/Nangle;
uW = cell((Nedge-1),(Nangle-1));
uN = cell((Nedge-1),(Nangle-1));
uS = cell((Nedge-1),(Nangle-1));
uE = cell((Nedge-1),(Nangle-1));
# For the cells, each row is a point on the boundary edge and each collumn is
# an angle of incidence.

ds = 1;

for i = 1:Nedge-1
    for j = 1:Nangle-1
        u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   # Only for left edge!
        uW[i,j] = squaregaussianrelation(u0,ds);
        u0 = [-1 + i*dl; -1; cos(j*dphi); sin(j*dphi)];   # Only for bottom edge!
        uS[i,j] = squaregaussianrelation(u0,ds);
        u0 = [1; -1+i*dl; cos(j*dphi + pi/2); sin(j*dphi + pi/2)];   # Only for right edge!
        uE[i,j]= squaregaussianrelation(u0,ds);
        u0 = [1 - i*dl; 1; cos(-j*dphi); sin(-j*dphi)];   # Only for top edge!
        uN[i,j] = squaregaussianrelation(u0,ds);
    end
end

return uW, uS, uE, uN;

end
