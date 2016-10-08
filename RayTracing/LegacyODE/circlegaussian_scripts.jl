using ODE
using PyPlot
include("RayTracing.jl")

Nrotate = 2^4;
Nangle = 2^5;
dphi = pi/Nangle;
dtheta = 2*pi/Nrotate;

### PART 1: Just for checking the ODE is setup and evolving correctly

for i = 1:Nrotate-1
    clf()
    for j = 1:Nangle-1
        u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)];
        # ~,u = ODE.ode45(gaussianmetric, u0, [0.0,1.0]);
        ~,u = ODE.ode45(gaussianmetric, u0, [0.0,1.0], stopevent = (s,u) -> ( u[1]^2 + u[2]^2 > 1 ) );
        # Noting how some rays don't make it out yet
        u1 = map(z->z[1],u);
        u2 = map(z->z[2],u);
        plot(u1,u2);
        ax = gca();
        if j == 1
          ax[:set_xlim]([-2,2]);
          ax[:set_ylim]([-2,2]);
          ax[:set_xlabel]("x-position");
          ax[:set_ylabel]("y-position");
        theta = 0:0.1:2*pi+0.1;
        x = cos(theta);
        y = sin(theta);
        plot(x,y);
        end
        pause(0.0000001);
    end
end

### PART 2:  This is for just extracting the exit data for each entrance data

ds = 1;  # Or can make it smaller/larger?
uExit = cell((Nrotate-1),(Nangle-1));
for i = 1:Nrotate-1
    for j = 1:Nangle-1
        u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)];
        uExit[i,j] = circlegaussianrelation(u0,ds);
    end
end

# Check with the last graphic:
print(uExit[end,:])

# To access the individual components:
# xexit = map(x -> x[1], uExit);  # Or also: xexit = [x[1] for x in uExit];
# yexit = map(x -> x[2], uExit);
# vxexit = map(x -> x[3], uExit);
# vyexit = map(x -> x[4], uExit);

### Part 3: Entire scattering relation
# Use the following:
uExitTot = CGscatteringrelation(Nrotate,Nangle,ds);

# Check with the prior one (need to pick NSEW)
print(uExit[end,:] - uExitTot[end,:])
