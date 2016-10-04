using ODE
using PyPlot
using Polynomials
include("RayTracing.jl")

Nrotate = 2^5;
Nangle = 2^5;
dphi = pi/Nangle;
dtheta = 2*pi/Nrotate;
ds=2.0;

### PART 1: Just for checking the ODE is setup and evolving correctly

for i = 1:Nrotate-1
    pause(1/Nrotate)
    clf()
    for j = 1:Nangle-1
        # ds = j*dphi*(Nangle-j)/Nangle;  # May not want to, so that we can see the whole evolution.
        u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)];
        ~,u = ode45(gaussianmetric, u0, [0.0,ds]);
        k = find( x -> x[1]^2 + x[2]^2 > 1, u[2:end]);
        if !isempty(k)
            u = u[1:k[1]+1];
        end

        # ~,u = ODE.ode45(gaussianmetric, u0, [0.0,ds], stopevent = (s,u) -> ( u[1]^2 + u[2]^2 > 1 ) );
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
    end
end

### PART 2:  This is for just extracting the exit data for each entrance data

ds = 1;  # Or can make it smaller/larger?
uExit = Array{Array}((Nrotate-1),(Nangle-1));
for i = 1:Nrotate-1
    for j = 1:Nangle-1
        ds = 2*j*dphi*(Nangle-j)/Nangle;
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
tic()
uExitTot = CGscatteringrelation(Nrotate,Nangle);

# Check with the prior one (need to pick NSEW)
print(uExit[end,:] - uExitTot[end,:])

xexit = map(x -> x[1], uExit);
yexit = map(x -> x[2], uExit);
display("Check that they lie on the boundary: x^2 + y^2 = 1?")
display(xexit.^2 + yexit.^2)
display("Check max across all entries:")
display(norm(xexit[:].^2+yexit[:].^2,Inf))
toc()
