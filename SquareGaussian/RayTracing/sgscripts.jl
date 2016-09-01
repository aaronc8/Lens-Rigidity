#using ODE
#using PyPlot
#include("gaussianmetric.jl");

Nedge = 2^4;   # number of partitions of sides.
Nangle = 2^5;   # number of partitions of angle.
dphi = pi/Nangle;   # This needs to be taken better care of
dl = 2/Nedge;
# options = odeset("Events",sgEventsFcn);

### PART 1: Just for checking the ODE is setup and evolving correctly

for i = 1:Nedge-2
    clf()
    for j = 1:Nangle-2
        # u0 = [-1;1-i*dl;cos(j*dphi-pi/2); sin(j*dphi-pi/2)];   # Only for left edge!
        # u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   # This is for bottom edge!
         u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   # Only for right edge!
        # u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   # This is for top edge!
        ~,u = ode45(gaussianmetric, u0, [0,1]);
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
        x = -1:0.1:1;
        y = ones(length(x),1);
        plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);
        end
        pause(0.000001);
    end
end

### PART 2:  This is for just extracting the exit data for each entrance data

ds = 1;  # Or can make it smaller/larger?
uNSEW = cell((Nedge-2),(Nangle-2));
for i = 1:Nedge-2
    for j = 1:Nangle-2
        # u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   # Only for left (West) edge!
        # u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   # This is for bottom (South) edge!
         u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   # Only for right (East) edge!
        # u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   # This is for top (North) edge!
        uNSEW[i,j] = squaregaussianrelation(u0,ds);
    end
end

# Check with the last graphic:
print(uNSEW[end,:])

# To access the individual components:
# xexit = map(x -> x[1], uNSEW);
# yexit = map(x -> x[2], uNSEW);
# vxexit = map(x -> x[3], uNSEW);
# vyexit = map(x -> x[4], uNSEW);

### Part 3: Entire scattering relation
# Use the following:
uW, uS, uE, uN = SGscatteringrelation(Nedge,Nangle);

# Check with the prior one (need to pick NSEW)
print(uE[end,:] - uNSEW[end,:])

#### Now ... how to imitate ODE's events in Matlab in Julia ... use norm??
# Like, the norm is dist(point, boundary) and if it is within a small tolerance
# it is considered at / out of the boundary.
