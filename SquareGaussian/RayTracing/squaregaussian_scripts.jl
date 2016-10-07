using ODE
using PyPlot
include("RayTracing.jl");

Nedge = 2^5;   # number of partitions of sides.
Nangle = 2^5;   # number of partitions of angle.
dphi = pi/Nangle;   # This needs to be taken better care of
dl = 2/Nedge;
ds = 1.25;
# options = odeset("Events",sgEventsFcn);


### Trying Events in ODE.jl:
# u0 = [-1,0,1/2^0.5,1/2^0.5];
# s0 = 0.0
# ~,u = ODE.ode45(gaussianmetric,u0,[0.,2.],stopevent = (s,u) -> (abs(u[1]) > 1 || abs(u[2]) > 1));
# u1 = map(z->z[1],u);
# u2 = map(z->z[2],u);
# plot(u1,u2);
# ax = gca();
#   ax[:set_xlim]([-2,2]);
#   ax[:set_ylim]([-2,2]);
#   ax[:set_xlabel]("x-position");
#   ax[:set_ylabel]("y-position");
# x = -1:0.1:1;
# y = ones(length(x),1);
# plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);
# # Event seems to not trigger root finding :(
#
# # Another option is to try using
# # isoutofdomain(y) = abs(y[1]) > 1 || abs(y[2]) > 1;  # adding isnan(y) makes it a vector...
#
# # Another approach is to define each approach to solving:
# ode=ODE.ExplicitODE(s0,u0,gaussianmetric)  # This is a difference - it needs gaussianmetric to be defined differently.
# stepper=ODE.RKStepperAdaptive{:rk45,Float64}()
# opts=ODE.Options{Float64}(tstop=2.0)
#
# for (s,u) in ODE.solve(ode,stepper,opts)    # Currently, not recognizing RHS of gaussianmetric?
#   if abs(u[1]) > 1 || abs(u[2]) > 1
#     print((s,u))
#     break
#   end
# end
#
# ## This is not working yet ... I tried with gaussianmetric to take 3 entries like s,x,v but didn't seem to work ##
#
# u1 = map(z->z[1],u);
# u2 = map(z->z[2],u);
# plot(u1,u2);
# ax = gca();
#   ax[:set_xlim]([-2,2]);
#   ax[:set_ylim]([-2,2]);
#   ax[:set_xlabel]("x-position");
#   ax[:set_ylabel]("y-position");
# x = -1:0.1:1;
# y = ones(length(x),1);
# plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);


### PART 1: Just for checking the ODE is setup and evolving correctly

for i = 1:Nedge-1
    pause(1/Nedge)
    clf()
    for j = 1:Nangle-1
        # u0 = [-1;1-i*dl;cos(j*dphi-pi/2); sin(j*dphi-pi/2)];   # Only for left edge!
        # u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   # This is for bottom edge!
        # u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   # Only for right edge!
        u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   # This is for top edge!

        ~,u = ode45(gaussianmetric, u0, [0.0,ds]);  # For official ODE.jl
        # ~,u = ODE.ode45(gaussianmetric, u0, [0.0,1.0], stopevent = (s,u) -> (abs(u[1]) > 1 || abs(u[2]) > 1));
        kx = find(x -> ( abs(x[1]) >= 1 ), u[2:end]);
        ky = find(x -> ( abs(x[2]) >= 1 ), u[2:end]);
        if isempty(kx) && !isempty(ky)
           k = ky[1]+1;
        end
        if isempty(ky) && !isempty(kx)
           k = kx[1]+1;
        end
        if !isempty(kx) && !isempty(ky)
           k = min(kx[1],ky[1]) + 1;
        end
        if !isempty(kx) || !isempty(ky)
           u = u[1:k]
        end

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
    end
end

### PART 2:  This is for just extracting the exit data for each entrance data

ds = 1;  # Or can make it smaller/larger?
uNSEW = Array{Array}((Nedge-1),(Nangle-1));
for i = 1:Nedge-1
    for j = 1:Nangle-1
        # u0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   # Only for left (West) edge!
        # u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   # This is for bottom (South) edge!
        # u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   # Only for right (East) edge!
        u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   # This is for top (North) edge!
        uNSEW[i,j] = squaregaussianrelation(u0,ds);
    end
end

# Check with the last graphic:
print(uNSEW[end,:])

# To access the individual components:
# xexit = map(x -> x[1], uNSEW);  # Or also xexit = [x[1] for x in uNSEW];
# yexit = map(x -> x[2], uNSEW);
# vxexit = map(x -> x[3], uNSEW);
# vyexit = map(x -> x[4], uNSEW);

### Part 3: Entire scattering relation
# Use the following:
tic();
uW, uS, uE, uN = SGscatteringrelation(Nedge,Nangle,ds);

# Check with the prior one (need to pick NSEW)
print(uN[end,:] - uNSEW[end,:])

xexit = map(x -> x[1], uNSEW);
yexit = map(x -> x[2], uNSEW);
display("x exits: ")
display(xexit)
display("y exits: ")
display(yexit)
display("Check every ray exited, max(xexit,yexit) in absolute value:")
display(max(abs(xexit),abs(yexit)))
display("Check the max across all entries:")
display(norm(max(abs(xexit),abs(yexit))[:] - 1,Inf))
toc();

#### Now ... how to imitate ODE's events in Matlab in Julia ... use norm??
# Like, the norm is dist(point, boundary) and if it is within a small tolerance
# it is considered at / out of the boundary.
