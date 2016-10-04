using ODE
using PyPlot
using Polynomials
include("C1RayTracing.jl")

Nrotate = 2^5;
Nangle = 2^5;
dphi = pi/Nangle;
dtheta = 2*pi/Nrotate;
ds=2.0;

circle(x,y) = x.^2 + y.^2 - 1;
circlegrad(x,y) = [2.*x, 2.*y];

## Everything has toggled on/off using gaussianmetric or gaussianmetrictheta! ##

### PART 1: Just for checking the ODE is setup and evolving correctly

for i = 1:Nrotate-1
    pause(1/Nrotate)
    clf()
    for j = 1:Nangle-1
        # ds = j*dphi*(Nangle-j)/Nangle;  # May not want to, so that we can see the whole evolution.
        # u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)];
        # ~,u = ode45(gaussianmetric, u0, [0.0,ds]);
        u0 = [cos(i*dtheta), sin(i*dtheta), i*dtheta + pi/2 + j*dphi];
        ~,u = ode45(gaussianmetrictheta, u0, [0.0,ds]);
        k = find( x -> circle(x[1],x[2]) > 0, u[2:end]);
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

## Everything has toggled on/off using gaussianmetric or gaussianmetrictheta! ##

### PART 2:  This is for just extracting the exit data for each entrance data

ds = 1;  # Or can make it smaller/larger?
uExit = Array{Array}((Nrotate-1),(Nangle-1));
tic()
for i = 1:Nrotate-1
    for j = 1:Nangle-1
        ds = 2*j*dphi*(Nangle-j)/Nangle;
        # u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)];
        # uExit[i,j] = scatteringRelation(gaussianmetric, circle, circlegrad, u0,ds);
        u0 = [cos(i*dtheta), sin(i*dtheta), i*dtheta + pi/2 + j*dphi];
        uExit[i,j] = scatteringRelation(gaussianmetrictheta,circle,circlegrad,u0,ds);
    end
end

# Check with the last graphic:
print(uExit[end,:])
xexit = map(x -> x[1], uExit);
yexit = map(x -> x[2], uExit);
display("Check that they lie on the boundary of domain:")
display(circle(xexit,yexit))
display("Max across all entries:")
display(norm(circle(xexit,yexit)[:],Inf))
toc()


############################################################################
## Try with an ellipse?

Nrotate = 2^5;
Nangle = 2^5;
dphi = pi/Nangle;
dtheta = 2*pi/Nrotate;

a = 2.0; b =0.5; ds = 2/min(a,b);
ellipse(x,y) = a^2.0.*x.^2. + b^2.0.*y.^2. - 1.;
ellipsegrad(x,y) = [2*a^2.*x, 2*b^2.*y];

## Everything has toggled on/off using gaussianmetric or gaussianmetrictheta! ##

### PART 1: Just for checking the ODE is setup and evolving correctly

for i = 1:Nrotate-1
    pause(.01)
    clf()
    for j = 1:Nangle-1
        # ds = j*dphi*(Nangle-j)/Nangle;  # May not want to, so that we can see the whole evolution.

        x = cos(i*dtheta)./a; y = sin(i*dtheta)./b;
        v = ellipsegrad(x,y);
        # Make sure we get the right angle of the normal:
        angle = atan(v[2]/v[1]);
        if sign(v[1]) > 0
          angle = angle+pi;
        end
        u0 = [x,y, cos(angle + pi/2 - j*dphi), sin(angle + pi/2 - j*dphi)];
        ~,u = ode45(gaussianmetric, u0, [0.0,ds]);
        # u0 = [x,y, angle + pi/2 - j*dphi];
        # ~,u = ode45(gaussianmetrictheta, u0, [0.0,ds]);
        k = find( x -> ellipse(x[1],x[2]) > 0, u[2:end]);
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
          ax[:set_xlim]([-ds/1.5,ds/1.5]);
          ax[:set_ylim]([-ds/1.5,ds/1.5]);
          ax[:set_xlabel]("x-position");
          ax[:set_ylabel]("y-position");
          theta = 0:0.1:2*pi+0.1;
          x = cos(theta)./a;
          y = sin(theta)./b;
          plot(x,y);
        end
    end
end

## Everything has toggled on/off using gaussianmetric or gaussianmetrictheta! ##

### PART 2:  This is for just extracting the exit data for each entrance data

ds = 1;  # Or can make it smaller/larger?
uExit = Array{Array}((Nrotate-1),(Nangle-1));
tic()
for i = 1:Nrotate-1
    for j = 1:Nangle-1
        ds = 2*j*dphi*(Nangle-j)/Nangle;
        x = cos(i*dtheta)./a; y = sin(i*dtheta)./b;
        # Make sure we get the correct normal angle:
        v = ellipsegrad(x,y);
        angle = atan(v[2]/v[1]);
        if sign(v[1]) > 0
          angle = angle+pi;
        end
        # u0 = [x,y, cos(angle + pi/2 - j*dphi), sin(angle + pi/2 - j*dphi)];
        # uExit[i,j] = scatteringRelation(gaussianmetric, ellipse, ellipsegrad, u0,ds);
        u0 = [x,y, angle + pi/2 - j*dphi];
        uExit[i,j] = scatteringRelation(gaussianmetrictheta,ellipse,ellipsegrad,u0,ds);
    end
end

# Check with the last graphic:
print(uExit[end,:])
xexit = map(x -> x[1], uExit);
yexit = map(x -> x[2], uExit);
display("Check that they lie on the boundary of domain:")
display(ellipse(xexit,yexit))
display("Max across all entries:")
display(norm(ellipse(xexit,yexit)[:],Inf))
toc()


#######################################################################
## Something more random - here's a graph from Stewart's Book:
# r = 2+cos(3Θ) -> cos(3A) = cos(2A)cos(A) - sin(2A)sin(A) = cos^3(A) - cos(A)sin^2(A) - 2sin^2(A)cos(A).
# cos(3A) = cos(A) (cos(A) - sin(A))^2.
# r = 2 + x/r * (x-y)^2/r^2 <=> r^4 = r^3 + x(x-y)^2.
## Or something simpler, r = 2 + cos(2t)? Maybe simpler? lol...
