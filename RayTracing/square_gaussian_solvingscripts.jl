## Scripts to solve the Lens Rigidity Problem on the square.
using ODE
using PyPlot
using Polynomials
using Interpolations
include("RayTracing.jl");
include("MetricSolver.jl");

N = 5;
Nedge = 2^N;   # number of partitions of sides.
Nangle = 2^N;   # number of partitions of angle.
dphi = pi/Nangle;   # This needs to be taken better care of
dl = 2/Nedge;
ds = 2;

###############################################################################
## Testing Development of the Integration functions for K_i^n
###############################################################################
## Generating the metric based on grid points. For now, just discretize the exact one.
## Then, all the things should look / behave close to the original things we've done.


# x = -2:1/Nedge:2;
# y = x;
# cxy = zeros(length(x),length(x));
# cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1] = exp(0.5.*x[Nedge+1 : 3*Nedge+1]'.^2 .+ 0.5.*y[Nedge+1 : 3*Nedge+1].^2);
# # For a circle - see which points are inside and then assign a value if it is inside.
# gradcxy,hesscxy = GradHessFFT(cxy,-2,2);
# cxy = cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1]

x = -1:2/Nedge:1;
y = x;
cxy = exp(0.5.*( x.^2 .+ y'.^2 ));
gradcxy,hesscxy = GradHessFinDiff(cxy);

# k = 2:Nedge;
# dcdx = zeros(Nedge+1, Nedge+1);
# dcdx[:,k] = (Nedge/4).*(cxy[:,k+1] - cxy[:,k-1]);
# dcdx[:,1] = (Nedge/4).*(-3.0.*cxy[:,1] + 4.0.*cxy[:,2] - cxy[:,3]);
# dcdx[:,end] = (Nedge/4).*(-3.0.*cxy[:,end] + 4.0.*cxy[:,end-1] - cxy[:,end-2]);




# c2xy = ones(Nedge,Nedge);

knots = ([x for x = -1:2/Nedge:1], [y for y = -1:2/Nedge:1]);
# metric,dmetric = generateMetric(knots,cxy);
# cspd,gradcspd,hesscspd = generateMetric(cxy,gradcxy,hesscxy);
cspd,gradcspd,hesscspd = generateMetric(knots,cxy,gradcxy,hesscxy);

x = -1:2/Nedge:1;   # For the actual grid ...
y = x;
cxy = exp(0.5.*(x.^2 .+ y'.^2));
# Can check surf(x,y,cxy) just to be safe....Also
display("Check the approximation agrees with the original metric:")
display(norm(cxy-cspd(x,y),Inf))

# If BSplines, The square is from [1/Nx, 1] x [1/Ny,1] so we need to rescale:
# gpre,dgpre = generateMetric(cxy);
# Nx,Ny = size(cxy);
# Dx = 1-1/Nx; Dy = 1-1/Ny; hx = 1/Nx; hy = 1/Ny;
# metric(x,y) = gpre( Dx*(x+hx+Dx/2)/2, Dx*(y+hy+Dx/2)/2 );

###############################################################################
# Construct the Hamiltonian system using the interpolated metric:
# How to get Julia to output a function?? If dH is done manually, works fine...
metric(x,y) = cspd(x,y).^2;
dmetric(x,y) = 2*cspd(x,y).*gradcspd(x,y);
dHtheta = makeHamiltonian(metric,dmetric,true);  # Or:
dH = makeHamiltonian(metric,dmetric,false);
## And checked its ray evolution is okay now

figure(1)
for i = 1:Nedge-1
    pause(1/Nedge)
    clf()
    for j = 1:Nangle-1
        # u0 = [-1;1-i*dl;cos(j*dphi-pi/2); sin(j*dphi-pi/2)];   # Only for left edge!
        # u0 = [-1;1-i*dl;j*dphi-pi/2];   # 3-vars
        # u0 = [-1+i*dl;-1;cos(j*dphi); sin(j*dphi)];   # This is for bottom edge!
        # u0 = [-1+i*dl;-1;j*dphi];   # 3-vars
         u0 = [1;-1+i*dl;cos(j*dphi+pi/2); sin(j*dphi+pi/2)];   # Only for right edge!
        # u0 = [1;-1+i*dl;j*dphi+pi/2];   # 3-vars
        # u0 = [1-i*dl;1;cos(-j*dphi); sin(-j*dphi)];   # This is for top edge!
        # u0 = [1-i*dl;1;-j*dphi];   # 3-vars

        # ~,u = ode45(gaussianmetric, u0, [0.0,ds]);  # For official ODE.jl
        ~,u = ode45(dH,u0,[0.0,ds]);
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

## Warning: May have to use the same form the whole way through! That way,
## dH will get properly overwritten every iteration?
## This is definitely necessary if we define dH manually and not thru the fn.

## Check working with scattering relation:
tic();
uWt, uSt, uEt, uNt = SGscatteringrelation(true,dHtheta,Nedge,Nangle,ds);
toc();  # Holy crap much slower  ... due to eval-ing interpolations??

tic();
uW, uS, uE, uN = SGscatteringrelation(false,dH,Nedge,Nangle,ds);
toc();  # and even slower ...  Also, some bugs if least squares fails (matrix gets singular)
# It is not as slow if the interpolant is better maybe because the adaptive stepping goes better.

tic();
uWexact, uSexact, uEexact, uNexact = SGscatteringrelation(true,gaussianmetrictheta,Nedge,Nangle,ds);
toc()

tic();
xtexit = map(x -> x[1], uWt);
ytexit = map(x -> x[2], uWt);
display("Check every ray exited, max(xexit,yexit) in absolute value, 3 Var:")
display(max(abs(xtexit),abs(ytexit)))
display("Check the max across all entries:")
display(norm(max(abs(xtexit),abs(ytexit))[:] , Inf))

xexit = map(x -> x[1], uW);
yexit = map(x -> x[2], uW);
display("Check every ray exited, max(xexit,yexit) in absolute value, 4 Var:")
display(max(abs(xexit),abs(yexit)))
display("Check the max across all entries:")
display(norm(max(abs(xexit),abs(yexit))[:] , Inf))

display("Check the max across all entries of differences between 3 var and 4 var:")
display(norm(max(abs(xexit - xtexit),abs(yexit-ytexit))[:], Inf))
toc();
## The diference is more noticable because we are using a linear approx to g?

xexact = map(x -> x[1], uWexact);
yexact = map(x -> x[2], uWexact);
display("Check the max across all entries with exact: (Off due to approx?)")
display(norm(max(abs(xexact - xtexit),abs(yexact-ytexit))[:], Inf))
display(norm(max(abs(xexit - xexact),abs(yexit-yexact))[:], Inf))

## The 3-Var is MUCH faster - maybe it is more prudent to evolve in 3 vars, then
## just take cos and sin of the theta component to get 4 vars.

###############################################################################
## Next part: compute the path interpolants: e.g. X_g(t) given X0.
## Just output s_exit because we're interpolating all points from 0 to exit.
# One example:
x0 = 1.0; y0 = -0.5; a0 = -5*pi/6;
Xg,sout = generatePath(dH, [x0,y0,cos(a0),sin(a0)], ds);
Xgtheta,souttheta = generatePath(dHtheta, [x0,y0,a0], ds);
Xexact,sexact = generatePath(gaussianmetric, [x0,y0,cos(a0),sin(a0)], ds);
Xtheta,stheta = generatePath(gaussianmetrictheta, [x0,y0, a0], ds);
display(" Check the path deviance, Without and With Theta perspective:" )
t1range = 0:0.01:sexact; t2range = 0:0.01:stheta;
display(norm(Xg(t1range) - Xexact(t1range), Inf))
display(norm(Xgtheta(t2range) - Xtheta(t2range), Inf))
## Also plotting the path was pretty clear :)
figure(2)
clf()
fourvar = Xg(t1range); threevar = Xgtheta(t2range);
fourexact = Xexact(t1range); threeexact = Xtheta(t2range);
ax = gca();
ax[:set_xlim]([-2,2]);
ax[:set_ylim]([-2,2]);
ax[:set_xlabel]("x-position");
ax[:set_ylabel]("y-position");
x = -1:0.1:1;
y = ones(length(x),1);
plot(x,y); plot(x,-y); plot(y,x); plot(-y,x);
plot(fourvar[1,:][1][:] - fourexact[1,:][1][:], fourvar[2,:][1][:] - fourexact[2,:][1][:], marker = "+"); # is a point lol.
plot(fourvar[1,:][1][:], fourvar[2,:][1][:], marker = "o", color = "r");
plot(threevar[1,:][1][:], threevar[2,:][1][:], marker = "x", color = "g");


## This lets us use it for the integral for K_n^i now.

###############################################################################
## Now, we have to compue the Hessian of H (called M), to be used for the Jacobian.
# But, how to check correctness...? Check with original? (._.) ...  :[  ... (.-.)
# M = HamiltonianHess(metric,dmetric,false);
# Mtheta = HamiltonianHess(metric,dmetric,true);
M = HamiltonianHess(cspd,gradcspd,hesscspd);
J0 = eye(4,4); # J0theta = eye(3,3);
Jacobian = geodesicJacobian(M,Xg,J0,sout);    # It is outputted as a 4x4, not a 16x1.
# JacobianTheta = geodesicJacobian(Mtheta,Xgtheta,J0theta,souttheta);
lambda(x,y) = 1;   # Just to see dV is working / computing something.
dlambda(x,y) = zeros(2,1);

dVdg = dVfrechet(cspd,gradcspd,lambda,dlambda);
display(dVdg([0,1,0.25,1]))   #Just to see it's working

## And applying it to find the mismatch:
Kni = linearMismatch(Jacobian,dVdg,Xg,sout);
display(Kni)


###################################################
# Put it all together to get all the Kni:

KW,KS,KE,KN = geodesicMismatch(cspd,gradcspd,hesscspd,lambda,dlambda,Nedge,Nangle,ds);
## Since doing with just dg = identity, the result should be very close still
## when the geodesics are short, i.e. in the corners of the matrix. Near Endpoints
## and only when theta is near 0 or pi. 
