## Scripts to solve the Lens Rigidity Problem on the square.
using ODE
using PyPlot
using Polynomials
using Interpolations
include("RayTracing.jl");
include("MetricSolver.jl");

Nedge = 2^5;   # number of partitions of sides.
Nangle = 2^5;   # number of partitions of angle.
dphi = pi/Nangle;   # This needs to be taken better care of
dl = 2/Nedge;
ds = 2;

###############################################################################
## Testing Development of the Integration functions for K_i^n
###############################################################################
## Generating the metric based on grid points. For now, just discretize the exact one.
## Then, all the things should look / behave close to the original things we've done.

x = linspace(-1,1,Nedge);
y = x';
c2xy = exp(x.^2 .+ y.^2);
# c2xy = ones(Nedge,Nedge);

knots = ([x for x = linspace(-1,1,Nedge)], [y for y = linspace(-1,1,Nedge)]);
metric,dmetric = generateMetric(knots,c2xy);
# Can check surf(x,y,cxy) just to be safe....Also
display(norm(c2xy-metric(x,y),Inf))

# If BSplines, The square is from [1/Nx, 1] x [1/Ny,1] so we need to rescale:
# gpre,dgpre = generateMetric(cxy);
# Nx,Ny = size(cxy);
# Dx = 1-1/Nx; Dy = 1-1/Ny; hx = 1/Nx; hy = 1/Ny;
# metric(x,y) = gpre( Dx*(x+hx+Dx/2)/2, Dx*(y+hy+Dx/2)/2 );

###############################################################################
# Construct the Hamiltonian system using the interpolated metric:
# How to get Julia to output a function?? If dH is done manually, works fine...
dHtheta = makeHamiltonian(metric,dmetric,true);  # Or:
dH = makeHamiltonian(metric,dmetric,false);
## And checked its ray evolution is okay now

## Warning: May have to use the same form the whole way through! That way,
## dH will get properly overwritten every iteration?
## This is definitely necessary if we define dH manually and not thru the fn.

## Check working with scattering relation:
tic();
uWt, uSt, uEt, uNt = SGscatteringrelation(true,dHtheta,Nedge,Nangle,ds);
toc();  # Holy crap I suck this is so slow compared to the original ... due to eval-ing interpolations??

tic();
uW, uS, uE, uN = SGscatteringrelation(false,dH,Nedge,Nangle,ds);
toc();  # and even slower .....

tic();
uWexact, uSexact, uEexact, uNexact = SGscatteringrelation(true,gaussianmetrictheta,Nedge,Nangle,ds);
toc()

tic();
xtexit = map(x -> x[1], uWt);
ytexit = map(x -> x[2], uWt);
display("Check every ray exited, max(xexit,yexit) in absolute value:")
display(max(abs(xtexit),abs(ytexit)))
display("Check the max across all entries:")
display(norm(max(abs(xtexit),abs(ytexit))[:] , Inf))

xexit = map(x -> x[1], uW);
yexit = map(x -> x[2], uW);
display("Check every ray exited, max(xexit,yexit) in absolute value:")
display(max(abs(xexit),abs(yexit)))
display("Check the max across all entries:")
display(norm(max(abs(xexit),abs(yexit))[:] , Inf))

display("Check the max across all entries of differences:")
display(norm(max(abs(xexit - xtexit),abs(yexit-ytexit))[:], Inf))
toc();
## The diference is more noticable because we are using a linear approx to g?

xexact = map(x -> x[1], uWexact);
yexact = map(x -> x[2], uWexact);
display("Check the max across all entries with exact: (Off due to approx?)")
display(norm(max(abs(xexact - xtexit),abs(yexact-ytexit))[:], Inf))
display(norm(max(abs(xexit - xexact),abs(yexit-yexact))[:], Inf))

###############################################################################
## Next part: compute the path interpolants: e.g. X_g(t,X0).
# One example:
Xg,sout = generatePath(dH, [-1,0,cos(pi/4),sin(pi/4)], ds);
Xgtheta,souttheta = generatePath(dHtheta, [-1,0,pi/4], ds);
Xexact,sexact = generatePath(gaussianmetric, [-1,0,cos(pi/4),sin(pi/4)], ds);
Xtheta,stheta = generatePath(gaussianmetrictheta, [-1, 0, pi/4], ds);
display(" Check the path deviance, Without and With Theta perspective:" )
t1range = 0:0.01:sexact; t2range = 0:0.01:stheta;
display(norm(Xg(t1range) - Xexact(t1range), Inf))
display(norm(Xgtheta(t2range) - Xtheta(t2range), Inf))
## Also plotting the path was pretty clear :)
## This lets us use it for the integral for K_n^i now.

###############################################################################
## Now, we have to compue the Hessian of H (called M), to be used for the Jacobian.
# But, how to check correctness...? Check with original? ._.  :|
M = HamiltonianHess(metric,dmetric,false);
Mtheta = HamiltonianHess(metric,dmetric,true);
J0 = eye(4,4); J0theta = eye(3,3);
Jacobian = geodesicJacobian(M,Xg,J0,sout);
JacobianTheta = geodesicJacobian(Mtheta,Xgtheta,J0theta,souttheta);
