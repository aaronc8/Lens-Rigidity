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


x = -2:1/Nedge:2;
y = x;
cxy = zeros(length(x),length(x));
cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1] = exp(0.5.*x[Nedge+1 : 3*Nedge+1]'.^2 .+ 0.5.*y[Nedge+1 : 3*Nedge+1].^2);
# For a circle - see which points are inside and then assign a value if it is inside.
gradcxy,hesscxy = GradHess(cxy,-2,2);
cxy = cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1]



# c2xy = ones(Nedge,Nedge);

# knots = ([x for x = -1:2/Nedge:1], [y for y = -1:2/Nedge:1]);
# metric,dmetric = generateMetric(knots,cxy);
cspd,gradcspd,hesscspd = generateMetric(cxy,gradcxy,hesscxy);

x = -1:2/Nedge:1;   # For the actual grid ...
y = x;
cxy = exp(0.5.*(x'.^2 .+ y.^2));
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

## Warning: May have to use the same form the whole way through! That way,
## dH will get properly overwritten every iteration?
## This is definitely necessary if we define dH manually and not thru the fn.

## Check working with scattering relation:
tic();
uWt, uSt, uEt, uNt = SGscatteringrelation(true,dHtheta,Nedge,Nangle,ds);
toc();  # Holy crap much slower  ... due to eval-ing interpolations??

tic();
uW, uS, uE, uN = SGscatteringrelation(false,dH,Nedge,Nangle,ds);
toc();  # and even slower ..... Also, some bugs if least squares fails (matrix gets singular)

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
# M = HamiltonianHess(metric,dmetric,false);
# Mtheta = HamiltonianHess(metric,dmetric,true);
M = HamiltonianHess(cspd,gradcspd,hesscspd);
J0 = eye(4,4); # J0theta = eye(3,3);
Jacobian = geodesicJacobian(M,Xg,J0,sout);
JacobianTheta = geodesicJacobian(Mtheta,Xgtheta,J0theta,souttheta);
