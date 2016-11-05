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

# Construct the Hamiltonian system using the interpolated metric:
# How to get Julia to output a function?? If dH is done manually, works fine...
dHtheta = makeHamiltonianTheta(metric,dmetric);  # Or:
dH = makeHamiltonian(metric,dmetric);
## And checked its ray evolution is okay now

## Warning: May have to use the same form the whole way through! That way,
## dH will get properly overwritten every iteration?
## This is definitely necessary if we define dH manually and not thru the fn.

## Check working with scattering relation:
tic();
uWt, uSt, uEt, uNt = SGscatteringrelation(1,dHtheta,Nedge,Nangle,ds);
toc();  # Holy crap I suck this is so slow compared to the original ... due to eval-ing interpolations??

tic();
uW, uS, uE, uN = SGscatteringrelation(0,dH,Nedge,Nangle,ds);
toc();  # and even slower .....

tic();
xtexit = map(x -> x[1], uWt);
ytexit = map(x -> x[2], uWt);
display("Check every ray exited, max(xexit,yexit) in absolute value:")
display(max(abs(xtexit),abs(ytexit)))
display("Check the max across all entries:")
display(norm(max(abs(xtexit),abs(ytexit))[:] - 1,Inf))

xexit = map(x -> x[1], uW);
yexit = map(x -> x[2], uW);
display("Check every ray exited, max(xexit,yexit) in absolute value:")
display(max(abs(xexit),abs(yexit)))
display("Check the max across all entries:")
display(norm(max(abs(xexit),abs(yexit))[:] - 1,Inf))

display("Check the max across all entries of differences:")
display(norm(max(abs(xexit - xtexit),abs(yexit-ytexit))[:], Inf))
toc();
## The diference is more noticable because we are using a linear approx to g?

## Next part: compute the path interpolants: e.g.
