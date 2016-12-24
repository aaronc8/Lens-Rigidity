## Scripts to solve the Lens Rigidity Problem on the square.
@everywhere using ODE
@everywhere using PyPlot
@everywhere using Polynomials
@everywhere using Interpolations
@everywhere include("RayTracing.jl");
@everywhere include("MetricSolver.jl");

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)));
square(x,y) =  - heaviside(x.+1.0).*heaviside(-x.+1.0).*heaviside(-y.+1.0).*heaviside(y.+1.0);
dxsquare(x,y) = -sign(x).*(4.0.*heaviside(x.+1.0).*heaviside(-x.-1.0)) - sign(x).*(4.0.*heaviside(x.-1.0).*heaviside(-x.+1.0))
dysquare(x,y) = -sign(y).*(4.0.*heaviside(y.+1.0).*heaviside(-y.-1.0)) - sign(y).*(4.0.*heaviside(y.-1.0).*heaviside(-y.+1.0))
gradsquare(x,y) = [ dxsquare(x,y) , dysquare(x,y) ];

@everywhere N = 5;
@everywhere Nedge = 2^N;   # number of partitions of sides.
@everywhere Nangle = 2^N;   # number of partitions of angle.
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

# Defining the physical domain in which the inversion will take place.
@everywhere x=-1:2/Nedge:1;
@everywhere y=x;

# Defining the wavespeed
@everywhere cxy=exp(0.5.*( x.^2 .+ y'.^2 ));
# Computing the gradient using finite differences
@everywhere gradcxy,hesscxy=GradHessFinDiff(cxy);


gradcxyExact = zeros(length(x),length(y),2);
gradcxyExact[:,:,1] = cxy.*(x   .+ 0*y');
gradcxyExact[:,:,2] = cxy.*(0*x'.+ y);
# surf(x,y,cxy)

# k = 2:Nedge;
# dcdx = zeros(Nedge+1, Nedge+1);
# dcdx[:,k] = (Nedge/4).*(cxy[:,k+1] - cxy[:,k-1]);
# dcdx[:,1] = (Nedge/4).*(-3.0.*cxy[:,1] + 4.0.*cxy[:,2] - cxy[:,3]);
# dcdx[:,end] = (Nedge/4).*(-3.0.*cxy[:,end] + 4.0.*cxy[:,end-1] - cxy[:,end-2]);




# c2xy = ones(Nedge,Nedge);

# defining the mesh
@everywhere knots=([xi for xi in x ], [yi for yi in y]);
# metric,dmetric = generateMetric(knots,cxy);
# cspd,gradcspd,hesscspd = generateMetric(cxy,gradcxy,hesscxy);

# building interpolation objects
@everywhere cspd,gradcspd,hesscspd=generateMetric(knots,cxy,gradcxy,hesscxy);

x = -1:2/Nedge:1;   # For the actual grid ...
y = x;
cxy = exp(0.5.*(x.^2 .+ y'.^2));
# Can check surf(x,y,cxy) just to be safe....Also
display("Check the approximation agrees with the original metric:")
display(norm(cxy-cspd(x,y),Inf))


xx = -1:0.5/Nedge:1;   # For the actual grid ...
yy = x;
cxy = exp(0.5.*(xx.^2 .+ yy'.^2));
println("Check the approximation error when interpolating")
display(norm(cxy-cspd(xx,yy),Inf))


# Can check surf(x,y,cxy) just to be safe....Also
println("Check the approximation of the derivative")
display(display(maximum(abs(gradcxyExact[:,:,1]-gradcspd(x,y)[1] ))))


# If BSplines, The square is from [1/Nx, 1] x [1/Ny,1] so we need to rescale:
# gpre,dgpre = generateMetric(cxy);
# Nx,Ny = size(cxy);
# Dx = 1-1/Nx; Dy = 1-1/Ny; hx = 1/Nx; hy = 1/Ny;
# metric(x,y) = gpre( Dx*(x+hx+Dx/2)/2, Dx*(y+hy+Dx/2)/2 );

###############################################################################
# Construct the Hamiltonian system using the interpolated metric:
# How to get Julia to output a function?? If dH is done manually, works fine...
@everywhere metric(x::Float64,y::Float64)=cspd(x,y).^2;
@everywhere dmetric(x::Float64,y::Float64)=2*cspd(x,y).*gradcspd(x,y);
@everywhere dHtheta=makeHamiltonian(metric,dmetric,true);  # Or:
@everywhere dH=makeHamiltonian(metric,dmetric,false);
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
uWexact, uSexact, uEexact, uNexact = SGscatteringrelation(false,gaussianmetric,Nedge,Nangle,ds);
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
Xg,sout = generatePath(square,dH, [x0,y0,cos(a0),sin(a0)], ds);
Xgtheta,souttheta = generatePath(square, dHtheta, [x0,y0,a0], ds);
Xexact,sexact = generatePath(square,gaussianmetric, [x0,y0,cos(a0),sin(a0)], ds);
Xtheta,stheta = generatePath(square,gaussianmetrictheta, [x0,y0, a0], ds);
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
xx = -1:0.1:1;
yy = ones(length(xx),1);
plot(xx,yy); plot(xx,-yy); plot(yy,xx); plot(-yy,xx);
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

# initial condition for the Jacobian
J0 = eye(4,4); # J0theta = eye(3,3);
Jacobian = geodesicJacobian(M,Xg,J0,sout);    # It is outputted as a 4x4, not a 16x1.

###############################################################################
## Now to set up the system to solve for the update:
## E.g.
RayData = linearMismatch(Jacobian,Xg,sout,cspd,gradcspd,Nedge)
## It will remove entries below a certain tolerance; like <= 1e-5 entries become zeros.
# goodindices = find(x -> x >= 1e-6, RayData)

# We had the exact data above already generated,
bEexact = zeros(4, (Nedge-1)*(Nangle-1));
sEout = zeros((Nedge-1)*(Nangle-1));
bNexact = zeros(4, (Nedge-1)*(Nangle-1));
sNout = zeros((Nedge-1)*(Nangle-1));
bWexact = zeros(4, (Nedge-1)*(Nangle-1));
sWout = zeros((Nedge-1)*(Nangle-1));
bSexact = zeros(4, (Nedge-1)*(Nangle-1));
sSout = zeros((Nedge-1)*(Nangle-1));

for k = 1 : (Nedge-1)*(Nangle-1)
  bEexact[:,k] = uEexact[k][1:4]; sEout[k] = uEexact[k][5];
  bNexact[:,k] = uNexact[k][1:4]; sNout[k] = uNexact[k][5];
  bWexact[:,k] = uWexact[k][1:4]; sWout[k] = uWexact[k][5];
  bSexact[:,k] = uSexact[k][1:4]; sSout[k] = uSexact[k][5];
end

bexact = append!(bEexact[:],bNexact[:]); bexact = append!(bexact,bWexact[:]); bexact = append!(bexact,bSexact[:]);
sout = append!(sEout,sNout); sout = append!(sout,sWout); sout = append!(sout, sSout);

A,b = MismatchSystem(square,cspd,gradcspd,hesscspd,Nedge,Nangle,bexact,sout);

# Just for fun:
update = A\b;
update = reshape(update,Nedge+1,Nedge+1);

figure(3)
x = -1:2/Nedge:1;
y = x;
cxy = exp(0.5.*( x.^2 .+ y'.^2 ));
surf(x,y,cxy)

figure(4)
cxy = update + cxy;
surf(x,y,cxy);


# ########################  Old mistake: #######################################
## Check that it's reasonable with a concrete example:
# JacobianTheta = geodesicJacobian(Mtheta,Xgtheta,J0theta,souttheta);
# lambda(x,y) = metric(x,y) - exp(x.^2 .+ y.^2);   # Just to see dV is working / computing something.
# dlambda(x,y) = dmetric(x,y) - [2.0.*x.*exp(x.^2 .+ y.^2), 2.0.*y.*exp(x.^2 .+ y.^2)];
#
# dVdg = dVfrechet(cspd,gradcspd,lambda,dlambda);
# display(dVdg([0,1,0.25,1]))   #Just to see it's working
#
# ## And applying it to find a sample mismatch:
# Kni = linearMismatch(Jacobian,dVdg,Xg,sout);
# display(Kni)
# ## Put it all together to get all the Kni:
# ## Since doing with just dg = identity, the result should be very close still
# ## when the geodesics are short, i.e. in the corners of the matrix. Near Endpoints
# ## and only when theta is near 0 or pi.
#
# tic()
# KW,KS,KE,KN = geodesicMismatch(square, cspd,gradcspd,hesscspd,lambda,dlambda,Nedge,Nangle,ds);
# toc()
