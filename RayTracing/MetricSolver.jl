## Functions used to solve for the Riemannian Metric.

function generateMetric(knots,cxy)
# Takes the computed values of cxy at (x,y) ∈ knots and interpolates.


itp = interpolate(knots, cxy, Gridded(Linear()));
metric(x,y) = itp[x,y];
gradmetric(x,y) = gradient(itp,x,y);  # More efficient to use "gradient!" ?
# itp = interpolate(cxy,BSpline(Quadratic(Flat())), OnGrid())
# Nx,Ny = size(cxy);
# # Renormalize to a unit interval in both dimensions, too:
# metric(x,y) = itp[x*Nx, y*Ny];
# gradmetric(x,y) = gradient(itp, x*Nx, y*Ny);

# Check the domain, but just return:
return metric,gradmetric;

end

###############################

function makeHamiltonianTheta(metric::Function,dmetric::Function)
# Takes any metric function (at least isotropic for now) and generates its
# corresponding Hamiltonian flow, dH.
  function dH(s,u)
    V = zeros(3);
    V[1] = cos(u[3]).*metric(u[1],u[2])^0.5;
    V[2] = sin(u[3]).*metric(u[1],u[2])^0.5;
    V[3] = 0.5.*dot( [sin(u[3]), -cos(u[3])] , dmetric(u[1],u[2]) )./metric(u[1],u[2]).^0.5 ;
    return V;
  end
  return dH;
end

## Because if I tried to put them into one method, both dH functions got created...

function makeHamiltonian(metric::Function,dmetric::Function)
  function dH(s,u)
    V = zeros(4);
    V[1] = u[3].*metric(u[1],u[2]);
    V[2] = u[4].*metric(u[1],u[2]);
    V[3] = -0.5.*(u[3].^2 + u[4].^2).*dmetric(u[1],u[2])[1];
    V[4] = -0.5.*(u[3].^2 + u[4].^2).*dmetric(u[1],u[2])[2];
    return V;
  end
  return dH;
end

# if angleTF == true
#   dH1(s,u) = cos(u[3]).*metric(u[1],u[2])^0.5;
#   dH2(s,u) = sin(u[3]).*metric(u[1],u[2])^0.5;
#   dH3(s,u) = 0.5.*dot( [sin(u[3]), -cos(u[3])] , dmetric(u[1],u[2]) )./metric(u[1],u[2]).^0.5 ;
#   dH(s,u) = [dH1(s,u), dH2(s,u), dH3(s,u)];
#   return dH;
# end
#
# ## If we don't use angle
# if angleTF == false
#   dH1(s,u) = u[3].*metric(u[1],u[2]);
#   dH2(s,u) = u[4].*metric(u[1],u[2]);
#   dH3(s,u) = -0.5.*(u[3].^2 + u[4].^2).*dmetric(u[1],u[2])[1];
#   dH4(s,u) = -0.5.*(u[3].^2 + u[4].^2).*dmetric(u[1],u[2])[2];
#   dH(s,u) = [dH1(s,u) dH2(s,u) dH3(s,u) dH4(s,u)];
# end

###############################

function generatePath(metric::Function,u0,ds)
# Generate the paths with any metric.
# Interpolate the Ray evolution so that we have the path as a parametric function.
# Similar to method "squaregaussianrelation" but we have to keep the whole path.
# We still need to solve for exit time s_out for upper bound of integral....

s,u = ode45(metric, u0, [0.0,ds]);
kx = find(x -> ( abs(x[1]) >= 1 ), u[2:end]);
ky = find(x -> ( abs(x[2]) >= 1 ), u[2:end]);

while isempty(kx) && isempty(ky)
    u0 = u[end];
    sadd,uadd = ode45(metric,u0,[0.0,ds]);
    kx = find(x -> ( abs(x[1]) >= 1 ), uadd[2:end]);
    ky = find(x -> ( abs(x[2]) >= 1 ), uadd[2:end]);
    sadd = sadd + s[end];
    shift = length(u) - length(uadd);
    s = append!(s,sadd);
    u = append!(u,uadd);
end

kshift = length(u) - shift;

if isempty(kx)
   k = ky[1]+1+kshift;
end
if isempty(ky)
   k = kx[1]+1+kshift;
end
if !isempty(kx) && !isempty(ky)
   k = min(kx[1],ky[1]) + 1 + kshift;
end

## Construct interpolant for exit time:
u1 = u[k-1]; u2 = u[k];
s1 = s[k-1][1]; s2 = s[k][1]; spread = s2-s1;

# u1 = u[end-1]; u2 = uf;
# s1 = s[end-1][1]; s2 = s[end][1]; sf = s2-s1;

# Need more data points (specifically 5) for 4th order:
dt = spread/4.0;
tspan = 0:dt:spread;
# sp,up = ODE.ode45(gaussianmetric,u1,tspan);
s,u = ode45(metric,u1,tspan);    # is still being adaptive???
up = u[1:5]; sp = s[1:5] + s1;
up[end] = u2; sp[end] = s2;
# Interpolate positions:
xp = map( v -> v[1], up);
xI = polyfit(sp,xp);
yp = map( v-> v[2], up);
yI = polyfit(sp, yp);

uf = u2; # Need to use for corner cases, if both |x|,|y| > 1.

# Because it's not root finding the exit point:
if abs(uf[1]) >= 1  # x-exit
  dom = xI - sign(u2[1]);
  dxI = polyder(xI);
  sout = newtonbisection(dom, dxI, s1, s2, 1e-5);
  uf[1] = xI(sout); uf[2] = yI(sout);
end


if abs(uf[2]) >= 1  # Also to account for if y exits before x.
  dom = yI - sign(uf[2]);
  dyI = polyder(yI);
  sout = newtonbisection(dom, dyI, s1, s2, 1e-5);
  uf[1] = xI(sout); uf[2] = yI(sout);
end

# Interpolate velocity:
for k = 1:length(uf)-2
  vp = map( t -> t[k+2], up);
  v = polyfit(sp,vp);
  uf[k+2] = v(sout);
end

return uf;

end

###############################

function HamiltonianHess(H::Function)
# Form the Hessian (it is almost but not quite the Hessian) of H to
# Compute/Solve for the Jacobian of the geodesic. (This is M)


end

################################

function geodesicJacobian(M,J0)
# Solve for the Jacobian matrix of geodesic w.r.t. initial condition.
# J0 is usually identity, but it could be something else (e.g. B if broken)


end

#################################

function linearMismatch()
# Compute the integral of the mismatch of geodesic evolution.

end

##################################
