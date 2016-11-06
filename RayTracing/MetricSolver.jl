## Functions used to solve for the Riemannian Metric.

function generateMetric(knots,c2xy)
# Takes the computed values of cxy at (x,y) ∈ knots and interpolates.

itp = interpolate(knots, c2xy, Gridded(Linear()));
metric(x,y) = itp[x,y];
gradmetric(x,y) = gradient(itp,x,y);  # More efficient to use "gradient!" ?

## BSplines Approach: (Needs rescaling.....)
# itp = interpolate(cxy,BSpline(Quadratic(Flat())), OnGrid())
# Nx,Ny = size(cxy);
# # Renormalize to a unit interval in both dimensions, too:
# metric(x,y) = itp[x*Nx, y*Ny];
# gradmetric(x,y) = gradient(itp, x*Nx, y*Ny);

# Check the domain, but just return:
return metric,gradmetric;

end

###############################

# function makeHamiltonianTheta(metric::Function,dmetric::Function)
# # Takes any metric function (at least isotropic for now) and generates its
# # corresponding Hamiltonian flow, dH.
#   function dH(s,u)
#     V = zeros(3);
#     V[1] = cos(u[3]).*metric(u[1],u[2])^0.5;
#     V[2] = sin(u[3]).*metric(u[1],u[2])^0.5;
#     V[3] = 0.5.*dot( [sin(u[3]), -cos(u[3])] , dmetric(u[1],u[2]) )./metric(u[1],u[2]).^0.5 ;
#     return V;
#   end
#   return dH;
# end

## Because if I tried to put them into one method, both dH functions got created...
## Fix: Name them two different things?

function makeHamiltonian(metric::Function,dmetric::Function, theta)

  function dH(s,u)
    V = zeros(4);
    V[1] = u[3].*metric(u[1],u[2]);
    V[2] = u[4].*metric(u[1],u[2]);
    V[3] = -0.5.*(u[3].^2 + u[4].^2).*dmetric(u[1],u[2])[1];
    V[4] = -0.5.*(u[3].^2 + u[4].^2).*dmetric(u[1],u[2])[2];
    return V;
  end

  function dHtheta(s,u)
    Vtheta = zeros(3);
    Vtheta[1] = cos(u[3]).*metric(u[1],u[2])^0.5;
    Vtheta[2] = sin(u[3]).*metric(u[1],u[2])^0.5;
    Vtheta[3] = 0.5.*dot( [sin(u[3]), -cos(u[3])] , dmetric(u[1],u[2]) )./metric(u[1],u[2]).^0.5 ;
    return Vtheta;
  end

  if theta == true
    return dHtheta;
  end
  if theta == false
    return dH;
  end

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
kshift = 0;

while isempty(kx) && isempty(ky)
    u0 = u[end];
    sadd,uadd = ode45(metric,u0,[0.0,ds]);
    kx = find(x -> ( abs(x[1]) >= 1 ), uadd[2:end]);
    ky = find(x -> ( abs(x[2]) >= 1 ), uadd[2:end]);
    sadd = sadd + s[end];
    shift = length(u) - length(uadd);
    kshift = length(u) - shift;
    s = append!(s,sadd);
    u = append!(u,uadd);
end

if isempty(kx)
   k = ky[1]+1+kshift;
end
if isempty(ky)
   k = kx[1]+1+kshift;
end
if !isempty(kx) && !isempty(ky)
   k = min(kx[1],ky[1]) + 1 + kshift;
end

# Now to give the interpolated path.
knots = (s[1:k],);
xval = map( t -> t[1], u[1:k]);
yval = map( t -> t[2], u[1:k]);
xitp = interpolate(knots, xval, Gridded(Linear()));
yitp = interpolate(knots, yval, Gridded(Linear()));
# Xg(s) = zeros(length(u0));
if length(u0) == 3   ## This is being overwritten by the other if statement still,,?
  thetaval = map(t -> t[3], u[1:k]);
  thetaitp = interpolate(knots,thetaval,Gridded(Linear()));
  Xgtheta(s) = [xitp[s], yitp[s], thetaitp[s]];
end

if length(u0) == 4
  vxval = map(t -> t[3], u[1:k]);
  vyval = map(t -> t[4], u[1:k]);
  vxitp = interpolate(knots, vxval, Gridded(Linear()));
  vyitp = interpolate(knots, vyval, Gridded(Linear()));
  Xg(s) = [xitp[s], yitp[s], vxitp[s], vyitp[s]];
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
s,u = ode45(metric,u1,tspan);    # tspan only ensures at least 5 pts
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
  # uf[1] = xI(sout); uf[2] = yI(sout);  # Don't need exit data exactly...
end

if abs(uf[2]) >= 1  # Also to account for if y exits before x.
  dom = yI - sign(uf[2]);
  dyI = polyder(yI);
  sout = newtonbisection(dom, dyI, s1, s2, 1e-5);
  # uf[1] = xI(sout); uf[2] = yI(sout); # Don't need exit data
end


if length(u0) == 3
  return Xgtheta,sout;
end

if length(u0) == 4
  return Xg,sout;
end

end

###############################

function HamiltonianHessForcing(g::Function,dg::Function, theta)
# Form the Hessian (it is almost but not quite the Hessian) of H to
# Use to compute/solve for the Jacobian of the geodesic. (This is M)

## BE Careful about the ordering!! Perhaps, go the collumnwise convention?

if theta == false
  function M(u)
    M = zeros(4,4);
    ## Build so that if we did reshape(M,4,4) it's the Hessian.

    # First 4 are H_x,v
    grad = dg(u[1],u[2]);
    M[1] = u[3].*grad[1];
    M[2] = u[4].*grad[1];
    M[5] = u[3].*grad[2];
    M[6] = u[4].*grad[2];
    # Next 4 are -H_x,x
    M[3] = 0;
    M[4] = -0.5.*(grad[1]+grad[2]).*(u[3].^2 + u[4].^2);
    M[7] = M[4];
    M[8] = 0;
    # Next 4 are H_v,v
    M[9] = g(u[1],u[2]);
    M[10] = 0;
    M[13] = 0;
    M[14] = M[9];
    # Last 4 are -H_v,x = -H_x,v transpose!
    M[11] = -M[1];
    M[12] = -M[5];
    M[15] = -M[2];
    M[16] = -M[6];
    return M;
  end
  #
  # function F(s,J)
  #   # Assume that J is still a matrix!
  #   return M(s,J)*J;
  # end

  return M;
end

if theta == true
  function Mtheta(u)
    Mtheta = zeros(3,3);
    c = g(u[1],u[2])^0.5;
    grad = dg(u[1],u[2]);
    xvel = cos(u[3]);
    yvel = sin(u[3]);
    # H_theta,x:
    Mtheta[1] = 0.5.*xvel.*grad[1];
    Mtheta[2] = 0.5.*yvel.*grad[1];
    Mtheta[4] = 0.5.*xvel.*grad[2];
    Mtheta[5] = 0.5.*yvel.*grad[2];
    # H_x,x
    g21 = grad[1]+grad[2];
    Mtheta[3] = -(xvel.*g21 + yvel./(2*c.^2).*grad[1].^2 - xvel.*grad[1].*grad[2]./(2*c.^2))./(2*c);
    Mtheta[6] = -(yvel.*g21 - xvel./(2*c.^2).*grad[2].^2 + yvel.*grad[1].*grad[2]./(2*c.^2))./(2*c);
    # H_theta,theta
    Mtheta[7] = -yvel.*c;
    Mtheta[8] = xvel.*c;
    # H_x,theta
    Mtheta[9] = -(xvel*grad[1] + yvel*grad[2])./(2*c);
    return Mtheta;
  end

  # function Ftheta(s,J)
  #   # Assume that J is still a matrix!
  #   return Mtheta(s,j)*J;
  # end
  #
  # return M;
end
#
# return M;
#
end


################################

function geodesicJacobian(M::Function,X::Function,J0,sout)
# Solve for the Jacobian matrix of geodesic w.r.t. initial condition.
# J0 is usually identity, but it could be something else (e.g. B if broken)
# Since we have the Group Property, I don't think we need to specify initial s?
function F(s,J)
  n = length(J0[1,:]);
  J = reshape(J,n,n);
  return (M(X(s))*J)[:];
end

s,J = ode45(F,J0[:],[0.0,sout]);
knots = (s,);

function Jacobian(t)
  Jac = zeros(size(J0));
  for k = 1:length(J0)
    Jcomp = map(a -> a[k], J);
    itpk = interpolate(knots,Jcomp,Gridded(Linear()));
    Jac[k] = itpk[t];
  end

  return Jac;
end

return Jacobian;

end

#################################

function linearMismatch()
# Compute the integral of the mismatch of geodesic evolution.

end

##################################
