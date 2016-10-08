# Since the Hamiltonian is zero on these trajectories:
# If we instead use just 3 variables where the velocity is all encapsulated in
# just the direction, theta = u[3].

function gaussianmetrictheta(s,u)
  # u(1) = x coordinate, u(2) = y coordinate,
  # u(3) = x velocity, u(4) = y velocity.
  # This has no bearing on the domain anyways!
      dH = zeros(3);
      # The positions:
      dH[1] = cos(u[3]).*exp(u[1].^2/2 + u[2].^2/2);
      dH[2] = sin(u[3]).*exp(u[1].^2/2 + u[2].^2/2);
      # The momenta:
      dH[3] = (u[1].*sin(u[3]) - u[2].*cos(u[3])).*exp(u[1].^2/2 + u[2].^2/2);
      return dH
end


function gaussianmetric(s,u)
# u(1) = x coordinate, u(2) = y coordinate,
# u(3) = x velocity, u(4) = y velocity.
# This has no bearing on the domain anyways!

    dH = zeros(4);
    # The positions:
    dH[1] = u[3].*exp(u[1].^2 + u[2].^2);
    dH[2] = u[4].*exp(u[1].^2 + u[2].^2);
    # The momenta:
    dH[3] = -(u[3].^2 + u[4].^2).*u[1].*exp(u[1].^2 + u[2].^2);
    dH[4] = -(u[3].^2 + u[4].^2).*u[2].*exp(u[1].^2 + u[2].^2);
    return dH
end

###############################################################################

function scatteringRelation(metric::Function,domain::Function,domaingrad::Function,u0,ds)
  # This will take an initial condition and evolve the ODE for the scattering
  # relation. It will also make sure to get when it exits, and return the
  # exit position and velocity.
  # The domain should be a function of x,y with two arguments.

# s,u = ODE.ode45(gaussianmetric, u0, [0.0,ds], stopevent = (s,u) -> (domain(u) >= 1) );
# uf = u[end];

s,u = ode45(metric, u0, [0.0,ds]);  # Non-legacy ODE.
k = find(x -> domain(x[1],x[2]) > 0, u[2:end]);  # So we don't get the starting pt.

scaling = max(ds,1/ds);
maxiter = 1;
while isempty(k) && maxiter < 1000*scaling   # In case of trapped rays?
# while uf[1]^2 + uf[2]^2 < 1
  maxiter = maxiter + 1;
  uf = u[end];   # We can just input u[end] into ode45?
  s,u = ode45(metric, uf, [0.0,ds]);
  k = find(x -> domain(x[1],x[2]) > 0, u[2:end]);
  #s,u = ODE.ode45(gaussianmetric,uf,[0.0,ds],stopevent = (s,u) -> ( u[1]^2 + u[2]^2 > 1) );
  #uf = u[end];
end
if maxiter >= 1000*scaling   # Suspect trapped rays.
  for k = 1:length(uf)
    uf[k]=Inf;
  end
  return uf;
end

k = k[1] + 1; # because it is index-skewed in the find!
u0 = u[k-1]; uf = u[k];   # Use u0,uf to save storage?
s1 = s[k-1][1]; s2 = s[k][1]; spread = s2-s1;

# Need more data points (specifically 5) for 4th order:
dt = spread/4.0;
tspan = 0:dt:spread;
s,u = ode45(metric,u0,tspan);    # is still being adaptive??
up = u[1:5]; sp = s[1:5] + s1;
up[end] = uf; sp[end] = s2;
xp = map( t -> t[1], up); yp = map( t-> t[2], up);
xI = polyfit(sp,xp); yI = polyfit(sp, yp);
dxI = polyder(xI); dyI = polyder(yI);

sdomain(s) = domain(xI(s),yI(s));
sdomaingrad(s) = dot(domaingrad(xI(s),yI(s)) , [dxI(s), dyI(s)] );
sout = newtonbisection(sdomain, sdomaingrad, s1, s2, 1e-5);

uf[1] = xI(sout); uf[2] = yI(sout);

# Interpolate velocities:
for k = 1:length(uf)-2
  vp = map( t -> t[k+2], up);
  v = polyfit(sp,vp);
  uf[k+2] = v(sout);
end

return uf;

end

# function scatteringRelation(metric::function,domain::function,domaingrad::function,u0,ds)
#   # The same function, but with LegacyODE commands....
#
# s,u = ODE.ode45(gaussianmetric, u0, [0.0,ds], stopevent = (s,u) -> (domain(u) >= 1) );
# uf = u[end];
#
#
# while domain(uf[1],uf[2]) < 1
#   s,u = ODE.ode45(metric,uf,[0.0,ds],stopevent = (s,u) -> ( u[1]^2 + u[2]^2 > 1) );
#   uf = u[end];
# end
#
# u0 = u[end-1];
# s1 = s[end-1][1]; s2 = s[end][1]; spread = s2-s1;
#
# # Need more data points (specifically 5) for 4th order:
# dt = spread/4.0;
# tspan = 0:dt:spread;
# s,u = ODE.ode45(metric,u0,tspan);
# up = u[1:5]; sp = s[1:5] + s1;
# up[end] = uf; sp[end] = s2;
# xp = map( v -> v[1], up); yp = map( v-> v[2], up);
# xI = polyfit(sp,xp); yI = polyfit(sp, yp);
# dxI = polyder(xI); dyI = polyder(yI);
#
# sdomain(s) = domain(xI(s),yI(s));
# sdomaingrad(s) = dot(domaingrad(xI(s),yI(s)) , [dxI(s), dyI(s)] );
# sout = newtonbisection(sdomain, sdomaingrad, s1, s2, 1e-5);
#
# uf[1] = xI(sout); uf[2] = yI(sout);
#
# # Interpolate velocities:
# for k = 1:length(uf)-2
#   vp = map( t -> t[k+2], up);
#   v = polyfit(sp,vp);
#   uf[k+2] = v(sout);
# end
#
# return uf;
#
# end

###########################################################################

# Is it possible to define an arbitrary method that gathers
# the full scattering relation? Because of finding boundary points....

function fullScatteringRelation(metric::Function,domain::Function,domaingrad::Function, Nedge, Nangle)
# Here, N should just be how many points on the boundary in general.

ds=1;
dphi = pi/Nangle;
uTotExit = Array{Any}(Nedge, Nangle);
# For the cells, Row ~ point on the boundary edge, Collumn ~ angle of incidence.

for i = 1:Nedge-1
    for j = 1:Nangle-1
      # 0) Define ds so that the evolution is fine enough at boundary?
      # 1) Find a way to parameterize the points on the boundary regardless of
      #    whatever boundary funtion is inputted....
      # 2) Evolve using those parameterizations (also need to get right angles!)


      # Example on Circle:
      # ds = 2*j*dphi*(Nangle-j)/Nangle;
      # u0 = [cos(i*dtheta), sin(i*dtheta), cos(i*dtheta + pi/2 + j*dphi), sin(i*dtheta + pi/2 + j*dphi)];
      # uTotExit[i,j] = scatteringRelation(metric,domain,domaingrad,u0,ds);
    end
end

return uTotExit;

end

###############################################################


function newtonbisection(f,df,a::Float64,b::Float64,tol)
# Newton solver for the boundary exit.
# Combines bisection too, since we know in the range of (s1,s2) there's a root.

maxiter = 100;
p=a; iter = 1;
while abs(f(p)[1]) > tol && iter < 100
  iter = iter + 1;
   p = p - f(p)/df(p);
   if p>b
       p=(a+b)/2; end;
   if p<a
       p=(a+b)/2; end;
   if f(p)*f(b)<0
       a=p;
   else
       b=p;
   end
end

return p;

end
