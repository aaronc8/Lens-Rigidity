## Functions used to solve for the Riemannian Metric.

function GradHessFFT(cxy::Array{Float64,2},p0,pf)
  # This function computes the gradient and Hessian of the wavespeedn cxy
  #
  ## TODO: what is the signature of this funtion?
  ## what's is p0 and pf?
## Solve for the derivatives of just c with fft. We can square c or do whatever
## to it from there.

  fftcxy = fft(cxy[1:end-1, 1:end-1]);
  fftcxy = fftshift(fftcxy);
  hy,hx = size(fftcxy);   ## We will only return the middle two quarters. WHY...?
  fftgradcxy = zeros(Complex{Float64},hy,hx,2);

# x,y both in [p0,pf]
  W = 0.5*hx; W = convert(Int64,W);
  for k = -W:1:W-1  # yrange
      for m = -W:1:W-1 # xrange
          fftgradcxy[k+W+1,m+W+1,1] = fftcxy[k+W+1,m+W+1]*2im*pi*m/(pf-p0);
          fftgradcxy[k+W+1,m+W+1,2] = fftcxy[k+W+1,m+W+1]*2im*pi*k/(pf-p0);
      end
  end
  fftgradcxy[:,:,1] = ifft(ifftshift(fftgradcxy[:,:,1]));
  fftgradcxy[:,:,2] = ifft(ifftshift(fftgradcxy[:,:,2]));

  ffthesscxy = zeros(Complex{Float64},hy,hx,2,2);
  for k = -W:1:W-1  # yrange
      for m = -W:1:W-1 # xrange
          ffthesscxy[k+W+1,m+W+1,1,1] = fftgradcxy[k+W+1,m+W+1,1]*2im*pi*m/(pf-p0);
          ffthesscxy[k+W+1,m+W+1,1,2] = fftgradcxy[k+W+1,m+W+1,1]*2im*pi*k/(pf-p0);
          ffthesscxy[k+W+1,m+W+1,2,1] = fftgradcxy[k+W+1,m+W+1,2]*2im*pi*m/(pf-p0);
          ffthesscxy[k+W+1,m+W+1,2,2] = fftgradcxy[k+W+1,m+W+1,2]*2im*pi*k/(pf-p0);
      end
  end
  ffthesscxy[:,:,1,1] = ifft(ifftshift(ffthesscxy[:,:,1,1]));
  ffthesscxy[:,:,1,2] = ifft(ifftshift(ffthesscxy[:,:,1,2]));
  ffthesscxy[:,:,2,1] = ifft(ifftshift(ffthesscxy[:,:,2,1]));
  ffthesscxy[:,:,2,2] = ifft(ifftshift(ffthesscxy[:,:,2,2]));

  r0 = convert(Int64, hx/4 + 1); rf = 3*r0 + -2;
  ## WHYYYY....???!!!

  return real(fftgradcxy[r0:rf, r0:rf, : ]), real(ffthesscxy[r0:rf, r0:rf, :, :]);
  # return fftgradcxy,ffthesscxy
end

#################################################

function GradHessFinDiff(cxy::Array{Float64,2})
  # Function to compute a finite difference approximation of the gradient and the
  # Hessian. (this should be low order)
# For derivatives, an alternative is to just use the interior points, too.

  ## We can add higher order derivates usien the Fornberg algoritun to produce the
  # stencils
  Nedge = length(cxy[:,1]) - 1;
  gradcxy = zeros(Nedge+1,Nedge+1,2);
  hesscxy = zeros(Nedge+1,Nedge+1,2,2);

  k = 2:Nedge;

  gradcxy[:,k,2] = (Nedge/4).*(cxy[:,k+1] - cxy[:,k-1]);
  gradcxy[:,1,2] = (Nedge/4).*(-3.0.*cxy[:,1] + 4.0.*cxy[:,2] - cxy[:,3]);
  gradcxy[:,end,2] = (-Nedge/4).*(-3.0.*cxy[:,end] + 4.0.*cxy[:,end-1] - cxy[:,end-2]);   #####
  gradcxy[k,:,1] = (Nedge/4).*(cxy[k+1,:] - cxy[k-1,:]);
  gradcxy[1,:,1] = (Nedge/4).*(-3.0.*cxy[1,:] + 4.0.*cxy[2,:] - cxy[3,:]);
  gradcxy[end,:,1] = (-Nedge/4).*(-3.0.*cxy[end,:] + 4.0.*cxy[end-1,:] - cxy[end-2,:]);   #####

  # For Hessian, just reuse the gradient this time since it uses cxy pts.
  hesscxy[:,k,1,2] = (Nedge/4).*(gradcxy[:,k+1,1] - gradcxy[:,k-1,1]);
  hesscxy[:,1,1,2] = (Nedge/4).*(-3.0.*gradcxy[:,1,1] + 4.0.*gradcxy[:,2,1] - gradcxy[:,3,1]);
  hesscxy[:,end,1,2] = (-Nedge/4).*(-3.0.*gradcxy[:,end,1] + 4.0.*gradcxy[:,end-1,1] - gradcxy[:,end-2,1]);
  hesscxy[k,:,1,1] = (Nedge/4).*(gradcxy[k+1,:,1] - gradcxy[k-1,:,1]);
  hesscxy[1,:,1,1] = (Nedge/4).*(-3.0.*gradcxy[1,:,1] + 4.0.*gradcxy[2,:,1] - gradcxy[3,:,1]);
  hesscxy[end,:,1,1] = (-Nedge/4).*(-3.0.*gradcxy[end,:,1] + 4.0.*gradcxy[end-1,:,1] - gradcxy[end-2,:,1]);

  hesscxy[:,k,2,2] = (Nedge/4).*(gradcxy[:,k+1,2] - gradcxy[:,k-1,2]);
  hesscxy[:,1,2,2] = (Nedge/4).*(-3.0.*gradcxy[:,1,2] + 4.0.*gradcxy[:,2,2] - gradcxy[:,3,2]);
  hesscxy[:,end,2,2] = (-Nedge/4).*(-3.0.*gradcxy[:,end,2] + 4.0.*gradcxy[:,end-1,2] - gradcxy[:,end-2,2]);
  hesscxy[k,:,2,1] = (Nedge/4).*(gradcxy[k+1,:,2] - gradcxy[k-1,:,2]);
  hesscxy[1,:,2,1] = (Nedge/4).*(-3.0.*gradcxy[1,:,2] + 4.0.*gradcxy[2,:,2] - gradcxy[3,:,2]);
  hesscxy[end,:,2,1] = (-Nedge/4).*(-3.0.*gradcxy[end,:,2] + 4.0.*gradcxy[end-1,:,2] - gradcxy[end-2,:,2]);

return gradcxy,hesscxy;

end

###################################################

function generateMetric(knots,cxy::Array{Float64,2},
                              dcxy::Array{Float64,3},
                              d2cxy::Array{Float64,4};
                              interpType::String="linear")
# Takes the computed values of cxy at (x,y) ∈ knots and interpolates.
# input : knots:  physical points in which the values of
#                 the functions are defined
#         cxy:    Array containing the value of the wavespeed
#                 at the points defined by knots
#         dcxy:   Array containing the value of the gradient of
#                 the wavespeed at the points defined by knots
#         d2cxy:  Array containing the value of the Hessian of
#                 the wavespeed at the points defined by knots
#         interpType: type of interpolation
# outputs: cspd:     function that interpolated the wavespeed
#          gradcspd: function that interpolated the gradient of
#                    the wavespeed
#          hesscspd: function that interpolated the hessian of
#                    the wavespeed

### For other domains, if BSplines, may need the length of box it is contained in
### due to the BSpline being off the knots.

  # gradmetric(x,y) = gradient(itp,x,y);
  if interpType == "linear"
    ## Knotted Approach:
    cInterp = interpolate(knots, cxy, Gridded(Linear()));
    # this gives us back an interpolation object
    cspd(x,y) = cInterp[x,y];

    # interpolation of the gradient
    dcdxInterp = interpolate(knots,dcxy[:,:,1],Gridded(Linear()));
    dcdyInterp = interpolate(knots,dcxy[:,:,2],Gridded(Linear()));
    gradcspd(x,y) = [dcdxInterp[x,y],dcdyInterp[x,y]];

    # interpolation of the Hessian
    d2cdx2Interp = interpolate(knots,d2cxy[:,:,1,1],Gridded(Linear()));
    d2cdxyInterp = interpolate(knots,d2cxy[:,:,1,2],Gridded(Linear()));
    d2cdyxInterp = interpolate(knots,d2cxy[:,:,2,1],Gridded(Linear()));
    d2cdy2Interp = interpolate(knots,d2cxy[:,:,2,2],Gridded(Linear()));
    hesscspd(x,y) = [d2cdx2Interp[x,y] d2cdxyInterp[x,y] ;
                     d2cdyxInterp[x,y] d2cdy2Interp[x,y]];

  elseif interpType == "cubic"
    # println("very experimental! beware things may not work as intended")
    # ##TODO: use knots with the information of the mesh and then
    # ## use the rescaling follwoing knots
    # c = interpolate(cxy,BSpline(Cubic(Natural())), OnGrid())
    # Nx,Ny = size(cxy);
    # # Renormalize to a unit interval in both dimensions, too:
    # cspd(x,y) = c[0.5*x*(Nx-1) + 0.5*(Nx+1) , 0.5*y*(Ny-1) + 0.5*(Ny+1)];

    # # interpolation of the gradient
    # dcdx = interpolate(dcxy[:,:,1],BSpline(Cubic(Natural())),OnGrid());
    # gradx(x,y) = dcdx[0.5*x*(Nx-1) + 0.5*(Nx+1), 0.5*y*(Ny-1) + 0.5*(Ny+1)]
    # dcdy = interpolate(dcxy[:,:,2],BSpline(Cubic(Natural())),OnGrid());
    # grady(x,y) = dcdy[0.5*x*(Nx-1) + 0.5*(Nx+1), 0.5*y*(Ny-1) + 0.5*(Ny+1)];

    # gradcspd(x,y) = [gradx(x,y),grady(x,y) ];

    # # interpolation of the Hessian
    # d2cdxx = interpolate(d2cxy[:,:,1,1],BSpline(Cubic(Natural())),OnGrid());
    # hessxx(x,y) = d2cdxx[0.5*x*(Nx-1) + 0.5*(Nx+1), 0.5*y*(Ny-1) + 0.5*(Ny+1)];
    # d2cdxy = interpolate(d2cxy[:,:,1,2],BSpline(Cubic(Natural())),OnGrid());
    # hessxy(x,y) = d2cdxy[0.5*x*(Nx-1) + 0.5*(Nx+1), 0.5*y*(Ny-1) + 0.5*(Ny+1)];
    # d2cdyx = interpolate(d2cxy[:,:,2,1],BSpline(Cubic(Natural())),OnGrid());
    # hessyx(x,y) = d2cdyx[0.5*x*(Nx-1) + 0.5*(Nx+1), 0.5*y*(Ny-1) + 0.5*(Ny+1)];
    # d2cdyy = interpolate(d2cxy[:,:,2,2],BSpline(Cubic(Natural())),OnGrid());
    # hessyy(x,y) = d2cdyy[0.5*x*(Nx-1) + 0.5*(Nx+1), 0.5*y*(Ny-1) + 0.5*(Ny+1)];

    # hesscspd(x,y) = [ hessxx(x,y) hessxy(x,y);
    #                   hessyx(x,y) hessyy(x,y)];
  end

  return cspd,gradcspd,hesscspd;

end

###############################


function makeHamiltonian(metric::Function,dmetric::Function, theta::Bool)
  # NOTE: in this case, I think that theta needs to be an optional parameter
  # with a default value
# function that gives back a the derivative of the Hamiltonian
# I'm not sure if this is fast though, I'm adding some typing
  @inline function dH(s::Float64,u::Array{Float64,1})
    V = zeros(Float64,4);
    cons = metric(u[1],u[2]);
    # position
    V[1:2] = u[3:4]*cons
    cons2 = -0.5.*(u[3].^2 + u[4].^2)
    # momentum
    V[3:4] = cons2*dmetric(u[1],u[2])
    return V;
  end

  function dHtheta(s::Float64,u::Array{Float64,1})
    Vtheta = zeros(Float64,3);
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

###############################
function generatePath(bdy::Function, metric::Function,u0,ds)
# Generate the paths with any metric.
# Interpolate the Ray evolution so that we have the path as a parametric function.
# Similar to method "squaregaussianrelation" but we have to keep the whole path.
# We still need to solve for exit time s_out for upper bound of integral....

    # solve the geodesic using the ODE solver
    s,u = ode45(metric, u0, [0.0,ds]);

    # we check if the ray is out of the domain
    k = find(x -> ( bdy(x[1],x[2]) >= -eps() ), u[2:end]);
    kshift = 0;

    # We make sure that the ray is out of the domain!
    scaling = max(ds,1/ds);
    maxiter = 1;

    # we interate until the ray if fully out of the domain
    while isempty(k) && maxiter < 1000*scaling ## In case of trapped rays
        u0 = u[end];
        sadd,uadd = ode45(metric,u0,[0.0,ds]);
        k = find(x -> ( bdy(x[1],x[2]) >= -eps() ), uadd[2:end]);
        sadd = sadd + s[end];
        #shift = length(u) - length(uadd);
        #kshift = length(u) - shift;
        kshift = length(u);
        s = append!(s,sadd);
        u = append!(u,uadd);
        maxiter = maxiter + 1;
    end


    if isempty(k)
      k = length(u);
    else
      k = k[1] + 1 + kshift;
    end

    ## Now to give the interpolated path.

    # define the interpolation values
    knots = (s[1:k],);

    # NOTE: I'm not sure that s are actually equidistant

    # extract the first and second componene of the phase
    # space vector
    xval = map( t -> t[1], u[1:k]);
    yval = map( t -> t[2], u[1:k]);

    # define the interpolation objects using the package
    xitp = interpolate(knots, xval, Gridded(Linear()));
    yitp = interpolate(knots, yval, Gridded(Linear()));

    # to interpolate the momentum, it will depend if we are in
    # cartesian of polar coordinates

    # Xg(s) = zeros(length(u0));
    # Polar coordinates
    if length(u0) == 3   ## This is being overwritten by the other if statement still,,?
      thetaval = map(t -> t[3], u[1:k]);
      thetaitp = interpolate(knots,thetaval,Gridded(Linear()));
      Xgtheta(s) = [xitp[s], yitp[s], thetaitp[s]];
    end

    # cartesian coordinates
    if length(u0) == 4
      vxval = map(t -> t[3], u[1:k]);
      vyval = map(t -> t[4], u[1:k]);
      vxitp = interpolate(knots, vxval, Gridded(Linear()));
      vyitp = interpolate(knots, vyval, Gridded(Linear()));
      Xg(s) = [xitp[s], yitp[s], vxitp[s], vyitp[s]];
    end

    ### WHY DO WE NEED THIS???

    ## Construct interpolant for exit time:    ## Maybe we don't need to? It should be
    ## given from the original ODE evolution...

    u1 = u[k-1]; u2 = u[k];
    s1 = s[k-1][1]; s2 = s[k][1]; spread = s2-s1; sout = s1+0.5*spread;

    if spread <= 1e-4 || norm(u1 - u2, 2) <= 1e-4
      if length(u1)  == 3
        return Xgtheta,sout;
      end
      if length(u1) == 4
        return Xg,sout;
      end
    end

    # u1 = u[end-1]; u2 = uf;
    # s1 = s[end-1][1]; s2 = s[end][1]; sf = s2-s1;

    # Need more data points (specifically 5) for 4th order:
    dt = spread/4.0;
    tspan = 0:dt:spread;
    # sp,up = ODE.ode45(gaussianmetric,u1,tspan);
    s,u = ode45(metric,u1,tspan,points=:specified);    # tspan only ensures at least 5 pts
    ## Need to find the times because otherwise, Least Squares is singular.
    # k2 = find(x -> x == tspan[2],s);
    # k3 = find(x -> x == tspan[3],s);
    # k4 = find(x -> x == tspan[4],s);
    # index = [1 k2[1] k3[1] k4[1] length(s)];

    # up = u[index]; sp = s[index] + s1;
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

function HamiltonianHess(cxy::Function, Dc::Function, D2c::Function)
# Form the Hessian (it is almost but not quite the Hessian) of H to
# Use to compute/solve for the Jacobian of the geodesic. (This is M)
# Didn't see a way to do this in theta perspective....

## BE Careful about the ordering!! Perhaps, go the collumnwise convention?
  function JacobianRHS(u)

    M = zeros(4,4);
    ## Build so that if we did reshape(M,4,4) it's the Hessian.
    spd = cxy(u[1],u[2]);
    grad = Dc(u[1],u[2]);
    hess = D2c(u[1],u[2]);

    # First 4 are H_x,v

    M[1] = 2*u[3].*grad[1].*spd;
    M[2] = 2*u[4].*grad[1].*spd;
    M[5] = 2*u[3].*grad[2].*spd;
    M[6] = 2*u[4].*grad[2].*spd;
    # Next 4 are -H_x,x
    M[3] = -(grad[1].^2 + spd.*hess[1,1]).*(u[3].^2 + u[4]).^2;
    M[4] = -(grad[1].*grad[2] + spd.*hess[2,1]).*(u[3].^2 + u[4]).^2;
    M[7] = -(grad[1].*grad[2] + spd.*hess[1,2]).*(u[3].^2 + u[4]).^2;  # M[4];  # Should equal...
    M[8] = -(grad[2].^2 + spd.*hess[2,2]).*(u[3].^2 + u[4]).^2;
    # Next 4 are H_v,v
    M[9] = spd.^2;
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

  return JacobianRHS;

end

################################

function geodesicJacobian(M::Function,X::Function,J0,sout)
# Solve for the Jacobian matrix of geodesic w.r.t. initial condition.
# J0 is usually identity, but it could be something else (e.g. "B" if broken)
# Since we have the Group Property, I don't think we need to specify initial s?
# Nor do we need J to be a function of X because it's inherent in M(X).

# Define the source term
@inline function F(s::Float64,J::Array{Float64,1})
  n = length(J0[1,:]);
  J = reshape(J,n,n);   # If already square, does nothing.
  return (M(X(s))*J)[:];   # This is the RHS of the ODE for J, also in vector form.
end

# run the ode solver
s,J = ode45(F,J0[:],[0.0,sout]);
knots = (s,);

#NOTE: this should be encapsulated in a type
# it makes a new interp object
@inline function Jacobian(t::Float64)
  Jac = zeros(size(J0));
  for k = 1:length(J0)   # Maybe Interpolations.jl can actually interpolate vector valued fns, but not sure, thought it failed before.
    Jcomp = map(a -> a[k], J);
    itpk = interpolate(knots,Jcomp,Gridded(Linear()));   # Just use knots + linear, simpler, o.o
    # IF we want to do a better interpolation:
    Jac[k] = itpk[t];
  end
  return Jac;
end

return Jacobian;

end

#################################

function lowerNeighbor(t,Nint,a,b)
# Assuming the interval [a,b] is divided into Nint+1 points, hence N intervals, find
# the lowest neighbor index and the weight/ratio of it.
# Used for interpolating the gridpoints for the update system.

  floatindex = (t-a)*Nint/(b-a) + 1;
  index = Int64(trunc(floatindex));
  weight = 1 - (floatindex - index);
  # This way, we'd have like x[index]*weight + x[index+1]*(1-weight) = t.

  answer = Array{Any}(2,1);
  answer[1] = index; answer[2] = weight;
  return answer;

end


function lowerNeighbor2D(x,y,Nxint,Nyint,a,b,c,d)
  # 2 dimensional version of the above, if convenient.
  answer = Array{Any}(2,2);
  answer[1,:] = lowerNeighbor(x,Nxint,a,b);
  answer[2,:] = lowerNeighbor(y,Nyint,c,d);
  return answer;

end


function rowDependence(kx,ky,wx,wy,Nedge)
# Construct the row in the system solving for the update.
# Sum of all the entries should be 1
    row = zeros(Nedge+1,Nedge+1);  # In case it starts on the Top or Right of square.
    if kx < 1  # Know it is really just along the West Edge.
      kx = 1; wx = 1;
    end
    if kx > Nedge # Know it is really along the East Edge.
      kx = Nedge; wx = 0;
    end
    if ky < 1 # Know it is along the South Edge.
      ky = 1; wy = 1;
    end
    if ky > Nedge  # Know it is along the North Edge.
      ky = Nedge; wy = 0;
    end
    row[kx,ky] = wx*wy;
    row[kx,ky+1] = wx*(1-wy);
    row[kx+1,ky] = (1-wx)*wy;
    row[kx+1,ky+1] = (1-wx)*(1-wy);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row,1,(Nedge+1)^2)
    return row;

end

## For derivatives: Check by sum going to 0.
function rowDependenceDx(kx,ky,wx,wy,Nedge,a,b)
# Construct the row in the system solving for the update.

  row = zeros(Nedge+2,Nedge+2);
  h = (b-a)/Nedge;
  if ky < 1 # Know it is along the South Edge.
    ky = 1; wy = 1;
  end
  if ky > Nedge  # Know it is along the North Edge.
    ky = Nedge; wy = 0;
  end

  if (kx > 1 && kx < Nedge)
    row[kx-1,ky] = -wx*wy/(2*h);
    row[kx+1,ky] = wx*wy/(2*h);
    row[kx,ky] = -(1-wx)*wy/(2*h);
    row[kx+2,ky] = (1-wx)*wy/(2*h);

    row[kx-1,ky+1] = -wx*(1-wy)/(2*h);
    row[kx+1,ky+1] = wx*(1-wy)/(2*h);
    row[kx,ky+1] = -(1-wx)*(1-wy)/(2*h);
    row[kx+2,ky+1] = (1-wx)*(1-wy)/(2*h);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row,1,(Nedge+1)^2)
    return row;
  end

  if kx <= 1
    kx=1;
    if kx < 1
      wx = 1;
    end
    row[kx,ky] = -3*wx*wy/(2*h);
    row[kx+1,ky] = 4*wx*wy/(2*h);
    row[kx+2,ky] = -wx*wy/(2*h);
    row[kx,ky] = row[kx,ky] - (1-wx)*wy/(2*h);
    row[kx+2,ky] = row[kx+2,ky] + (1-wx)*wy/(2*h);

    row[kx,ky+1] = -3*wx*(1-wy)/(2*h);
    row[kx+1,ky+1] = 4*wx*(1-wy)/(2*h);
    row[kx+2,ky+1] = -wx*(1-wy)/(2*h);
    row[kx,ky+1] = row[kx,ky+1] - (1-wx)*(1-wy)/(2*h);
    row[kx+2,ky+1] = row[kx+2,ky+1] + (1-wx)*(1-wy)/(2*h);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row,1,(Nedge+1)^2)
    return row;
  end

  if kx >= Nedge
    kx = Nedge;
    if kx > Nedge
      wx = 0;
    end
    row[kx+1,ky] = -3*wx*wy/(-2*h);
    row[kx,ky] = 4*wx*wy/(-2*h);
    row[kx-1,ky] = -wx*wy/(-2*h);
    row[kx,ky] = row[kx,ky] - (1-wx)*wy/(2*h);
    row[kx-1,ky] = row[kx-1,ky] + (1-wx)*wy/(2*h);

    row[kx+1,ky+1] = -3*wx*(1-wy)/(-2*h);
    row[kx,ky+1] = 4*wx*(1-wy)/(-2*h);
    row[kx-1,ky+1] = -wx*(1-wy)/(-2*h);
    row[kx,ky+1] = row[kx,ky+1] - (1-wx)*(1-wy)/(2*h);
    row[kx-1,ky+1] = row[kx-1,ky+1] + (1-wx)*(1-wy)/(2*h);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row,1,(Nedge+1)^2)
    return row;
  end

end

function rowDependenceDy(kx,ky,wx,wy,Nedge,a,b)
# Construct the row in the system solving for the update.

  if kx < 1  # Know it is really just along the West Edge.
    kx = 1; wx = 1;
  end
  if kx > Nedge # Know it is really along the East Edge.
    kx = Nedge; wx = 0;
  end

  storek = kx; kx = ky; ky = storek;  # Because we can just reuse the idea before.
  storew = wx; wx = wy; wy = storew;

  row = zeros(Nedge+2,Nedge+2);
  h = (b-a)/Nedge;

  if (kx > 1 && kx < Nedge)
    row[kx-1,ky] = -wx*wy/(2*h);
    row[kx+1,ky] = wx*wy/(2*h);
    row[kx,ky] = -(1-wx)*wy/(2*h);
    row[kx+2,ky] = (1-wx)*wy/(2*h);

    row[kx-1,ky+1] = -wx*(1-wy)/(2*h);
    row[kx+1,ky+1] = wx*(1-wy)/(2*h);
    row[kx,ky+1] = -(1-wx)*(1-wy)/(2*h);
    row[kx+2,ky+1] = (1-wx)*(1-wy)/(2*h);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row',1,(Nedge+1)^2)
    return row;
  end

  if kx <= 1
    kx = 1;
    if kx < 1
      wx = 1;
    end
    row[kx,ky] = -3*wx*wy/(2*h);
    row[kx+1,ky] = 4*wx*wy/(2*h);
    row[kx+2,ky] = -wx*wy/(2*h);
    row[kx,ky] = row[kx,ky] - (1-wx)*wy/(2*h);
    row[kx+2,ky] = row[kx+2,ky] + (1-wx)*wy/(2*h);

    row[kx,ky+1] = -3*wx*(1-wy)/(2*h);
    row[kx+1,ky+1] = 4*wx*(1-wy)/(2*h);
    row[kx+2,ky+1] = -wx*(1-wy)/(2*h);
    row[kx,ky+1] = row[kx,ky+1] - (1-wx)*(1-wy)/(2*h);
    row[kx+2,ky+1] = row[kx+2,ky+1] + (1-wx)*(1-wy)/(2*h);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row',1,(Nedge+1)^2)
    return row;
  end

  if kx >= Nedge
    kx = Nedge;
    if kx > Nedge
      wx = 0;
    end
    row[kx+1,ky] = -3*wx*wy/(-2*h);
    row[kx,ky] = 4*wx*wy/(-2*h);
    row[kx-1,ky] = -wx*wy/(-2*h);
    row[kx,ky] = row[kx,ky] - (1-wx)*wy/(2*h);
    row[kx-1,ky] = row[kx-1,ky] + (1-wx)*wy/(2*h);

    row[kx+1,ky+1] = -3*wx*(1-wy)/(-2*h);
    row[kx,ky+1] = 4*wx*(1-wy)/(-2*h);
    row[kx-1,ky+1] = -wx*(1-wy)/(-2*h);
    row[kx,ky+1] = row[kx,ky+1] - (1-wx)*(1-wy)/(2*h);
    row[kx-1,ky+1] = row[kx-1,ky+1] + (1-wx)*(1-wy)/(2*h);

    row = row[1:Nedge+1, 1:Nedge+1];
    row = reshape(row',1,(Nedge+1)^2)
    return row;
  end

end

#
# #################################
#
# function dVfrechet(cspd, dcspd, lambda, dlambda)
# # Compute the frechet derivative of the Hamiltonian flow
# function dVdg(X)
#   dV = zeros(4,1);
#   c = cspd(X[1],X[2]);
#   dc = dcspd(X[1],X[2]);
#   L = lambda(X[1],X[2]);
#   dL = dlambda(X[1],X[2]);
#
#   dV[1] = 2.0.*c.*L.*X[3];
#   dV[2] = 2.0.*c.*L.*X[4];
#   dV[3] = -(X[3].^2 + X[4].^2).*(L.*dc[1] + c.*dL[1]);
#   dV[4] = -(X[3].^2 + X[4].^2).*(L.*dc[2] + c.*dL[2]);
#   return dV;
# end
# return dVdg;
# end

#################################

function linearMismatch(J,Xg,tf,c,gradc,Nedge)
# Compute the 4 rows of the matrix system, corresponding to the ray.

s = 0:tf/2^4:tf; h = tf/2^4; n = 2^4;   # This way 2^(-16) ~ 1e-5 error with Simpson.
# Kni = J(tf-0)*dV(Xg(0)) + J(0)*dV(Xg(tf));  # Endpoints of Simpson.
# step = Nedge^2;  # Useful if we solve for grad(c) too.

Rows = zeros(4,(Nedge+1)^2);  # To solve for c, dcdx, and dcdy on the grid, its 3*(Nedge+1)^2
# There's 4 rows, one for each component of ΔXg.
RowA = zeros(4,(Nedge+1)^2);
RowB = zeros(4,(Nedge+1)^2);
RowC = zeros(4,(Nedge+1)^2);

M = J(tf);

for j = 1:2:2^4
  A = inv(J(s[j]));    A = M*A*h/3;   uA = Xg(s[j]);
  B = inv(J(s[j+1])); B = 4*M*B*h/3; uB = Xg(s[j+1]);
  C = inv(J(s[j+2])); C = M*C*h/3;   uC = Xg(s[j+2]);

  Acorner = lowerNeighbor2D(uA[1],uA[2],Nedge,Nedge,-1,1,-1,1);
  Bcorner = lowerNeighbor2D(uB[1],uB[2],Nedge,Nedge,-1,1,-1,1);
  Ccorner = lowerNeighbor2D(uC[1],uC[2],Nedge,Nedge,-1,1,-1,1);

  colA = rowDependence(Acorner[1],Acorner[2],Acorner[3],Acorner[4],Nedge)[:];
  DxcolA = rowDependenceDx(Acorner[1],Acorner[2],Acorner[3],Acorner[4],Nedge,-1,1)[:];
  DycolA = rowDependenceDx(Acorner[1],Acorner[2],Acorner[3],Acorner[4],Nedge,-1,1)[:];

  colB = rowDependence(Bcorner[1],Bcorner[2],Bcorner[3],Bcorner[4],Nedge)[:];
  DxcolB = rowDependenceDx(Bcorner[1],Bcorner[2],Bcorner[3],Bcorner[4],Nedge,-1,1)[:];
  DycolB = rowDependenceDx(Bcorner[1],Bcorner[2],Bcorner[3],Bcorner[4],Nedge,-1,1)[:];

  colC = rowDependence(Ccorner[1],Ccorner[2],Ccorner[3],Ccorner[4],Nedge)[:];
  DxcolC = rowDependenceDx(Ccorner[1],Ccorner[2],Ccorner[3],Ccorner[4],Nedge,-1,1)[:];
  DycolC = rowDependenceDx(Ccorner[1],Ccorner[2],Ccorner[3],Ccorner[4],Nedge,-1,1)[:];

  for r = 1:4    ## Follow the multiplication of J*dV ## (Same for all 3 pieces though)
    cA = c(uA[1],uA[2]);  DcA = gradc(uA[1],uA[2]); vA = uA[3]^2 + uA[4]^2;
    RowA[r,:] = (A[r,1]*2*cA*uA[3] + A[r,2]*2*cA*uA[4] - vA*(A[r,3]*DcA[1] - A[r,4]*DcA[2]) )*colA;
    RowA[r,:] = RowA[r,:] - vA*cA*( A[r,3]*DxcolA + A[r,4]*DycolA );

    cB = c(uB[1],uB[2]);  DcB = gradc(uB[1],uB[2]); vB = uB[3]^2 + uB[4]^2;
    RowB[r,:] = (B[r,1]*2*cB*uB[3] + B[r,2]*2*cB*uB[4] - vB*(B[r,3]*DcB[1] - B[r,4]*DcB[2]) )*colB;
    RowB[r,:] = RowB[r,:] - vB*cB*( B[r,3]*DxcolB + B[r,4]*DycolB );

    cC = c(uC[1],uC[2]);  DcC = gradc(uC[1],uC[2]); vC = uC[3]^2 + uC[4]^2;
    RowC[r,:] = (C[r,1]*2*cA*uA[3] + C[r,2]*2*cC*uC[4] - vC*(C[r,3]*DcC[1] - C[r,4]*DcC[2]) )*colC;
    RowC[r,:] = RowC[r,:] - vC*cC*( C[r,3]*DxcolC + C[r,4]*DycolC );
  end
  Rows = Rows + RowA + RowB + RowC;
end

badindex = find(x -> abs(x) <= 1e-5, Rows);
Rows[badindex] = 0;

return Rows;

end


#####################################
## Put all the Kni together.

function MismatchSystem(bdy,cn::Function,gradcn::Function,
                            hesscn,Nedge,Nangle,bexact,sout)
# Gather the exit data mismatch for the square problem.
# The Ray Data should be in order E,N,W,S (counterclockwise).

metric(x,y) = cn(x,y).^2; dmetric(x,y) = 2.0.*cn(x,y).*gradcn(x,y);

dH = makeHamiltonian(metric,dmetric,false);

M = HamiltonianHess(cn,gradcn,hesscn);
J0 = eye(4,4);

dl = 2/Nedge;
dphi = pi/Nangle;
rowcount = 1;

# To solve A ̃g = b.
A = zeros(16*(Nedge-1)*(Nangle-1), (Nedge+1)^2);
b = zeros(16*(Nedge-1)*(Nangle-1));  # Will reshape to length of M's column.

# For the cells, each row is a point on the boundary edge and each collumn is
# an angle of incidence.

  for i = 1:Nedge-1
      ds = dl*2*i*(Nedge-i)/Nedge;   # For bad points near corners.
      for j = 1:Nangle-1
          uE0 = [1; -1+i*dl; cos(j*dphi + pi/2); sin(j*dphi + pi/2)];   # Only for right edge!
          XE,sEout = generatePath(bdy, dH, uE0, ds);
          JE = geodesicJacobian(M,XE,J0,sEout);
          RayData = linearMismatch(JE,XE,sout[rowcount],cn,gradcn,Nedge);
          L2err = norm(XE(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount],2)^2;
          # the residual is being weghted in here!
          A[4*(rowcount-1)+1:4*rowcount, :] = RayData./L2err;
          b[4*(rowcount-1)+1:4*rowcount] = (XE(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount])./L2err;
          rowcount = rowcount + 1;
      end
  end
  for i = 1:Nedge-1
      ds = dl*2*i*(Nedge-i)/Nedge;   # For bad points near corners.
      for j = 1:Nangle-1
          uN0 = [1 - i*dl; 1; cos(-j*dphi); sin(-j*dphi)];   # Only for top edge!
          XN,sNout = generatePath(bdy, dH, uN0, ds);
          JN = geodesicJacobian(M,XN,J0,sNout);
          RayData = linearMismatch(JN,XN,sout[rowcount],cn,gradcn,Nedge);
          L2err = norm(XN(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount],2)^2;
          A[4*(rowcount-1)+1:4*rowcount, :] = RayData./L2err;
          b[4*(rowcount-1)+1:4*rowcount] = (XN(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount])./L2err;
          rowcount = rowcount + 1;
      end
    end
    for i = 1:Nedge-1
        ds = dl*2*i*(Nedge-i)/Nedge;   # For bad points near corners.
        for j = 1:Nangle-1
          uW0 = [-1; 1-i*dl; cos(j*dphi - pi/2); sin(j*dphi - pi/2)];   # Only for left edge!
          XW,sWout = generatePath(bdy, dH, uW0, ds);
          JW = geodesicJacobian(M,XW,J0,sWout);
          RayData = linearMismatch(JW,XW,sout[rowcount],cn,gradcn,Nedge);
          L2err = norm(XW(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount],2)^2;
          A[4*(rowcount-1)+1:4*rowcount, :] = RayData./L2err;
          b[4*(rowcount-1)+1:4*rowcount] = (XW(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount])./L2err;
          rowcount = rowcount + 1;
        end
    end
    for i = 1:Nedge-1
        ds = dl*2*i*(Nedge-i)/Nedge;   # For bad points near corners.
        for j = 1:Nangle-1
          uS0 = [-1 + i*dl; -1; cos(j*dphi); sin(j*dphi)];   # Only for bottom edge!
          XS,sSout = generatePath(bdy, dH, uS0, ds);
          JS = geodesicJacobian(M,XS,J0,sSout);
          RayData = linearMismatch(JS,XS,sout[rowcount],cn,gradcn,Nedge);
          L2err = norm(XS(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount],2)^2;
          A[4*(rowcount-1)+1:4*rowcount, :] = RayData./L2err;
          b[4*(rowcount-1)+1:4*rowcount] = (XS(sout[rowcount]) - bexact[4*(rowcount-1)+1:4*rowcount])./L2err;
          rowcount = rowcount + 1;
      end
  end

return A,b;

end
