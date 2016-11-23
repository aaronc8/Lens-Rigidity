## FFT Practice in Julia Now
display("Absolute Errors of derivatives.")

## Basic Case:
n=3; N = 2^n; W = 2^(n-1);
xx = 0:1/N:1;
yy = sin(2*pi*xx);
sinfft = fft(yy[1:end-1]);
shiftsinfft = fftshift(sinfft);
dsinshift = shiftsinfft;
for k = -W:1:W-1
    dsinshift[k+W+1] = shiftsinfft[k+W+1]*2*1im*pi*k;
end
dsin = ifft(ifftshift(dsinshift));
exactdsin = 2*pi*cos(2*pi*xx);
display("Basic Case on [0,1]:")
display(dsin-exactdsin[1:end-1])


##################################################

## Slight 1-D Generalization:
n=3; N = 2^n; W = 2^(n-1); a = -2; b = 2;
xx = a:(b-a)/N:b;    # xx = xx./(b-a);
yy = cos(2*pi*xx);
sinfft = fft(yy[1:end-1]);
shiftsinfft = fftshift(sinfft);
dsinshift = shiftsinfft;

for k = -W:1:W-1
    dsinshift[k+W+1] = shiftsinfft[k+W+1]*2*1im*pi*k/(b-a);
# for k = W:-1:-W+1
#     dsinshift(-k+N/2+1) = shiftsinfft(-k+N/2+1)*2*1i*pi*k/(b-a);
end
dsin = ifft(ifftshift(dsinshift));
exactdsin = -2*pi*sin(2*pi*xx);
display("Slight Generalization to any [a,b]")
display(dsin-exactdsin[1:end-1])


###########################################
## 2-D case: (Basic for now)(Set up to be able to change box or pts)
n = 5; N = 2^n; W = 2^(n-1);
x1 = -5.0; x2 = 5.0; y1 = -7.0; y2 = 7.0;
xx = x1:(x2-x1)/N:x2; # xx = xx./(x2-x1);
yy = y1:(y2-y1)/N:y2; # yy = yy./(y2-y1);
# [xgrid,ygrid] = meshgrid(xx,yy);
zz = cos(3*pi*xx').*sin(2*pi*yy);  # Need to transpose x to preserve that y goes up / down
# zz = exp(xx'.^2 .+ yy.^2);

# prodfft = rfft(zz[1:end-1,1:end-1]);  ## Check phases....
# prodfft = fft(fftshift(fft(zz[1:end-1,1:end-1])).').';
prodfft = fft(zz[1:end-1,1:end-1]);
shiftprodfft = fftshift(prodfft); #.';
h,w = size(shiftprodfft);
dprodshift = zeros(Complex{Float64},h,w,2);

for k = -W:1:W-1  # yrange
    for m = -W:1:W-1 # xrange
        dprodshift[k+W+1,m+W+1,1] = shiftprodfft[k+W+1,m+W+1]*2im*pi*m/(x2-x1);
        dprodshift[k+W+1,m+W+1,2] = shiftprodfft[k+W+1,m+W+1]*2im*pi*k/(y2-y1);
    end
end
# dprod = ifft(ifftshift(ifft(ifftshift(dprodshift)).')).';
# dprod = ifft(ifftshift(dprodshift));
dprod = zeros(Complex{Float64},h,w,2);
# dprod[:,:,1] = ifft(ifftshift(ifft(ifftshift(dprodshift[:,:,1])).')).';  # Double checked, this does invert it.
# dprod[:,:,2] = ifft(ifftshift(ifft(ifftshift(dprodshift[:,:,2])).')).';
dprod[:,:,1] = ifft(ifftshift(dprodshift[:,:,1]));
dprod[:,:,2] = ifft(ifftshift(dprodshift[:,:,2]));

dxprod(x,y) = -3*pi*sin(3*pi*x).*sin(2*pi*y);
dyprod(x,y) = 2*pi*cos(3*pi*x).*cos(2*pi*y);
# dxprod(x,y) = 2*x.*exp(x.^2 .+ y.^2);
# dyprod(x,y) = 2*y.*exp(x.^2 .+ y.^2);

exactgrad = zeros(h+1,w+1,2);
exactgrad[:,:,1] = dxprod(xx',yy);
exactgrad[:,:,2] = dyprod(xx',yy);
#  display(dprod - exactgrad[1:end-1,1:end-1,:])

display("Slight generalization in 2D with (x,y) in [a,b] x [c,d]")
display(norm(dprod[:,:,1]-exactgrad[1:end-1,1:end-1,1],Inf))
display(norm(dprod[:,:,2]-exactgrad[1:end-1,1:end-1,2],Inf))



####################################################
## Practice + Storage of fft code in Ray Tracing:

## Check without extending the domain:
N = 5;
Nedge = 2^N;   # number of partitions of sides.

x = -1:1/Nedge:1;
y = x;
L = length(x);
cxy = ones(L,L);
# cxy= exp(0.5.*x'.^2 .+ 0.5.*y.^2);
cxy = cos(2*pi*x').*sin(5*pi*y);
# For a circle - see which points are inside and then assign a value.
dcxy = ones(L,L,2);
# dcxy[:,:,1] = x'.*exp(0.5.*x'.^2 .+ 0.5.*y.^2);
# dcxy[:,:,2] = y.*exp(0.5.*x'.^2 .+ 0.5.*y.^2);
dcxy[:,:,1] = -2*pi*sin(2*pi*x').*sin(5*pi*y);
dcxy[:,:,2] = 5*pi*cos(2*pi*x').*cos(5*pi*y);
d2cxy = ones(L,L,2,2);   ## If i want to test out the hessian. Should be fine if gradient works.


fftcxy = fft(cxy[1:end-1, 1:end-1]);
fftcxy = fftshift(fftcxy);
hy,hx = size(fftcxy);
fftgradcxy = ones(Complex{Float64},hy,hx,2);

x2 = 1.; y2=1.; x1=-1.; y1 = -1.; W = 0.5*(length(x)-1); W = convert(Int64,W);
for k = -W:1:W-1  # yrange
    for m = -W:1:W-1 # xrange
        fftgradcxy[k+W+1,m+W+1,1] = fftcxy[k+W+1,m+W+1]*2im*pi*m/(x2-x1);
        fftgradcxy[k+W+1,m+W+1,2] = fftcxy[k+W+1,m+W+1]*2im*pi*k/(y2-y1);
    end
end
fftgradcxy[:,:,1] = ifft(ifftshift(fftgradcxy[:,:,1]));
fftgradcxy[:,:,2] = ifft(ifftshift(fftgradcxy[:,:,2]));

ffthesscxy = ones(Complex{Float64},hy,hx,2,2);
for k = -W:1:W-1  # yrange
    for m = -W:1:W-1 # xrange
        ffthesscxy[k+W+1,m+W+1,1,1] = fftgradcxy[k+W+1,m+W+1,1]*2im*pi*m/(x2-x1);
        ffthesscxy[k+W+1,m+W+1,1,2] = fftgradcxy[k+W+1,m+W+1,1]*2im*pi*k/(y2-y1);
        ffthesscxy[k+W+1,m+W+1,2,1] = fftgradcxy[k+W+1,m+W+1,2]*2im*pi*m/(x2-x1);
        ffthesscxy[k+W+1,m+W+1,2,2] = fftgradcxy[k+W+1,m+W+1,2]*2im*pi*k/(y2-y1);
    end
end
ffthesscxy[:,:,1,1] = ifft(ifftshift(ffthesscxy[:,:,1,1]));
ffthesscxy[:,:,1,2] = ifft(ifftshift(ffthesscxy[:,:,1,2]));
ffthesscxy[:,:,2,1] = ifft(ifftshift(ffthesscxy[:,:,2,1]));
ffthesscxy[:,:,2,2] = ifft(ifftshift(ffthesscxy[:,:,2,2]));

display("If we don't splice into the square:")
display(norm(dcxy[1:end-1,1:end-1,:][:] - real(fftgradcxy[:]),Inf));

## Check when we do splice into square:
N = 5;
Nedge = 2^N;   # number of partitions of sides.

x = -2:1/Nedge:2;
y = x;
L = length(x);
cxy = ones(L,L);
# cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1] = exp(0.5.*x[Nedge+1: 3*Nedge+1]'.^2 .+ 0.5.*y[Nedge+1 : 3*Nedge+1].^2);
cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1] = cos(2*pi*x[Nedge+1: 3*Nedge+1]').*sin(5*pi*y[Nedge+1 : 3*Nedge+1]);
# For a circle - see which points are inside and then assign a value.
dcxy = ones(L,L,2);
# dcxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1,1] = x[Nedge+1: 3*Nedge+1]'.*exp(0.5.*x[Nedge+1: 3*Nedge+1]'.^2 .+ 0.5.*y[Nedge+1 : 3*Nedge+1].^2);
# dcxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1,2] = y[Nedge+1: 3*Nedge+1].*exp(0.5.*x[Nedge+1: 3*Nedge+1]'.^2 .+ 0.5.*y[Nedge+1 : 3*Nedge+1].^2);
dcxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1,1] = -2*pi*sin(2*pi*x[Nedge+1: 3*Nedge+1]').*sin(5*pi*y[Nedge+1 : 3*Nedge+1]);
dcxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1,2] = 5*pi*cos(2*pi*x[Nedge+1: 3*Nedge+1]').*cos(5*pi*y[Nedge+1 : 3*Nedge+1]);
d2cxy = ones(L,L,2,2);   ## If i want to test out the hessian. Should be fine if gradient works.


fftcxy = fft(cxy[1:end-1, 1:end-1]);
fftcxy = fftshift(fftcxy);
hy,hx = size(fftcxy);
fftgradcxy = ones(Complex{Float64},hy,hx,2);

x2 = 2.; y2=2.; x1=-2.; y1 = -2.; W = 0.5*(length(x)-1); W = convert(Int64,W);
for k = -W:1:W-1  # yrange
    for m = -W:1:W-1 # xrange
        fftgradcxy[k+W+1,m+W+1,1] = fftcxy[k+W+1,m+W+1]*2im*pi*m/(x2-x1);
        fftgradcxy[k+W+1,m+W+1,2] = fftcxy[k+W+1,m+W+1]*2im*pi*k/(y2-y1);
    end
end
fftgradcxy[:,:,1] = ifft(ifftshift(fftgradcxy[:,:,1]));
fftgradcxy[:,:,2] = ifft(ifftshift(fftgradcxy[:,:,2]));

ffthesscxy = ones(Complex{Float64},hy,hx,2,2);
for k = -W:1:W-1  # yrange
    for m = -W:1:W-1 # xrange
        ffthesscxy[k+W+1,m+W+1,1,1] = fftgradcxy[k+W+1,m+W+1,1]*2im*pi*m/(x2-x1);
        ffthesscxy[k+W+1,m+W+1,1,2] = fftgradcxy[k+W+1,m+W+1,1]*2im*pi*k/(y2-y1);
        ffthesscxy[k+W+1,m+W+1,2,1] = fftgradcxy[k+W+1,m+W+1,2]*2im*pi*m/(x2-x1);
        ffthesscxy[k+W+1,m+W+1,2,2] = fftgradcxy[k+W+1,m+W+1,2]*2im*pi*k/(y2-y1);
    end
end
ffthesscxy[:,:,1,1] = ifft(ifftshift(ffthesscxy[:,:,1,1]));
ffthesscxy[:,:,1,2] = ifft(ifftshift(ffthesscxy[:,:,1,2]));
ffthesscxy[:,:,2,1] = ifft(ifftshift(ffthesscxy[:,:,2,1]));
ffthesscxy[:,:,2,2] = ifft(ifftshift(ffthesscxy[:,:,2,2]));


cxy = cxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1];   ## Do same for gradients?
display("Splicing the square into a larger square, filling undefined spaces with 1s:")
display(norm(dcxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1,:][:] - real(fftgradcxy[Nedge+1:3*Nedge+1, Nedge+1:3*Nedge+1,:][:]),Inf));
