clear all ; clc ; close all
tic

%%--%% Parameters %%--%%

%Constants
hbar = 1.0545718e-34; %Planck
kbol = 1.38064852e-23; %Boltzmann
mu0 = 1.25663706e-6; %Permeability
muB = 9.274009994e-24; %Bohr magneton
a0 = 5.3e-11; %Bohr radius
m = 2.6983757e-25; %164Dy

%Number of grid points
Nx=256.;
Ny=256.;

%Length of box (dimensionless)
Lx=32;
Ly=32;

%Trap frequencies
wx = 2*pi*50;
wy = 2*pi*50;
wz = 2*pi*250;

%Units
w0 = wx;
l = sqrt(hbar/(m*wx)); %HO length
lz = sqrt(hbar/(m*wz)); %HO length


% Simulation inputs
as = 150 *a0/l; %scattering length
add = 130.8 *a0/l; %Dy
edd = add / as;

Om=0.8; %Dimensonless rotation frequency

Norm=20000; %atom number
dt=0.001; %timestep
tol=1e-10; %imaginary time tolerance

alpha=90 /360*2*pi; % Dipole angle (radians)
eta=0;% eta = 0 dipoles in z-x plane, pi/2 dipoles in z-y plane
sigma=lz/l; %ratio of lengthscales

%Trap gamma
gx=wx^2/w0^2;
gy=wy^2/w0^2;

%%--%% Calculating matrices %%--%%

Npx=Nx/2;
Npy=Ny/2;

bta=as*Norm;
hx=Lx./(2.*Npx-1); %x step size
deltakx=2*pi/(Lx); %momentum step size
x=(-Npx:1:(Npx-1)).*hx; %x array
kx=(-Npx:1:(Npx-1)).*deltakx; %momentum array
wavx=fftshift(kx);

hy=Ly./(2.*Npy-1); %y step size
deltaky=2*pi/(Ly); %momentum step size
y=(-Npy:1:(Npy-1)).*hy; %x array
ky=(-Npy:1:(Npy-1)).*deltaky; %momentum array
wavy=fftshift(ky);

wav_x=zeros(Nx,Ny);wav_y=wav_x;pot=wav_x;
ksx=fftshift(kx);
ksy=fftshift(ky);
[KSY,KSX]=meshgrid(ksy,ksx);
k=(KSX.^2 + KSY.^2).^(1/2)*sigma;
kssd=(cos(eta)*KSX+sin(eta)*KSY).^2*sigma^2;

[Y,X]=meshgrid(y,x);
%initialise
for i=1:Nx
    for j=1:Ny
        wav_x(i,j)=ksx(i);
        wav_y(i,j)=ksy(j);
        pot(i,j)=0.5*(gx*x(i)^2+gy*y(j)^2);
    end
end

% %Here are the 2D dipolar potentials required for the calculations
Udd_perp=2 - 3/sqrt(2)*sqrt(pi)*k.*erfcx(k/sqrt(2));
Udd_par=-1 + 3*sqrt(pi/32).*kssd./k.*erfcx(k/sqrt(2));

%solves issue with 0/0...
Udd_par(1,1)=-1;

Udd=cos(alpha)^2 * Udd_perp + sin(alpha)^2 * Udd_par;

%Initial vortex guess
u=sqrt(X.^2+Y.^2)./sqrt(X.^2+Y.^2 + 1/1.4).*exp(-pot/(6^2)).*exp(2*pi*1i*rand(size(pot)) + 1i*atan2(Y,X));

cur_norm=trapz(abs(u(:)).^2)*(hx*hy);
u=u/sqrt(cur_norm);

%%--%% Imaginary time

g=2*sqrt(2*pi)*bta/sigma;
gdd=2*sqrt(2*pi)*bta*edd/sigma;

[u_imag,mu]=ssfm_imag(u,x,y,dt,wav_x,wav_y,Udd,Y,X,pot,Om,g,gdd,tol);
