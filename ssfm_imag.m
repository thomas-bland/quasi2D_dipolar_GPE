function [u,mu]=ssfm_imag(u,x,y,dt,wav_x,wav_y,Udd,Y,X,pot,Om,g,gd,tol)

%Uses a split-step Fourier method to solve ground state quasi-2d dipolar eGPE

Nx=length(x);
Ny=length(y);
hx=x(2)-x(1);
hy=y(2)-y(1);

j=0;

mu = 0;
mu_err=1;
mu_old=1;

dt=-1i*abs(dt); %get rid of this line to make a real time simulation!
Phi=zeros(Nx,Ny);

while mu_err>tol    

    %kin x
    u=fft(u,Nx,1)/Nx;
    u=u.*exp(-0.25*1i*dt*wav_x.^2-0.5*dt*1i*Om*Y.*wav_x);
    u=ifft(u,Nx,1)*Nx;
     
    %kin y
    u=fft(u,Ny,2)/Ny;
    u=u.*exp(-0.25*1i*dt*wav_y.^2+0.5*dt*1i*Om*X.*wav_y);
    u=ifft(u,Ny,2)*Ny;
      
    %pot
    if gd~=0
        frho=fftn(abs(u).^2)/(Nx*Ny);
        Phi=real(ifftn(frho.*Udd))*(Nx*Ny);
    end
    u=u.*exp(-1i*dt*(pot + g*abs(u).^2 + gd*Phi - mu));

    %kin y
    u=fft(u,Ny,2)/Ny;
    u=u.*exp(-0.25*1i*dt*wav_y.^2+0.5*dt*1i*Om*X.*wav_y);
    u=ifft(u,Ny,2)*Ny;
   

    %kin x
    u=fft(u,Nx,1)/Nx;
    u=u.*exp(-0.25*1i*dt*wav_x.^2-0.5*dt*1i*Om*Y.*wav_x);
    u=ifft(u,Nx,1)*Nx;
    
    %renorm
    cur_norm=trapz(abs(u(:)).^2)*(hx*hy);
    u=u/sqrt(cur_norm);
    
    %calc mu error
    if gd~=0
        frho=fftn(abs(u).^2)/(Nx*Ny);
        Phi=real(ifftn(frho.*Udd))*(Nx*Ny);
    end
    Ekin = real(trapz(trapz(0.5*conj(u).*ifft2((wav_x.^2+wav_y.^2).*fft2(u))))*hx*hy);
    Eg = 0.5*trapz(trapz(g*abs(u).^4))*hx*hy;
    Edd = 0.5*trapz(trapz(gd*abs(u).^2.*Phi))*hx*hy;
    Erot = real(trapz(trapz(Om*conj(u).*(Y.*ifft2(wav_x.*fft2(u)) - X.*ifft2(wav_y.*fft2(u)))))*hx*hy);
    
    mu = Ekin + 2*Eg + 2*Edd + Erot;
    mu_err=abs(mu-mu_old)/abs(mu);
    mu_old=mu;
    
    %plot
    if mod(j,25)==0
        subplot(1,2,1)
        pcolor(abs(u).^2);shading interp;axis off;axis square;
        subplot(1,2,2)
        pcolor(heaviside(abs(u).^2-0.05*mean(abs(u(:)).^2)).*angle(u));shading interp;axis off;axis square;
        title(mu_err)
        drawnow
    end
    if any(isnan(u(:))==1)
        disp('NaNs...')
        break
    end
    j=j+1;
end