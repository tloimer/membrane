function nu=nu(T,a)
%NU(T,A)    Effective kinematic viscosity of the 2ph-mixture, annular [m2/s].
%
%  NU = rho^-2 *( (1-x^2)/(mul*rhol) + x^2/(mug*rhog) )^-1

rhol = rho(T);
vg = v(T,ps(T));
rho1 = a./vg+(1-a).*rhol;
rho2 = rho1.*rho1;
x1 = a./(vg.*rho1);
x2 = x1.*x1;

nu = 1./( rho2.* ( (1-x2)./(mul(T).*rho1) + x2.*vg./mug(T) ));

%%%
%p = ps(T);
%ng=nug(T,p);
%nl=nul(T);

%nu=1./( x.^2.*(1./ng+1./nl) + (1-2*x)./nl + 2*x.*(1-x)./(v(T,p).*nl) );
