function x=x(T,a)
%X(T,A)     Mass fraction of the vapor.
%  X = A*rhog/rho, where A is the volume fraction of the vapor.

% x = a*rhog/rho = a*rhog/(a*rhog + (1-a)*rhol) = a/(a + (1-a)rhol*vg)
x = a./( a + (1-a).*v(T,ps(T)).*rho(T) );

% the reverse:
% a = x./((1-x)./(v(T,ps(T)).*rho(T))+x)
