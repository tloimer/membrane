function lambda=lambda(T,a)
%LAMBDA(T,x) Mobility, plug.
%  lambda=kappag*nu/nug.

% plug:
%p=ps(T);
% lambda = a*rhog/rho = a*rhog/(a*rhog + (1-a)*rhol) = a/(a + (1-a)rhol*vg)
lambda = a./( a + (1-a)*v(T,ps(T)).*rho(T) );
