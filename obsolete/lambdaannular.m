function lambda=lambda(T,x)
%LAMBDA(T,x) Mobility, annular.
%  lambda=kappag*nu/nug.

lambda = 1-((1-x) + x.*(x-1)).*nu(T,x)./nul(T);
