function lambda=lambda(T,x)
%LAMBDA(T,x) Mobility, bundle.
%  lambda=kappag*nu/nug.

lambda = x.*nul(T)./(x.*nul(T)+(1-x).*nug(T,ps(T)));
