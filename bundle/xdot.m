function xdot=xdot(T,a)
%XDOT(T,A)  Mobility, bundle.
%  Mobility = Vapor mass flux over total mass flux.
%  XDOT = a*nu/nug.

xdot= (a.*mug(T)+(1-a).*mul(T)).*a./ (mug(T).*(a + (1-a).*rho(T).*v(T,ps(T)) ));
