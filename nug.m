function nug=nug(T,p)
%NUG(T,P)     Kinematic viscosity of the vapor [m2/s].
%
%  Calls MUG, V.

nug = mug(T).*v(T,p);
