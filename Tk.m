function Tk = Tk(p)
%TK(P)      Saturation temperature at a curved interface [K].
%
%  See PS. See also DPSDT.
%
%  Calls DPKDT, KELV, PS, TS.

% an estimate, accurate to first order:
%   (pk(Ts(p))-p) / (Ts(p) - Tk) = dpkdT(Ts)

ts = Ts(p);
Tk = ts - (kelv(ts)-1).*p./dpkdT(ts);
ts = Tk-ts; % i reuse a variable!

options =optimset(optimset('fzero'),'TolX',1e-9);
Tk = fzero(@pres,Tk+ts*[-.1 .1],options,p);

function dt = pres(T,p)
dt = p - kelv(T).*ps(T);
