function [ps dps] = ps(T)
%PS(T)      Vapor pressure [Pa].
%
%  PS(T) returns the vapor pressure.
%
%  [P DP] = PS(T) returns the vapor pressure and the derivative of the vapor
%  pressure with respect to T.
%
%  Calls ANTOINEEQ.
%  See also TS.

A = Antoineeq(T);
[Aa Ab Ac T0 Tc n E F] = deal(A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8));

ps = 10.^(Aa-Ab./(Ac+T));

% If A(4)(=T0) == 0 ANTOINEEQ returns the classical Antoine equation, otherwise
% coefficients for the extended Antoine eq. are returned. See Landolt-Börnstein,
% New Series, Group IV: Physical Chemistry.  Vapor Pressure of Chemicals, vol.
% 20A: J. Dykyj, J. Svoboda, R.C. Wilhoit, M. Frenkel, K.R. Hall (1999).
if T0 ~= 0
  chi = (T-T0)/Tc;
  ps = ps.*10.^(0.43429.*chi.^n + E.*chi.^8 + F.*chi.^12);
end

% The derivative. Mathematica says:
% In[1]:= D[10.^(Aa-Ab/(Ac+T)+0.43429*chi^n+E*chi^8+F*chi^12)/.
%   chi->(T-T0)/Tc,T]/.T0->T-chi*Tc
% Out[2]//InputForm= 
% 2.302585092994046*10.^(Aa + 0.43429*chi^n + chi^8*E + chi^12*F - Ab/(Ac + T))*
%  (Ab/(Ac + T)^2 + (8*chi^7*E)/Tc + (12*chi^11*F)/Tc + 
%   (0.43429*chi^(-1 + n)*n)/Tc)

if nargout==2
  %ln10 = log(10);
  ln10 = 2.302585092994046;
  dps = ln10.*ps.*Ab./(Ac + T).^2;
  if T0 ~= 0
    dps = dps+ln10*ps.*(8*chi.^7.*E+12*chi.^11.*F+0.43429.*n.*chi.^(n-1))./Tc;
  end
end
