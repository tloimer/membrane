function Ts = Ts(p)
%TS(P)      Saturation temperature [K].
%
%  TS(P) returns the saturation temperature by inverting the Antoine equation or
%  by iteratively solving the extended Antoine equation.
%
%  See PS.
%  Calls ANTOINEEQ, PS, NEWTON.

% The Antoine eq. underpredicts vapor pressure. Therefore, the temperature
% returned by the Antoine eq., Ta, is higher than the saturation temperature
% that is obtained by using the extended Antoine eq., Ts. Matlab does not have
% Newtonian iteration, therefore we program it ourselves. Otherwise, use fzero
% with the start interval [T0, Ta].
%
% Ps |  ext. Ant.    Ant. 
%    |         /   ´ 
%    | -------/--´
%    |       /:´:
%    |      /´: :
%    |  __-´  : :
%    |________________
%         T0 Ts Ta    Ts

%log(10) = 2.302585092994046;

% Get all Antoine coeffs. by using a temperature that is for sure higher than
% any T0.
A = Antoineeq(1000);
[Aa Ab Ac T0] = deal(A(1),A(2),A(3),A(4));

%ps = 10.^(Aa-Ab./(Ac+T));
Ts = Ab./(Aa-log(p)/2.302585092994046)-Ac;

if Ts>T0 
  % We are in the range of the extended Antoine eq.
  Ts = newton(@ps,Ts,p,1e-2);
end
% RES in NEWTON: With ps ~ 1e5, ps is solved to 7 digits accuracy. Since Ts is
% then accurate to RES/dps, and dps > 1e3, Ts is accurate to 1e-5. 
