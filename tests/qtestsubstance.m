function qtestsubstance
%QTESTSUBSTANCE Quickly test a few data points of each substance

% Data retrieved 2014-09-19
% Simply export data form http://webbook.nist.gov/chemistry/fluid to tab
% delimited data, columns:
% For saturation properties - temperature increments
% 1 Temperature (K),  2 Pressure (MPa),  3 Density (l, kg/m3),
% 4 Volume (l, m3/kg),  5 Internal Energy (l, kJ/kg),  6 Enthalpy (l, kJ/kg),
% 7 Entropy (l, J/g*K),  8 Cv (l, J/g*K),  9 Cp (l, J/g*K),
% 10 Sound Spd. (l, m/s),  11 Joule-Thomson (l, K/MPa),  12 Viscosity (l, Pa*s),
% 13 Therm. Cond. (l, W/m*K),  14 Surf. Tension (l, N/m), 15 Density (v, kg/m3),
% 16 Volume (v, m3/kg),  17 Internal Energy (v, kJ/kg),  18 Enthalpy (v, kJ/kg),
% 19 Entropy (v, J/g*K),  20 Cv (v, J/g*K),  21 Cp (v, J/g*K),
% 22 Sound Spd. (v, m/s),  23 Joule-Thomson (v, K/MPa),  24 Viscosity (v, Pa*s),
% 25 Therm. Cond. (v, W/m*K)

fprintf('Values NIST (matlab) unit\n\n');

try
s = substance('propane');
C = [300.00 0.99780 489.48 0.0020430 268.36 270.40 1.2421 1.6799 2.7482 ...
  707.76 0.017115 9.5302e-05 0.092868 0.0067563 21.626 0.046242 556.86 603.0 ...
  2.3508 1.5905 2.0406 214.72 21.268 8.3399e-06 0.019236];
printinfo;
catch
fprintf('Propane not implemented.\n');
end

fprintf('\n');
s = substance('isobutane');
C = [300.00 0.37000 548.32 0.0018237 262.82 263.50 1.2203 1.6903 2.4422 ...
  810.25 -0.22061 0.00014822 0.088602 0.0098917 9.6096 0.10406 551.87 590.37 ...
  2.3099 1.5781 1.8100 197.74 24.641 7.5457e-06 0.017024];
printinfo;

fprintf('\n');
s = substance('butane');
C = [300.00 0.25760 570.68 0.0017523 263.54 264.00 1.2225 1.7293 2.4512 ...
  890.88 -0.27027 0.00015563 0.10393 0.011603 6.5164 0.15346 584.05 623.58 ...
  2.4211 1.6018 1.8111 202.15 26.968 7.4364e-06 0.016784];
printinfo;

fprintf('\n');
s = substance('nitrogen');
C = [80 0.13687 793.94 0.0012595 -116.75 -116.58 2.9028 1.0691 2.0555 824.36...
  -0.32322 0.00014505 0.14020 0.0082740 6.0894 0.16422 56.622 79.099 5.3487 ...
  0.77733 1.1449 176.72 25.490 5.6413e-06 0.0077802];
printinfo;

% A nested function - can access C and s without passing it as input arguments.
function printinfo()
T = C(1);
p = C(2)*1e6;
fprintf(['%s at T = %.0f K:\n  p = %g (%g) Pa;  '...
  'rhol = %g (%g) kg/m3;  cpl = %g (%g) J/kgK\n  '...
  'mul = %g (%g) Pas;  kl = %g (%g) W/mK;  sigma = %g (%g) N/m\n  '...
  'v = %g (%g) m3/kg;  cpg = %g (%g) J/kgK\n  mug = %g (%g) Pas;  ' ...
  'kg = %g (%g) W/mK;\n'], s.name, ...
  T, p, s.ps(T), C(3),s.rho(T), C(9)*1e3,s.cpl(T), C(12),s.mul(T),...
  C(13),s.kl(T), C(14),s.sigma(T), C(16),s.v(T,p), C(21)*1e3,s.cpg(T,p),...
  C(24),s.mug(T), C(25),s.kg(T));
end

end % end qtestsubstance
