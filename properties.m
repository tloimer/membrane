function properties(T)
%PROPERTIES(T)   Prints the material properties.
%
%  Calls COSTHETA, CPG, CPL, CURV, EPSILON, K, KG, KL, KAPPA, MOLM, MUG,
%  MUL, PS, R, RHO, SIG, V.

% fluid properties
[R M]=molm;
disp(['Substance: ' substance]);

% thermic properties
disp(sprintf(['molecular mass: %g kg/kmol\n' ...
  'saturation pressure: %g bar\nliquid density: %g kg/m3\n'...
  'vapor specific volume: %g m3/kg'],M,ps(T)./1e5,rho(T),v(T,ps(T))));

% caloric properties
disp(sprintf(['liquid heat capacity cp: %g J/kgK\n'...
  'vapor heat capacity cp: %g J/kgK\nenthalpy of vaporization: %g kJ/kg'],...
  cpl(T),cpg(T,ps(T)),r(T)*1e-3));

% transport properties
disp(sprintf(['liquid viscosity: %g mPas\nvapor viscosity: %g muPas\n'...
  'liquid thermal conductivity: %g W/mK\n'...
  'vapor thermal conductivity: %g W/mK\nsurface tension: %g mN/m'],...
  mul(T).*1e3,mug(T)*1e6,kl(T),kg(T),sig(T).*1e3));

% membrane properties
disp(sprintf('\nwetting angle: %g°, cos(theta)=%g',...
  acos(costheta)*180/pi,costheta));
if costheta/sqrt(1.5*kappa/epsilon)==curv
  por = 'cylindrical pores';
else
  por = 'slit pores';
end
disp(sprintf([por ', capillary pressure: %g bar'],curv.*sig(T)/1e5));
disp(sprintf(['void fraction: %g\n'...
  'thermal conductivity of vapor-filled membrane: %g W/mK'],...
  epsilon,k(T,1)))
