%% Test of expression for mass flux in mnum.m, 1.4.2009.
function testmnum()

%% Determine the type of flow according to the flow map, Fig. 1.
%% Wetting cases.
% Missing cases (cap. condens. with liq. film and non-wetting with liq. flow
% were hardly reachable.

f = fmodel('plug');
s = substance('isobutane'); 

p1=3.1e5; Troom=21;
T1 = Troom + 273.15;

disp(sprintf('\n2ph - vapor'));
mdia = 180e-9; p2 = 1e5; theta = 0;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);

disp(sprintf('\nliq.flow - vapor'));
mdia = 70e-9;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);

disp(sprintf('\nliq.film - liq.flow - vapor'));
mdia = 10e-9; theta = 60;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);

disp(sprintf('\nliq.flow'));
p1=3.1038e5; p2 = 3e5;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);

disp(sprintf('\nliq.film - liq.flow: missing'));

disp(sprintf('\n2ph - liq.flow'));
p1 = 3.103e5; mdia = 200e-9; p2 = 3.0684e5; theta = 0;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);
p1 = 3.1e5;

disp(sprintf('\nNon - wetting'));

disp(sprintf('liq.film - vapor'));
mdia = 20e-9; p2 = 1e5; theta = 120;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);

disp(sprintf('\nliq.film - liq.flow - vapor'));
mdia = 72e-9; p2 = 1e4; theta = 100;
s = substance('ethanol'); 
T1 = 376; Troom = T1-273.15;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);

disp(sprintf('\nliq.film - 2ph - vapor'));
mdia = 120e-9; p2 = 3e4;
analysis(mdia,T1,p1,p2,Troom,theta,s,f);


%%% SUBFUNCTIONS %%%

function analysis(mdia,T1,p1,p2,Troom,theta,s,f)
meps = 0.38; mkm=1.38; mtopology = 'porousround';
mtau = 1; mbeta = 8.1; mL = 1e-3;

mem = membrane(mdia,meps,mkm,mtopology,mtau,mbeta,mL);

% lineare Berechnung und Ausgabe
flin = mlinear(p1,p2,Troom,theta,s,mem,f);
disp(sprintf('MLINEAR: %u Strömungszustände, %s',flin.sol.len,...
  [flin.flow(1:flin.sol.len).color]));
disp(sprintf('  p1 = %.2f bar, p2 = %5.3f bar, T1 = %4.1f °C.', ...
  flin.info.p0/1e5, flin.sol.pe/1e5, flin.info.T0-273.15));
disp(sprintf('  massflux (TL, 2007): %g kg/m2s',flin.lin.m));
merror = 1-flin.info.m/flin.lin.m;
if merror > 1e-8
disp(sprintf(['  Difference between mass flux mlin from formula copied from'...
   ' TL (2007)\n  and mass flux m calculated in file, 1-m/mlin = %g'],merror));
end
disp(sprintf('  kappa/kappaK = %f, CC = %g', ...
  flin.info.kap/flin.info.kapc,flin.info.C));

%numerische Berechnung
try
[m fnum] = mnum(T1,p1,p2,theta,s,mem,f,'m',flin.lin.m);
disp(sprintf('MNUM: %u Strömungszustände, %s = %s',-fnum.sol.len,...
  [fnum.flow(1:-fnum.sol.len).color], fnum.sol.states));
disp(sprintf('  T1 = %5.2f °C, massflux: %g kg/m2s',...
  fnum.info.T1-273.15,fnum.sol.m));
catch ME
for i = 1:length(ME.stack)
disp(sprintf('Fehler in %s, Zeile %u.',ME.stack(i).name,ME.stack(i).line));
end
disp(ME.identifier); disp(ME.message);
end
