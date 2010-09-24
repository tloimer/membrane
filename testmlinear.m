%% Test of expression for mass flux in mlinear.m, 4.4.2008.
function testmlinear()

%% Determine the type of flow according to the flow map, Fig. 1.
%% Wetting cases.
% Missing cases (cap. condens. with liq. film and non-wetting with liq. flow
% were hardly reachable.

%% 2ph - vapor
fl = mlinear(3.1e5,1e5,23,150e-9,0,'isobutane');
analysis(fl);

%% liq.flow - vapor
fl = mlinear(3.1e5,1e5,23,70e-9,0,'isobutane');
analysis(fl);

%% liq.film - liq.flow - vapor
fl = mlinear(3.1e5,1e5,23,10e-9,60,'isobutane');
analysis(fl);

%% liq.flow
fl = mlinear(3.1e5,3e5,23,10e-9,60,'isobutane');
analysis(fl);

%% 2ph - liq.flow
fl = mlinear(3.1e5,3.0525e5,23,150e-9,0,'isobutane');
analysis(fl);

%% Non - wetting

%% liq.film - vapor
fl = mlinear(3.1e5,1e5,23,20e-9,120,'isobutane');
analysis(fl);

%% liq.film - 2ph - vapor
fl = mlinear(3.1e5,1e4,23,72e-9,100,'ethanol');
analysis(fl);

function analysis(fl)
disp(sprintf('  p1 = %4.2f bar, T1 = %4.1f °C.', ...
  fl.info.p0/1e5, fl.info.T0-273.15));
disp(sprintf(['  Difference between mass flux mlin from formula copied from'...
   'TL\n  and mass flux m calculated in file, 1-m/mlin = %g'], ...
   1-fl.info.m/fl.lin.m));
disp(sprintf('  kappa/kappaK = %f, CC = %g', ...
  fl.info.kap/fl.info.kapc,fl.info.C));
