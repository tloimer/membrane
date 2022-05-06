%% Simply produce figure windows and show plots
% See nice_plots.m for export to pdf-files.

if ~exist('substance.m')
    %addpath('~/matlab/program');
    addpath('/j/a/program');
end

data = readmatrix('AAO_40nm_isobutane_T.csv', 'NumHeaderLines', 1);
% Columns:
% 1 - Feed pressure, bar;  2 - Permeate pressure, bar;
% 3 - Saturation pressure, bar;  4 - Flux, ml/min;
% 5 - Permeance, m3/(m2 bar h);  6 - Temperature, °C;
desc = table2struct(readtable('AAO_40nm_isobutane_desc.csv'));
% desc.file, .year, .month, .day,
%    .material, .substance, .poredia, .porosity, .area, .thickness

% thermal conductivity
if (strcmp(lower(desc.material), 'aao'))
    km = 1.22;
end
f = fmodel('parallel');
s = substance(desc.substance);

%% Auxiliary functions

% Convert volume flow rate Q to permeance.
% Arguments: Q (ml/min), p1 (Pa), p2 (Pa), area (m2).
% Result:    permeance (m3/m2 bar h)
Qtopermeance = @(Q,p1,p2,area)  6 * Q ./ (area * (p1 - p2));

% Convert mass flux m to permeance.
% Assume, permeance is cast to pressure at atmospheric conditions!
% (= s.v(20 °C, 1 bar))
% Arguments: m (kg/m2 s), p1 (Pa), p2 (Pa).
% Returns:   permeance (m3/m2 bar h)
mtopermeance = @(m,p1,p2)  3600 * m * s.v(273.15+20, 1e5) * 1e5 / (p1 - p2);


% Set membrane properties
% membrane(pore dia., epsilon, thermal cond., tname, tau, beta, thickness)
mem = membrane(desc.poredia, desc.porosity, km, 'tube', 1, 9.0541,...
            desc.thickness);
% mstackstruct(theta, {{layer1 layer2} {second membrane}}, fmodel)
ms = mstackstruct(0, {{mem}}, f);

len = size(data,1);
madi = zeros(len,1);
Tadi = madi;
T1 = madi;
ma = madi;
Ta = madi;
mb = madi;
Tb = madi;
md = madi;
Td = madi;
me = madi;
Te = madi;
mf = madi;
Tf = madi;

% Remember, data.data(data_points,columns):
% (:,1)  (:,2)  (:,3)    (:,4)      (:,5)                 (:6)
% p1/bar p2/bar psat/bar Q/(ml/min) Permeance m3/m2 bar h T_2/°C
for i = 1 : size(data,1)
    p1 = data(i,1) * 1e5;
    p2 = data(i,2) * 1e5;
    T1(i) = s.Ts(data(i,3) * 1e5);

    [m, mr] = mnumadiabat(T1(i), p1, p2, s, ms);
    Tadi(i) = mr.T2;
    madi(i) = mtopermeance(m, p1, p2);

    [m, mr] = mT1eqT2(T1(i), p1, p2, s, ms);
    Td(i) = mr.T2;
    md(i) = mtopermeance(m, p1, p2);

    [m, mr] = mnumTdiff(0, T1(i), p1, p2, s, ms);
    Ta(i) = mr.T2;
    ma(i) = mtopermeance(m, p1, p2);

try
    [m, mr] = mnumTdiff(1, T1(i), p1, p2, s, ms);
    Tb(i) = mr.T2;
    mb(i) = mtopermeance(m, p1, p2);
catch
    Tb(i) = Td(i);
    mb(i) = md(i);
fprintf("mnumTdiff(1,..), error: i = %d\n", i);
end

    [m, mr] = mnumtest(0.3, T1(i), p1, p2, s, ms);
    Te(i) = mr.T2;
    me(i) = mtopermeance(m, p1, p2);

    [m, mr] = mnumTdiff(0.3, T1(i), p1, p2, s, ms);
    Tf(i) = mr.T2;
    mf(i) = mtopermeance(m, p1, p2);
end

pr = data(:,1)./data(:,3);
% plot temperature difference versus  p1/psat
%figure()
plot(pr,T1 - Tadi,'xb', pr,T1 - Ta,'ob', pr,T1 - Td,'+r', pr,T1 - Tb,'or',...
	pr,T1-Te,'xk', pr,T1-Tf,'ok');
legend({'adiabat', 'fac = 0', 'diabat', 'fac = 1',...
		'mnumtest(0.3)','mnumTdiff(0.3)'}, 'Location', 'northwest');
legend('boxoff');
xlabel('p_1/p_{sat}');
ylabel('T_1 - T_2 (°C)');
ylim([-0.1 7])

figure()
plot(pr,madi,'xb', pr,ma,'ob', pr,md,'+r', pr,mb,'or', pr,me,'xk', pr,mf,'ok');
legend({'adiabat', 'fac = 0', 'diabat', 'fac = 1', ...
		'mnumtest(0.3)','mnumTdiff(0.3)'}, 'Location', 'northwest');
legend('boxoff');
xlabel('p_1/p_{sat}');
ylabel('permeance  (m^3/m^2 bar h)');
ylim([0 150]);
