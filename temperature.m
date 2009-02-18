%     Temperaturvergleich:
% Kalibration:
%   T1 [C,K] = Troom [C,K] - T1 [muV] / 41
%   T2 [C,K] = Troom [C,K] - T2 [muV] / 41
%   T1 -T2 = (T1-T2) [muV]/(-41)
% Es entspricht  T1tc = T1 [muV], T2tc = T2 [muV], T12tc = (T1-T2) [muV].
%
% Data:
% [substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc] ...
%   = deal(data{:});
% [T1, Troom] = C; [p1, p2] = bar; [Vflow] = ml/min; 

%____Aufgabe 1: Vergleiche T1tc mit T1 (geschätzt)
% if T1Tc exists & Troom exists T1est = T1, T1meas = Troom - T1tc/41
k = isfinite(Troom) & isfinite(T1tc);
T1est = T1(k);
T1meas = Troom(k) - T1tc(k)./41;
disp('___Vergleich T1 mit T1 durch thermocouple gemessen');
disp(sprintf('%u Datenpunkte, davon haben %u Troom und T1 [mKV] gesetzt.',...
  size(T1,1),size(T1est,1)));
trange = [10 30]; xlim(trange); ylim(trange);
plot(T1est-273.15,T1meas-273.15,'k.',trange,trange,'k-');
xlabel('T_1 estimated, T_{room} or from Ts(p1)')
ylabel('T_1 by thermocouple')

%___Aufgabe 2: Vergleiche T2 - T1 [uV] (Thermocouple-Messung) mit Delta T über
% Joule-Thomson Effekt.
indT1aT2 = isfinite(T1tc) & isfinite(T2tc);
indT12 = isfinite(T12tc);
disp('___Vergleich T1 - T2 mit JT*(p1-p2)___');
disp(sprintf('%u Datenpunkte mit T1 und T2 gemessen',sum(indT1aT2)));
disp(sprintf('%u Datenpunkte mit T1 - T2 gemessen',sum(indT12)));
disp(sprintf('%u Datenpunkte sind doppelt gemessen: T1, T2 und T1-T2.',...
  sum(indT1aT2 & indT12)));
T1aT2meas = (T1tc(indT1aT2)-T2tc(indT1aT2))/(-41);
T12meas = T12tc(indT12)/(-41);
% Aber substance(..) ist nicht vektorisierbar:
% Arrays für Resultat allokieren
T1aT2calc = T1aT2meas;
T12calc = T12meas;
% das ging nicht:
% T1aT2calc(k) = T1tmp(k)-substance(sname{k}).intjt(T1tmp(k),p1tmp(k),p2tmp(k));
% ??? Undefined variable "substance" or class "substance".
% Bug, oder was?
% Bauen wir also ein struct f:
f.butane = substance('butane');
f.isobutane = substance('isobutane');
f.nitrogen = substance('nitrogen');
f.ethanol = substance('ethanol');
% Benötigte Daten zuweisen
ind = indT1aT2; numT1aT2 = find(ind);
T1tmp = T1(ind); p1tmp = p1(ind); p2tmp = p2(ind); sname = {substancename{ind}};
for k = 1:sum(ind)
  T1aT2calc(k) = T1tmp(k)-f.(sname{k}).intjt(T1tmp(k),p1tmp(k),p2tmp(k));
end
ig = T1aT2meas > 0.8*T1aT2calc & T1aT2meas < 1.25*T1aT2calc;

ind = indT12; %numT12 = find(ind);
T1tmp = T1(ind); p1tmp = p1(ind); p2tmp = p2(ind); sname = {substancename{ind}};
for k = 1:sum(ind)
  T12calc(k) = T1tmp(k)-f.(sname{k}).intjt(T1tmp(k),p1tmp(k),p2tmp(k));
end

% Plot
figure;
trange = [min([T1aT2calc;T12calc]) max([T1aT2calc;T12calc])];
plot(T12calc,T12meas,'k.',T1aT2calc,T1aT2meas,'r.',trange,trange,'k-');
title('K: \Delta T gemessen (K.), rot: T1 und T2 gemessen (R.)')
xlabel('T1 - T2 [K] by JT-coefficient');
ylabel('T1 - T2 [K] measured');
disp('So, und nun xlim([-1 5]); ylim([-1 5]);');

disp(sprintf('%u .gute. Punkte',sum(ig)));
disp('Drucken: z.B. T1(numT1aT2(ig))');
figure;
plot(T1aT2calc(ig),T1aT2meas(ig),'k.',trange,trange,'k-');
