% readdata füllt folgende Arrays:
% [memname,poredia,L,memdia,eps,model,tau,beta] = deal(data{:});
% [substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc]...
% Temperaturleitfähigkeit hard-coded = 1.38 W/m2K.

% Arrays allokieren:
kappa = zeros(size(p1));
ishphil = logical(kappa);
%[kappac psat1 pK1 T2calc mexp mgas miso m90 mcalc] = deal(kap_kc);

% Die Daten werden der Reihe nach durchgegangen, Daten mit
% gleicher Membrane zusammengefasst. %% und gleicher Substanz zusammengefasst.

f = fmodel('homogeneous');
%Tnorm = 25 + 273.15;
first = 1;
lendata = size(p1,1);
while first < lendata
  mname = datamemname{first};
  % gleiche Substanz
  % sname = substancename{first};
  next = first + 1; % first muss nicht mit sich selbst verglichen werden
  while next <= lendata && ...      % gleiche Substanz __
    strcmp(mname,datamemname{next}) % && strcmp(sname,substancename{next})
    next = next + 1;
  end
  last = next - 1;

  % Der Spass beginnt. 
  % s = substance(sname);
  ismem = strcmp(mname,memname); % ein logisches Array zur Indizierung;
  % doppelt vorkommende Membranen sollten einen Fehler geben.
  if sum(ismem) > 1, error(['Mehrere gleiche Membranen ' mname '!']); end
  mem = membrane(poredia(ismem),eps(ismem),1.38,model{ismem},...
    tau(ismem),beta(ismem),L(ismem));
  kappa(first:last) = mem.kappa;
  % hphil oder hphob wurde in readdata nicht eingelesen; hier wird aus dem
  % Membrannamen auf die Benetzungseigenschaft geschlossen
  ishphil(first:last) = isempty(strfind(mname,'HP')) &  ...
    isempty(strfind(mname,'U'));

%  disp(sprintf('%8s: kappa = %9.3e m2,  poredia = %3.0f nm,  eps = %4.2f',...
%    mname,mem.kappa,mem.dia*1e9,mem.epsilon));

  first = next;
end
