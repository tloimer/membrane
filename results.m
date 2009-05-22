% readdata füllt folgende Arrays:
% [memname,poredia,L,memdia,eps,model,tau,beta] = deal(data{:});
% [substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc]...
% Temperaturleitfähigkeit hard-coded = 1.38 W/m2K.

% Arrays allokieren:
kap_kc = zeros(size(p1));
[kappac psat1 pK1 T2calc mexp mgas miso m90 mcalc] = deal(kap_kc);
% kap_kc = kappa/kappa_c;
% kappaK = kappaK(kappac), nicht genau kappaK! Schreit bei mehr als 1% Differenz
% mexp = experimentell gemessener Massenfluss,
% mgas = Massenfluss, falls der Dampf ein idealas Gas wäre,
% miso = isotherm, m90 = Kontaktwinkel von 90° (WS),
% mcalc = 0° Kontaktwinkel,
% mjt0 = (T1 = T2) aber nicht isotherm (Joule-Thomson Koeffizient = 0);
% unphysikalisch, wegen q2 = 0, nicht gerechnet.

% Um die Raumtemperatur zuzuschreiben, die manchmal fehlt, werden die Daten nach
% der Reihe durchgegangen, Daten mit gleicher Membrane und gleicher Substanz
% zusammengefasst.

f = fmodel('homogeneous');
Tnorm = 25 + 273.15;
first = 1;
lendata = size(p1,1);
while first < lendata
  mname = datamemname{first};
  sname = substancename{first};
  next = first + 1; % first muss nicht mit sich selbst verglichen werden
  while next <= lendata && ...
    strcmp(mname,datamemname{next}) && strcmp(sname,substancename{next})
    next = next + 1;
  end
  last = next - 1;
  % Nun haben wir zwischen [first:last] gleiche Substanz und gleiche Membran.
  % Raumtemperatur berechnen, falls nötig
  noTr = isnan(Troom(first:last));
  if any(noTr)
    % mind. 1 Datenpunkt hat keine Raumtemperatur
    if all(noTr)
      % kein einziger Punkt hat eine Raumtemperatur; Wir setzen 25°C;
      Troom(first:last) = Tnorm;
    else
      % einzelne Punkte ohne Raumtemperatur wird der Durchschnitt des restlichen
      % Experimentes zugewiesen.
      Troom(first-1+find(noTr)) = mean(Troom(first-1+find(~noTr)));
    end
  end
  % Ende Raumtemperatur berechnen / zuweisen.

  % Der Spass beginnt. 
  s = substance(sname);
  ismem = strcmp(mname,memname); % ein logisches Array zur Indizierung;
  % doppelt vorkommende Membranen sollten einen Fehler geben.
  if sum(ismem) > 1, error(['Mehrere gleiche Membranen ' mname '!']); end
  mem = membrane(poredia(ismem),eps(ismem),1.38,model{ismem},...
    tau(ismem),beta(ismem),L(ismem));
  flow_flux = area(ismem)*60*1e6; % mass flow [mg/min] over mass flux [kg/m2s]
  betakn_nu = mem.beta * 3*sqrt(pi/(8*s.R))/mem.dia;
  kap_L = mem.kappa/mem.L;
  for i = first:last
    % gemessener Massenstrom
    mexp(i) = Vflow(i)/(flow_flux*s.v(Troom(i),101300));
    % berechnet, numerisch, theta = 90°.
    try
      [m90(i) fl] = mnum(T1(i),p1(i),p2(i),90,s,mem,f);
      fl90{i} = fl;
      % f_handles müssen ohnehin die Funktionen finden; Vars. verkleinern.
      fl90{i}.info.substance = sname; fl90{i}.info.membrane = mname;
      fl90{i}.info.fmodel = fl90{i}.info.fmodel.name;
    catch
      disp(sprintf('Fehler bei Punkt %u: m90',i));
    end
    % berechnet, numerisch, theta = 0°.
    try
      [mcalc(i) fl] = mnum(T1(i),p1(i),p2(i),0,s,mem,f); 
      T2calc(i) = fl.sol.T2;
      fl0{i} = fl;
      % f_handles müssen ohnehin die Funktionen finden; Vars. verkleinern.
      fl0{i}.info.substance = sname; fl0{i}.info.membrane = mname;
      fl0{i}.info.fmodel = fl0{i}.info.fmodel.name;
    catch
      disp(sprintf('Fehler bei Punkt %u: mcalc',i));
    end
    % Gasströmung, isotherm
    corrKn = betakn_nu/sqrt(T1(i));
    mgas(i) = kap_L*(p1(i)-p2(i))*( 1/s.nug(T1(i),(p1(i)+p2(i))/2) + corrKn );
    % Kondensation ist nur bei Dämpfen möglich
    if ~strcmp(sname,'nitrogen')
      [psat1(i) tmp] = s.ps(T1(i)); %tmp ist dps/dT
      kappac(i) = s.nul(T1(i))*f.kmliq(mem.epsilon,mem.km,s.kl(T1(i)))...
        /(tmp*s.hvap(T1(i)));
      kap_kc(i) = mem.kappa/kappac(i);
      % kappak, Ccc weiss ich nicht - in mnum versteckt
      % isotherm, eventuell mit Kondensation
      pcap = mem.fcurv(1)*s.sigma(T1(i));
      pK1(i) = psat1(i)*exp(-pcap/(s.R*s.rho(T1(i))*T1(i)));
      if p1(i) >= pK1(i)
        pv_nu = (pK1(i)-p2(i))*( 1/s.nug(T1(i),(pK1(i)+p2(i))/2) + corrKn );
        pl_nu = (p1(i)-pK1(i)+pcap*(p1(i)-pK1(i))/(psat1(i)-pK1(i)))...
	  /s.nul(T1(i));
        miso(i) = kap_L*(pv_nu+pl_nu);
      else
        miso(i) = mgas(i);
      end
    else % nitrogen
      miso(i) = mgas(i);
    end %end if ~nitrogen
  end
  first = next;
end

file = 'results0902';
dlmwrite([file '.tsv'],...
  [kap_kc kappac psat1 pK1 T2calc mexp mgas miso mcalc m90],'\t');
save(file,'kap_kc','kappac','psat1','pK1','T2calc','mexp','mgas','m90',...
  'miso','mcalc','fl0','fl90');
