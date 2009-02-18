% Data:
%[substancename datamemname exp_id T1 p1 p2 Vflow Troom] = deal(data{:});
% [T1, Troom] = C; [p1, p2] = bar; [Vflow] = ml/min; 
% Now cycle through the membranes and display nitrogen data:
sname = 'nitrogen';
lendata = size(p1,1);
s = substance(sname); f = fmodel('homogeneous');
%topology = 'porousround';
topology = 'tube';
tau = 1; beta = 8.1;
plots = 0;
dothis = true;
for m = 1:size(memname,1)
  %dothis = ...
  %    input(sprintf('Display membrane %s?  1(yes)/0 [1]: ',memname{m}));
  if isempty(dothis), dothis = true; end
  if ~dothis, continue; end

  % Collect the nitrogen data for this membrane
  indices = []; % Preallocation is necessary.
  for k = 1:lendata
    if strcmp(memname{m},datamemname{k}) && strcmp(sname,substancename{k})
      indices = [indices k];
    end
  end

  if isempty(indices)
    disp(sprintf('No nitrogen data found for membrane %s.',memname{m}));
    continue;
  end

  % Assign pred [-] and pmean [Pa] , we probably need it anyway.
  lenfound = size(indices,2);
  pred = p1(indices);
  pmean = (pred + p2(indices)) / 2;
  psat = T1(indices); % allocate psat plus reuse as  T1 [K] below:
  for k = 1:lenfound
    psat(k) = s.ps(psat(k)); % ps is not vectorizable
  end
  if ~any(~isfinite(psat)) % if no infinite value is encountered in psat
    pred = pred/psat;
    xlabelstr = 'p_{red} [-]';
  else
    pred = pred*1e-5;
    xlabelstr = 'p_1 [bar]';
  end
  % else, pred stays p1 [bar].

  % Sort after pred
  % or comment out to sort after pmean.
   pred = pmean*1e-5; xlabelstr = 'p_{mean} [bar]';
  [pred sorted] = sort(pred);
  pmean = pmean(sorted);
  sorted = indices(sorted);  % Reused sorted, pmean.

  Tr = Troom(sorted);
  % Set room temperature, if a value is missing
  indices = isnan(Troom(sorted)); % Reuse indices!
  if any(indices) % Some Troom-values are NaN
    if ~any(~indices) % All values are NaN!
      %disp('Room temperature set to 25 C.');
      disp(sprintf('Membrane %s: min(T1) = %0.2f C, max(T1) = %0.2f C.',...
        memname{m},min(T1(sorted))-273.15,max(T1(sorted))-273.15));
      %Tr = input(sprintf(['Set room temperature [C] for membrane %s,'...
      %	' exp. %s: '],memname{m},exp_id{sorted(1)}));
      Tr = mean(T1(sorted));
      disp(sprintf('Room temperature set to %0.2f C.',Tr-273.15));
      Tr = ones(lenfound,1) * (Tr);%+273.15);
    else
      % Calculate mean
      Tr(find(indices)) = mean(Tr(find(~indices)))
    end
  end

%  disp(sprintf('%s',memname{m}));
%  disp(' #   p1     p2     pmean  T1     mu(T1)     nu(T1,pmean)');

  % And calculate for each data point
  % Because del p/L = massflux*nu/kappa, the normalized permeability should
  % ideally be one:
  %   massflux * L*nu/(del p*kappa) = 1.
  % Commonly, one sets nu = mu*R*T/pmean and plots the permeance, expecting it
  % to be linear in pmean:
  % massflux*L/(del p) = pmean * kappa/(mu*R*T)
  mem = membrane(poredia(m),eps(m),1.38,topology,1,0,L(m));
  flow_flux = area(m)*60*1e6; % mass flow [mg/min] over mass flux [kg/m2s]
  L_kappa = L(m)/mem.kappa;
  % delete old calc1, preallocate array of correct length
  perm = pred; nu = pred; kn_nu = pred;
  %data = pred; calc1 = pred; shouldbe1 = pred; calc2 = pred;
  for k = 1:lenfound
    ind = sorted(k);
    delp = p1(ind) - p2(ind);  vol = s.v(Tr(k),101300);
    %one_mflux = s.nug(T1(ind),pmean(k)) * L_kappa / delp;
    %disp(sprintf('p1 = %0.3f, p2 = %0.3f bar: v = %0.2f m3/kg.',...
    %  p1(ind)*1e-5,p2(ind)*1e-5,vol));
    %shouldbe1(k) = Vflow(ind) * one_mflux / (vol*flow_flux);
    %permeance(k) = Vflow(ind) * L(m) / (s.v(Tr(k),101300)*((p1(ind)-p2(ind)));
    perm(k) = Vflow(ind) * L(m) / (flow_flux*vol*delp);
    nu(k) = s.nug(T1(ind),pmean(k));
    kn_nu(k) = 3*sqrt(pi/(8*s.R)) / (sqrt(T1(ind))*mem.dia);
    %mflux = mnum(T1(ind),p1(ind),p2(ind),0,s,mem,f);
    %calc1(k) = mflux*one_mflux;
%    calc2(k) = mnum(T1(ind),p1(ind),p2(ind),0,s,mem,f)*L(m)/delp;
%    disp(sprintf(' %2u: %0.3f  %0.3f  %0.3f  %2.2f  %0.4g  %0.4g',...
%      k,p1(ind)*1e-5,p2(ind)*1e-5,pmean(k)*1e-5,T1(ind)-273.15,...
%      s.mug(T1(ind)),s.nug(T1(ind),pmean(k))));
  end

% Calculate tau and beta by least square fit, scale mflux*L/(p1-p2) to 1.
c1c2 = lscov([1./(nu.*perm) kn_nu./perm],ones(lenfound));
%c1c2 = lsqnonneg([1./(nu.*perm) kn_nu./perm],ones(lenfound));
tau = mem.kappa/c1c2(1);
beta = c1c2(2)/c1c2(1);
disp(sprintf('Membrane %s: tau = %0.4f, beta = %0.4f.',memname{m},tau,beta));

end
