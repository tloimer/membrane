function calctaubeta(mem,data)
%CALCTAUBETA Calculate tau and beta.
%  CALCTAUBETA(MEM,DATA) computes the tortuosity tau and the molecular flow
%  correction factor beta for the homogeneous membrane MEM, given the
%  permeance data DATA. The membrane MEM must contain a field 'area'. The
%  results for tau and beta are displayed. The permeance data as well as
%  the 95% confidence intervals for the regression line and for the
%  permeance data are plotted.
%
%  The regression line is constructed for
%    permeance/(kn_nu*mem.kappa) = beta/tau + (1/tau)/(kn_nu*nu)
%  with permeance = massflux * L / (p1 - p2), nu_app = nu / (1 + beta*Kn),
%  kn_nu = Kn/nu. Hence, the regression line and especially the confidence
%  intervals are valid for beta/tau and 1/tau, not for tau and beta.
%  (Usually, the differences are small).
%
%  See also MEMBRANE, READDATA.

% Copied and reworked from 11jms/matlab3/membranes11.m. Membranes11.m,
% in turn, derived from 11jms/matlab/calctaubeta.m.
% History of documentation below.

% CALCTAUBETA Berechne tau und beta durch lineare Regression an N2-Daten.
% Mit Darcy'schem Gesetz, massflux = (kappa/nu_app) (del p/L), del p = p1 - p2;
% und nu_app = nu/(1 + beta Kn), Kn = lambda/d = 3nu\sqrt(pi/(8RT))/d,
% siehe Gl. (3), Loimer (2007), sowie perm = massflux*L/(delp),
%   perm = mem.kappa/(tau nu) + (beta*mem.kappa/tau)*kn_nu.
% Zuerst wurde für die Bestimmung von tau und beta lscov verwendet
%   1 = (mem.kappa/tau)/(perm*nu) + (beta*mem.kappa/tau)*(kn_nu/perm),
% mit c1c2(1) = mem.kappa/tau, c1c2(2) = beta*c1c2(1).
% Einfacher ist aber
%  1/nu = tau*(perm/mem.kappa) - beta*kn_nu.
%
% Lineare Regression verlangt besser
%  perm/(kn_nu*mem.kappa) = (1/tau)/(kn_nu*nu) + (beta/tau)
% Mit der Konstante mem.kappa, fast einer Konstante kn_nu und den bei jeder
% Messung veränderlichen perm und nu, mit nu proportional zu 1/pmean. Damit
% entspricht in der üblichen Statistik-Notation
%  hat alpha (oder hat beta_0) = (beta/tau),
%  hat beta (oder hat beta_1) = (1/tau).

% Literature: A.V. Metcalfe, Statistics in Engineering, 1994.
% See pages 209ff for linear regression, pages 163ff for the
% Student-distribution

len = size(data.substance,1);

% Get the substance. Check, if only one substance is given.
% data is a struct arranged in a column-vector of values
s = substance(data.substance{1});
if ~all(strcmp(s.name,{data.substance{2:len}}))
  error('All data must be measured with the same substance.');
end

if mem.tau ~= 1
  new = membrane(mem.dia,mem.epsilon,mem.km,mem.tname,1,8.1,mem.L);
  new.area = mem.area;
  mem = new;
end

% Sort after (p1 + p2) / 2
pmean = 0.5 * (data.p1 + data.p2);
[pmean ind] = sort(pmean,1);

% Set the mean temperature within the membrane.
if isfield(data,'T1')
 if isfield(data,'T2')
   Tmean = 0.5 * (data.T1(ind) + data.T2(ind));
 else
   Tmean = data.T1(ind);
 end
else
 if isfield(data,'Troom') && all(isfinite(data.Troom))
   Tmean = data.Troom(ind);
 else
   error('Neither upstream nor room temperature given.');
 end
end

% Set the room temperature
if isfield(data,'Troom')
  Tr = data.Troom(ind);
  if any(~isfinite(Tr))
    Tinf = find(~isfinite(Tr));
    Tr(Tinf) = Tmean(Tinf);
  end
else
  Tr = Tmean;
end

% s.v is vectorizable, with column vectors Troom and proom
permeance = data.volume(ind) .* mem.L ./ (data.duration(ind) .* mem.area ...
		.* s.v(Tr,data.proom(ind)) ...
		.* (data.p1(ind) - data.p2(ind)));
nu = s.nug(Tmean,pmean);
kn_nu = 3*sqrt(pi./(8*s.R*Tmean)) / mem.dia;

% Regression:
%  perm/(kn_nu*mem.kappa) = beta/tau + (1/tau)/(kn_nu*nu)
xi = 1./(nu.*kn_nu);
[b,bint,r,rint,stats] = ...
  regress(permeance./(mem.kappa.*kn_nu), [ones(len,1) xi]);
% [b,bint,r,rint,stats] = regress(y,[x1 x2]) computes
%   y = b(1)*x1 + b(2)*x2,
%   bint - confidence intervals vor b,
%   r - errors, i.e., deviation from the regress line,
%   rint - 95% confidence intervals for the data, centered at r_i, should
%          therefore contain 0, otherwise this point is an outlier,
%   stats - stats(1) - R^2 statistics, stat(2) - F statistic,
%           stat(3) - p-value,
%   stats(4) - estimate of std. dev^2 = sum r^2 / (n - 2),
%              n - number of observations
%  The confidence interval for the regress line, at point x1 = 1 and x2 = xi,
%  is given by
%    tinv(0.975,n-2) * std.dev * sqrt(1/n + (xi-mean(x))^2 / sum(xj-mean(x))^2).
%  For the confidence interval for the data, at x1 = 1 and x2 = xi, 1 must be
%  added to the expression under the square root,
%    tinv(0.975,n-2) * std.dev
%          * sqrt(1 + 1/n + (xi-mean(x))^2 / sum(xj-mean(x))^2).

% To compute the mass flux, 1/tau and beta/tau are used, hence regress on those.
fprintf('regression for 1/tau, beta/tau, 95%% confidence interval\n');
tr = 1/b(2); trlo = 1/bint(2,2); trhi = 1/bint(2,1);
br = b(1)*tr; brlo = bint(1,1)*trlo; brhi = bint(1,2)*trhi;
fprintf(' tau = %0.4f + %0.4f - %0.4f    (1/tau = %0.4f +/- %0.4f)\n',...
  tr,trhi-tr,tr-trlo, b(2),bint(2,2)-b(2));
fprintf(' beta = %0.4f + %0.4f - %0.4f (beta/tau = %0.4f +/- %0.4f)\n',...
  br,brhi-br,br-brlo, b(1),bint(1,2)-b(1));

% Statistics for the output
rstd = sqrt(stats(4));	% stats(4) = sum of errors^2 / (n - 2)
rsumsq = var(xi,1)*len;	% var(X,1) = sum (xi - mean(x))^2 / n
rxmean = mean(xi);
fprintf(['  statistical parameters\n' ...
  '  std. dev. %.4f,  samples %d,  mean %.4g,  sum sq. err %.4g\n'], ...
  rstd, len, rxmean, rsumsq);

% Better use correct statistics.
% This would be, I believe, the values of tau and beta, if not used in the
% combination 1/tau, beta/tau. Just print for informational purposes.
varb2 = stats(4)/rsumsq;
taur = tr*(1+varb2*tr^2); betar = br*(1+2*varb2*tr^2);
fprintf('regression for tau and beta: tau = %0.4f, beta = %0.4f.\n',...
  taur,betar);

% Plot the result.
memtb = membrane(mem.dia,mem.epsilon,mem.km,mem.tname,tr,br,mem.L);
regressplot(permeance, pmean, Tmean, s, memtb, rstd, len, rxmean, rsumsq);


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function regressplot(perm, pmean, Tmean, s, mem, rstd, len, rxmean, rsumsq)
%REGRESSPLOT Visualize the result of a linear regression.

% pmean is already sorted in ascending order.
gran = 0.2e5;		% granularity for display: 0.2 bar
t = pmean(len) - pmean(1);
pmin = floor((pmean(1) - 0.1*t) / gran) * gran;
pmax = ceil((pmean(len) + 0.1*t) / gran) * gran;

% Calculate expected permeance data for the mean of all measured temperatures.
T = mean(Tmean);
cp = [pmin pmin + [1:10] * 0.1 * (pmax-pmin)]';
ckn_nu = 3*sqrt(pi/(8*s.R*T)) / mem.dia;
cmu = s.mug(T);
cnu = s.mug(T) * s.v(T,cp);
cperm = mem.kappa * (1./cnu + mem.beta*ckn_nu);

% Compute the 95% confidence interval for the permeance data and
% for tau and beta, respectively.
cerrfac = 1/len + (1./(cnu*ckn_nu) - rxmean).^2/rsumsq;
student = tinv(0.975, len-2);
% For the regression, tau was set to 1 and the regression was done on
% perm/(kappa*kn_nu). The membrane mem available here has tau and beta set from
% the results of the regression. Hence, mem.kappa = the original kappa/tau. For
% the confidence interval, mem.kappa must be multiplied by tau to get the
% original kappa.
crstd = student * mem.kappa * mem.tau * ckn_nu * rstd * sqrt(cerrfac);
cmstd = student * mem.kappa * mem.tau * ckn_nu * rstd * sqrt(1+cerrfac);

cx = [cp;NaN;cp]*1e-5;	% Use NaN or Inf to put a gap into a line.
plot(pmean*1e-5,perm,'k*', cp*1e-5,cperm,'k-',...
     cx,[cperm+crstd;NaN;cperm-crstd],'--',...
     cx,[cperm+cmstd;NaN;cperm-cmstd],':');
xlabel('p_{mean} [bar]');
ylabel('permeance [sec]');
