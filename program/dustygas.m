function mems = dustygas(mem,data,Tin,sname)
%DUSTYGAS   Compute membrane properties from mass flux data.
%  DUSTYGAS(MEMBRANES,DATA) computes the resistance properties of the
%  individual MEMBRANEs in the cell array MEMBRANES from nitrogen permeation
%  DATA. The array DATA consists of 4 or 5 columns,
%   1  number of membranes
%   2  upstream pressure / Pa
%   3  downstream pressure / Pa
%   4  mass flux / kg/m2/s
%   5  Temperature / K
%  Column 5 may be omitted, if the temperature for all measurements is given
%  as the third argument to DUSTYGAS, see below. DUSTYGAS returns a cell
%  array of membranes, with TAU and BETA set from the calculations. TAU and
%  BETA are really computed from the dusty gas parameters K0 and B0, which
%  are found from linear regression on the mass flow data.
%
%  MEMBRANES must be given in reverse, if lines with % Forward are
%  un-commented in the source. With MEMBRANES = {MEM1 MEM2 MEM3}, the data
%  through one membrane corresponds to the experiment p1 - MEM1 - p2,
%  through two membranes to p1 - MEM2 MEM1 - p2, and through three membranes
%  to p1 - MEM3 MEM2 MEM1 - p2.
%
%  If lines with % Backflow are un-commented, then
%  with MEMBRANES = {MEM1 MEM2 MEM3}, the data through one membrane
%  corresponds to the experiment  p1 - MEM1 - p2, through two membranes to
%  p1 - MEM1 MEM2 - p2, and finally to p1 - MEM1 MEM2 MEM3 - p2.
%
%  If the global variable VERBOSE is greater than 0, the data and the best
%  fit for TAU and BETA is printed.
%
%  DUSTYGAS(MEMBRANES,DATA,T) uses the default Temperature T if it is
%  not given in the DATA array.
%
%  DUSTYGAS(MEMBRANES,DATA,T,NAME) computes the membrane properties
%  from permeation data with the substance NAME. See SUBSTANCE for possible
%  NAMEs.
%
%  See also MEMBRANE, SUBSTANCE.

global VERBOSE;

if nargin < 4
    sname = 'nitrogen';
end
s = substance(sname);

if nargin > 2 && not(isempty(Tin))
    if size(data,2) < 5
	data(:,5) = Tin;
    elseif size(data,2) == 5
	data(isempty(data(:,5)),5) = Tin;
    end
end

if not(iscell(mem))
    error('The membrane must be given as a cell array of membrane(s).');
end

n = length(mem);

% count backwards, to preallocate cell arrays - although this only
% yields a minor performance advantage
for i = n:-1:1
    thislayer = data(:,1) == i;
    p1{i} = data(thislayer, 2);
    p2{i} = data(thislayer, 3);
    m{i} = data(thislayer, 4);	% kg/m2/s
    T{i} = data(thislayer, 5);
    % auxiliary variables; we have substance.R = R/M.
    a{i} = 4*sqrt(8./(pi*s.R*T{i}))/3;
    b{i} = 1./(s.mug(T{i}).*s.R.*T{i});
end

% With the mass flow modeled from the dusty gas model (Hussain et. al, J.
% Porous Media 12, 2009; Thomas et al., Catalysis Today 67,2001),
%
%              1       4          8  M   dp
%   mflow = -( -- B0 + - K0 sqrt(------) --,
%              nu      3         pi R T  dz
%
% the one-dimensional mass flow through a plane, homogeneous porous body is
% given by
%
%   mflow L     M        p1 + p2   4          8 M
%   ------- = ------- B0 ------- + - K0 sqrt(------),
%   p1 - p2   eta R T       2      3         pi R T
%
% where eta is the dynamic viscosity and the ideal gas law, p/rho = R T/M,
% was used.
%
% With the abbreviations a1 and b1,
%
%        4       8 M               M
%   a1 = - sqrt(------),    b1 = -------,
%        3      pi R T           eta R T
%
% for one layer, regress on K0 and B0 for the equation
%
%   mflow L1                   p1 + p2
%   -------- = a1 K01 + b1 B01 -------.					(1)
%   p1 - p2                       2
%
% Here a1, b1 and L1 and the data tuples (mflow, p1, p2) are known.
%
% For the next layer,
%
%   mflow L1 = a1 K01 p1 - a1 K01 p2 + b1 B01 p1^2/2 - b1 B01 p2^2/2,
%
%   0.5 b1 B01 p2^2 + a1 K01 p2 + mflow L1 - a1 K01 p1 - 0.5 b1 B01 p1^2 = 0,
%
% hence p2 is given by
%
%          a1 K01          a1 K01      2   2 mflow L1
%   p2 = - ------ + sqrt( (------ + p1)  - ---------- ).		(2)
%          b1 B01          b1 B01            b1 B01
%
% With the data tuples (mflow, p2, p3), get the unknwon properties K02 and
% B02 for the second layer by regressing on
%
%   mflow L2                   p2 + p3
%   -------- = a2 K02 + b2 B02 -------.
%   p2 - p3                       2
%
% For forward flow, if the downstream-most layer is the support layer,
%
%          a1 K01          a1 K01      2     2 mflow L1
%   p1 = - ------ + sqrt( (------ + p2)   +  ---------- ).		(3)
%          b1 B01          b1 B01              b1 B01


% p1 ... upstream pressure,
% pu ... pressure upstream of the layer to be determined
% p2 ... downstream pressure
%   p1  ...properties known ... pu               p2
%      | layer 1 | layer 2 | ... |  last layer |
%
% Forward flow
% pu ... pressure _downstream_ of the layer to be determined
%    p1            pu  ...properties known ...   p2
%      | last layer | ... | layer 2 |  layer 1 |

K0 = zeros(1,n);
B0 = zeros(1,n);
for j = 1:n
%fprintf('\n\n #### Membrane %d ####\n\n', j);		% DEBUG
%    pu = p1{j};		% Backward
    pu = p2{j};			% Forward
    for i = 1:j-1
	c1 = a{j}.*K0(i)./b{j}/B0(i);
%	pu = -c1 + sqrt((c1 + pu).^2 - 2*m{j}.*mem{i}.L./(b{j}.*B0(i))); % Back
	pu = -c1 + sqrt((c1 + pu).^2 + 2*m{j}.*mem{i}.L./(b{j}.*B0(i))); % Forw
    end

%    pmean = 0.5*(pu + p2{j});				% Backward
%    [hat,~,mse] = lscov([a{j} b{j}.*pmean], m{j}*mem{j}.L./(pu-p2{j}));
    pmean = 0.5*(pu + p1{j});				% Forward
    [hat,~,mse] = lscov([a{j} b{j}.*pmean], m{j}*mem{j}.L./(p1{j}-pu));
    K0(j) = hat(1);
    B0(j) = hat(2);
%fprintf(' K0(%d) = %.6g, B0(%d) = %.6g\n', j, K0(j), j, B0(j));	% DEBUG
    if VERBOSE > 0
	pdelta = 0.05*(max(pmean) - min(pmean));
	pm = (min(pmean)-pdelta:pdelta:max(pmean)+1.1*pdelta)';
	Tmean = mean(T{j});
	%mhat = a{j}*K0(j) + b{j}.*B0(j).*pm;
	mhat = 4*sqrt(8/(pi*s.R*Tmean))*K0(j)/3 ...
	    + B0(j)*pm/(s.mug(Tmean)*s.R*Tmean);
	N = size(pm,1);
	meanp = mean(pm);
	varp = var(pm);
	conf_mean = @(xp) sqrt(mse*(1/N + (xp-meanp).^2/(N-1)/varp));
	conf_ind = @(xp) sqrt(mse*(1 + 1/N + (xp-meanp).^2/(N-1)/varp));
	figure('Name',sprintf('Layer %d',j));
% Backward
%	plot(pmean*1e-5, m{j}*mem{j}.L./(pu-p2{j}), 'ok', pm*1e-5, mhat,'-', ...
% Forward
	plot(pmean*1e-5, m{j}*mem{j}.L./(p1{j}-pu), 'ok', pm*1e-5, mhat,'-', ...
	    [pm;NaN;pm]*1e-5, [mhat;NaN;mhat] + tinv(0.975,N-2)...
		* [conf_mean(pm); NaN; -conf_mean(pm)], '--', ...
	    [pm;NaN;pm]*1e-5, [mhat;NaN;mhat] + tinv(0.975,N-2)...
		* [conf_ind(pm); NaN; -conf_ind(pm)], '--');
	xlabel('p_{mean} / bar');
	ylabel('mflux*L/(p_1-p_2) / s');
    end
end

% Write the results, cf. file notes/B0k0_vs_taubeta.pdf,
%
%   kappa = B0,  tau = ftau(kappa) (mem.dia, mem.epsilon),
%   beta = 32 k0 dia / (9 pi B0)
%
% Re-initialize the membranes, because the function definitions retain
% values that were in force when these functions were defined.
for i = n:-1:1
    mems{i} = membrane(mem{i}.dia, mem{i}.epsilon,mem{i}.km, mem{i}.tname,...
		mem{i}.ftau(B0(i)), 32*K0(i)*mem{i}.dia/(9*pi*B0(i)), mem{i}.L);
end
