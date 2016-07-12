% Test DUSTYGAS by creating mass flow data, artificially applying an random
% error and computing tau and beta from the noisy signal.

if ~exist('substance.m','file'), addpath('../program'); end

T1 = 293.15;
s = substance('nitrogen');

% From program/mem_data.m case hermsdorf2014
l1 = membrane(2.5e-6,	0.6,	30,	'tube',	4.8,	7.1,	1e-3);
l2 = membrane(0.8e-6,	0.5,	30,	'tube',	3.05,	9.5,	25e-6);
l3 = membrane(0.4e-6,	0.5,	30,	'tube',	2.66,	9.1,	25e-6);
l4 = membrane(0.1e-6,	0.4,	30,	'tube',	1.6,	8.3,	25e-6);
l5 = membrane(30e-9,	0.3,	11.8,	'tube',	3,	8.1,	5e-6);

mems = {l1 l2 l3 l4 l5}; % Forward
%mems = {l5 l4 l3 l2 l1}; % backward
n = length(mems);

f = fmodel('emt');
theta = 0;

mf{5} = mstackstruct(theta, {{l5 l4 l3 l2 l1}}, f);
mf{1} = mstackstruct(theta, {{l1}}, f);
mf{2} = mstackstruct(theta, {{l2 l1}}, f);
mf{3} = mstackstruct(theta, {{l3 l2 l1}}, f);
mf{4} = mstackstruct(theta, {{l4 l3 l2 l1}}, f);

mb{5} = mstackstruct(theta, {{l1 l2 l3 l4 l5}}, f);
mb{1} = mstackstruct(theta, {{l1}}, f);
mb{2} = mstackstruct(theta, {{l1 l2}}, f);
mb{3} = mstackstruct(theta, {{l1 l2 l3}}, f);
mb{4} = mstackstruct(theta, {{l1 l2 l3 l4}}, f);

%pm = [1:0.4:3]*1e5;
pm = [2 3]*1e5;
dp = 0.5e5;
N = size(pm,2);
mnum = zeros(n,N);
mr = zeros(n,N);

for j = 1:n % or parfor?
    for i = 1:N
	[mnum(j,i) mb{j}] = mgaseous(T1, pm(i)+0.5*dp, pm(i)-0.5*dp, s, mb{j});
	mb{1}.printsolution(mb{j});
    end
end

% Now, randomize the data.
u = -1 + 2*rand(n,N);
% Or, for debugging, rather not
u = zeros(size(u));

for j = 1:n
    mstd = mnum(j,1)*0.1;
    mr(j,:) = mnum(j,:) + mstd * u(j,:);
end

data = zeros(n*N,4);
for j = 1:n
   range = (j-1)*N+1:j*N;
   data(range,1) = j;
   data(range,2) = pm' + 0.5*dp;
   data(range,3) = pm' - 0.5*dp;
   data(range,4) = mr(j,:)';
end

global VERBOSE;
VERBOSE=1;

mem = dustygas(mems,data,T1);

VERBOSE = 0;

% Print the result.
for j = 1:n
    fprintf(['Membrane %d, L =%4.0f mum, dia =%4.0f nm, tau = %4.2f, '...
	'beta = %3.1f.\n  Given: kappa = B0 = %.5g, K0 = %.5g\n'],...
	j, mems{j}.L*1e6, mems{j}.dia*1e9, mems{j}.tau, mems{j}.beta,...
	mems{j}.kappa, mems{j}.beta*9*pi*mems{j}.kappa/(32*mems{j}.dia));
    fprintf(...
	'  Calc.: kappa = B0 = %.5g, K0 = %.5g,  tau = %.3f, beta = %.3f\n',...
	mem{j}.kappa, mem{j}.beta*9*pi*mem{j}.kappa/(32*mem{j}.dia),...
	mem{j}.tau, mem{j}.beta);
    fprintf(...
      'Changes: kappa = %.3g%%, K0 = %.3g%%, tau = %.3g%%, beta = %.3g%%.\n',...
	(mem{j}.kappa/mems{j}.kappa - 1) * 100,...
	((mem{j}.beta*9*pi*mem{j}.kappa/(32*mem{j}.dia)) ...
	 / (mems{j}.beta*9*pi*mems{j}.kappa/(32*mems{j}.dia)) - 1) * 100, ...
	(mem{j}.tau/mems{j}.tau - 1) * 100, (mem{j}.beta/mems{j}.beta - 1)*100);
end
