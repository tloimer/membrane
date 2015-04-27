%TESTCALCTAUBETA Compare data from membranes11.tsv with current calculation.

if ~exist('substance.m'), addpath('../program'); end

% read1103 only works within a script, not in a function body
read1103;
% Populates following arrays by reading .../data/data0905.tsv:
%   substancename, datamemname, exp_id, T1, p1, p2, Vflow, Troom, T1tc,
%   T2tc, T12tc,
% and by reading .../data/membranes11.tsv:
%   memname, poredia, L, memdia, eps, model, tau, beta, rstd, rlen, rxmean,
%   rsumsq
n2 = 'nitrogen';

for i = 1:length(memname)
  mem = membrane(poredia(i),eps(i),1.38,model{i},1,8.1,L(i));
  mem.area = area(i);
  isdata = strcmp(memname{i},datamemname) & strcmp(n2,substancename);
  lendata = sum(isdata);
  if lendata == 0
    fprintf(['No ' n2 ' data found for membrane ' memname{i} '.\n']);
    continue;
  end
  fprintf(['Membrane ' memname{i} ':\n']);
  onesdata = ones(size(T1(isdata)));
  [sarray{1:lendata,1}] = deal(n2);
  data = struct('T1',T1(isdata),'p1',p1(isdata),'p2',p2(isdata),...
		'Troom',Troom(isdata), 'proom',101300*onesdata,...
		'volume',Vflow(isdata)*1e-6,'duration',60*onesdata,...
		'substance',{sarray});
  clear('sarray');
  fprintf('Original values, tau = %.2f, beta = %.2f.\n', tau(i), beta(i));
  figure('Name',['Membrane ' memname{i}]);
  calctaubeta(mem,data);
  fprintf('\n');
end
