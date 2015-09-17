%   readdata schreibt folgende Variablen:
% memname,poredia,L,memdia,eps,model,tau,beta
% substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc

fid = fopen('data/membranes11.tsv');
%fid = fopen('../data/membranes');
% membranes:
% membrane material wetting poredia[nm] L[mm] memdia[mm] eps[%] model tau beta
% To ignore remaining lines in a field: %*[^\n]
%data = textscan(fid,'%s%*s%*s%n%n%n%n%s%n%n%*[^\n]',...
data = textscan(fid,'%s%n%n%n%n%s%n%n%n%n%n%n',...
  'HeaderLines',1,'ReturnOnError',0,'Delimiter','\t');
fclose(fid);
[memname,poredia,L,memdia,eps,model,tau,beta,rstd,rlen,rxmean,rsumsq] ...
   = deal(data{:});
%data:
% substance membrane exp_id T1[C] p1 p2 Q[ml/min] Troom[C]
fid = fopen('data/data0905.tsv');
data = textscan(fid,'%s%s%s%n%n%n%n%n%n%n%n',...
  'HeaderLines',1,'ReturnOnError',0,'Delimiter','\t');
fclose(fid);
[substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc]...
  = deal(data{:});

% scale membrane properties to base SI-units:
% membrane
poredia = 1e-9*poredia;
L = 1e-3*L;
area = 1e-6*memdia.^2*pi/4;
eps = 1e-2*eps;
% data points
T1 = T1 + 273.15;
p1 = p1*1e5;
p2 = p2*1e5;
Troom = Troom + 273.15;
% Vflow remains ml/min

clear fid data;
