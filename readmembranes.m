%   readdata schreibt folgende Variablen:
% memname,poredia,L,memdia,eps,model,tau,beta
% substancename datamemname exp_id T1 p1 p2 Vflow Troom T1tc T2tc T12tc

fid = fopen('../data/membranesguessed');
%fid = fopen('../data/membranes');
% membranes:
% membrane material wetting poredia[nm] L[mm] memdia[mm] eps[%] model tau beta
% To ignore remaining lines in a field: %*[^\n]
%data = textscan(fid,'%s%*s%*s%n%n%n%n%s%n%n%*[^\n]',...
data = textscan(fid,'%s%*s%*s%n%n%n%n%s%n%n',...
  'HeaderLines',1,'ReturnOnError',0,'Delimiter','\t');
fclose(fid);
[memname,poredia,L,memdia,eps,model,tau,beta] = deal(data{:});

% scale membrane properties to base SI-units:
% membrane
poredia = 1e-9*poredia;
L = 1e-3*L;
area = 1e-6*memdia.^2*pi/4;
eps = 1e-2*eps;
% data points

clear fid data;
