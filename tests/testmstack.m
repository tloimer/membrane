function testmstack
%TESTMSTACK Test MSTACK with 10 randomly selected configurations.

%global VERBOSE

% Read the membrane data, copied from read1103.m.
fid = fopen('/home/tloimer/projects/11jms/data/membranes11.tsv');
data = textscan(fid,'%s%n%n%n%n%s%n%n%n%n%n%n',...
  'HeaderLines',1,'ReturnOnError',0,'Delimiter','\t');
fclose(fid);
[memname,poredia,L,memdia,eps,model,tau,beta,rstd,rlen,rxmean,rsumsq] ...
   = deal(data{:});
% scale membrane properties to base SI-units:
poredia = 1e-9*poredia;
L = 1e-3*L;
area = 1e-6*memdia.^2*pi/4;
eps = 1e-2*eps;
clear fid data;

% Set some conditions.
T1 = 295;
theta = 0;
f = fmodel('homogeneous');
s(2) = substance('butane'); psat(2) = s(2).ps(T1);
s(1) = substance('isobutane'); psat(1) = s(1).ps(T1);
p12 = min(0.6*psat(2),1e5);

% The stack has three membranes.
n = 3;
% Already distribute theta and f, otherwise it is distributed in mstack.
theta(1,1:n) = theta;
f(1,1:n) = f;

% Test a selected case which threw an errors.
% T4 in flow12>flow45 was empty, because p > pcrit
% 1    Membranes  5  7 10:  11,  55, 100 nm: HP11A, HP55B, DTN05K14.
mems = [5 7 10];
for i = 3:-1:1
  mm = mems(i);
  mem(i) = membrane(poredia(mm),eps(mm),1.38,model{mm},tau(mm),beta(mm),L(mm));
end
for i = 1:2
  %disp(s(i).name);
  m =  mstack(T1,psat(i),psat(i)-p12,theta,s(i),mem,f);
end

pred = [0.6:0.05:0.85 0.9:0.02:1];

allmembranes = size(poredia,1); % The number of all available membranes.
m = [pred;pred]; % Pre-allocate the mass flux.
lenpred = size(pred,2);
mems = zeros(lenpred,n);
for j = 1:lenpred
  % find a random membrane configuration of three different membranes
  mems(j,:) = randperm(allmembranes,n);
  [~,mind] = sort(poredia(mems(j,:)));
  mems(j,:) = mems(j,mind(:)); % Sort membrane indices by pore diameter.
% DEBUG
%fprintf('%-4.2g Membranes',pred(j)); fprintf(' %2u',mems(j,:));
%fprintf(': %3.0f',poredia(mems(j,1))*1e9);
%fprintf(', %3.0f',poredia(mems(j,2:n))*1e9);
%fprintf(' nm: %s',memname{mems(j,1)}); fprintf(', %s',memname{mems(j,2:n)});
%fprintf('.\n');
% END DEBUG
  for i = n:-1:1
    mm = mems(j,i);
    mem(i) =membrane(poredia(mm),eps(mm),1.38,model{mm},tau(mm),beta(mm),L(mm));
  end
  for i = 1:2
    p1 = pred(j)*psat(i);
    [m(i,j) mguess] = mstack(T1,p1,p1-p12,theta,s(i),mem,f);
    m(i,j) = m(i,j)/mguess;
  end
end

% Plot it.
figure('Name',mfilename);
plot(pred,m(2,:),'kx',pred,m(1,:),'k+');
legend('butane','isobutane','Location','Best'); legend('boxoff');
xlabel('p_1/p_{sat}');
ylabel('mflux/mgas');
xlim([0.55 1.05]);

% Print information on membranes.
fprintf('pred  Membrane stack\n');
for j = 1:lenpred
  fprintf('%-4.2g  %3.0f',pred(j),poredia(mems(j,1))*1e9);
  fprintf(', %3.0f',poredia(mems(j,2:n))*1e9);
  fprintf(' nm: %s',memname{mems(j,1)}); fprintf(', %s',memname{mems(j,2:n)});
  fprintf('.\n');
end

print('-deps2',[mfilename '.eps']);
