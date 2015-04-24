%ENTROPY Compare entropy calculated in substance.m with tabulated values from
%  NIST, http://webbook.nist.gov/chemistry/fluid .
if ~exist('substance.m')
  addpath('../program');
end
% Test entropy
butgsT300 = 1e3*[2.9120 2.8119 2.7529 2.7108 2.6779 2.6509 2.6278 2.6077 ...
	2.5899 2.5738 2.5591 2.5456 2.5331 2.5215 2.5105 2.5002 2.4905 ...
	2.4812 2.4724 2.4639 2.4558 2.4480 2.4405 2.4333 2.4262 2.4211];

butlsT300 = 1e3*[1.2225 1.2224 1.2223 1.2221 1.2219 1.2218 1.2216 1.2214 ...
	1.2212 1.2210 1.2209 1.2207 1.2205 1.2203 1.2201 1.2200 1.2198];

butgpT300 = 1e5*[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 ...
	1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.576];

butlpT300 = 1e5*[ 2.576 2.6 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 ...
	9.0 9.5 10];

b = substance('butane');
psat300 = b.ps(300);
sgref = b.s(300,psat300,1,300);
slref = b.s(300,psat300,0,300);

close all;

figure('Name',mfilename);
plot(butgpT300, butgsT300-butgsT300(end),'-+',...
     butgpT300, b.s(300,butgpT300,1,300)-sgref,'-o');
ylabel('s [J/kgK]');
xlabel('p [bar]');
legend('Nist','matlab');
title('butane, 300 K, vapor phase');

figure('Name',mfilename);
plot(butlpT300, butlsT300-butlsT300(1),'-+',...
     butlpT300, b.s(300,butlpT300,0,300)-slref,'-o');
ylabel('s [J/kgK]');
xlabel('p [bar]');
legend('Nist','matlab');
title('butane, 300K, liquid phase');
