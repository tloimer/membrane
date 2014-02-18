
if ~exist('substance.m'), addpath('../program'); end

T1 = 293.15;
%but = substance('butane');
s = substance('isobutane');
%n2 = substance('nitrogen');
f = fmodel('plug');

pu1 = membrane(10e-9,0.6,36,'tube',3,8.1,20e-6);
pu2 = membrane(100e-9,0.6,36,'tube',3,8.1,150e-6);
pu3 = membrane(6e-6,0.6,36,'tube',3,8.1,2e-3);

psmems = {{pu1 pu2 pu3}}; pms = mstackstruct(0,psmems,f); pmsorig = pms;
prmems = {{pu3 pu2 pu1}}; pmr = mstackstruct(0,prmems,f); pmrorig = pmr;

psat1 = s.ps(T1); % approx. 2.7e

prel = 0.95;
p1 = prel*psat1;
[mf,pms] = mnumadiabat(T1,p1,p1-0.5e5,s,pms);
pms.plotsolution(pms);
xlabel('z/L');
ylabel('p [bar]');

pms.plotT(pms);
%print('-depsc2','figure1p.eps');
%setaxes = @(ah) set(ah,'Units','points','Position',[33 28 216 128],...
%			'FontName','Times','FontSize',7);
setaxes = @(ah) set(ah,'FontName','Times','FontSize',11);
setfig = @(fh) set(fh,'PaperUnits','points','PaperPosition',[0 0 254 174],...
			'PaperSize',[254 174],'OuterPosition',[100 100 330 280]);

% finally, for the pressure plot
% >> set(gcf,'PaperPosition',[0 0 254 144],'PaperSize',[254 144])
% >> set(gca,'Position',[0.62 0.13 0.35 0.84])
% >> set(gca,'Position',[0.13 0.13 0.21 0.84])
% >> set(gca,'Position',[0.34 0.13 0.28 0.84])
