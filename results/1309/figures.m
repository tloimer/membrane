function figures1309()
%FIGURES1309 Create figures for presumably J. Membr. Sci. 2014.

if ~exist('substance.m'), addpath('../../program'); end

T1 = 293.15;
but = substance('butane');
%iso = substance('isobutane');
%n2 = substance('nitrogen');
f = fmodel('plug');

pu1 = membrane(10e-9,0.6,1.38,'tube',3,8.1,20e-6);
pu2 = membrane(100e-9,0.6,1.38,'tube',3,8.1,400e-6);
pu3 = membrane(6e-6,0.6,1.38,'tube',3,8.1,2e-3);

psmems = {{pu1 pu2 pu3}}; pms = mstackstruct(0,psmems,f); pmsorig = pms;
prmems = {{pu3 pu2 pu1}}; pmr = mstackstruct(0,prmems,f); pmrorig = pmr;

psat1 = but.ps(T1); % approx. 2.7e

poben = [1.1:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

% Plot with p1 - p2 = 0.5 bar.
for i = 1:length(poben)
  mf(i) = mnumadiabat(T1,poben(i),poben(i)-0.5e5,but,pms);
  mr(i) = mnumadiabat(T1,poben(i),poben(i)-0.5e5,but,pmr);
end

twoplots(1);

poben = [1.1:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

for i = 1:length(poben)
  mf(i) = mnumadiabat(T1,poben(i),poben(i)-1e5,but,pms);
  mr(i) = mnumadiabat(T1,poben(i),poben(i)-1e5,but,pmr);
end

twoplots(3);

% Plot with p1 - p2 = 1.5 bar.

poben = [1.6:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.6:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

for i = 1:length(poben)
  mf(i) = mnumadiabat(T1,poben(i),poben(i)-1.5e5,but,pms);
  mr(i) = mnumadiabat(T1,poben(i),poben(i)-1.5e5,but,pmr);
end

twoplots(5);

% Default PaperUnits are centimeters
fprintf(['In the eps-files, change\n'...
	'0 cap  to 1 cap (approx. line 156),\n'...
	'and approx. at line 68, change\n'...
	'/DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef  to\n'...
	'/DO { [.1 dpi2point mul 2 dpi2point mul] 0 setdash } bdef\n']);

status = unix(['sed -i ''20,220 {/^0 cap/s//1 cap/;/DO { \[\.5 dpi2point mul 4 dpi2point/'...
	       's//DO { [.1 dpi2point mul 2 dpi2point/}'' figure?.eps']);
if ~status, fprintf('Changed the eps-files!'); end

%- NESTED FUNCTIONS --------------------------------------- NESTED FUNCTIONS ---

function twoplots(i) %------------------------------------------------- twoplots
% the figure is always less wide than given! add 6 mm, instead of
% PaperPosition [0 0 8.9 5]; Perfectly 16x9 would be: 8.8 4.95
% Plot with p1 - p2 = 1 bar.
figure('PaperPosition',[0 0 9.5 5],'PaperSize',[9.5 5]);
hl = plot(poben/psat1,mf,'k-',poben/psat1,mr,'k:');
set(gca,'FontName','Times','FontSize',7,'OuterPosition',[0 0 1 1]);
set(hl(1),'LineWidth',0.3);
set(hl(2),'LineWidth',0.5);
xlabel('{\it p}_{\fontsize{6}1}/{\it p}_{\fontsize{6}sat}');
ylabel('mass flux [gm^{\fontsize{6}-2}s^{\fontsize{6}-1}]');
%ylim([0.004 0.016]);
legend('separation layer upstream (flow direction A)',...
       'separation layer downstream (flow direction B)')
legend('Location','NorthWest'); legend('boxoff');
print('-deps2',sprintf('figure%u.eps',i));

figure('PaperPosition',[0 0 9.5 5],'PaperSize',[9.5 5]);
hl = plot(poben/psat1,mf./mr,'k-');
set(gca,'FontName','Times','FontSize',7);
set(hl(:),'LineWidth',0.3);
xlabel('{\it p}_{\fontsize{6}1}/{\it p}_{\fontsize{6}sat}');
ylabel('mass flux ratio flow direction A/B');
print('-deps2',sprintf('figure%u.eps',i+1));
end %-------------------------------------------------------------- end twoplots

end %%% END FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FIGURES %%%
