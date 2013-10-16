function mflow()
%MFLOW       Create figures for presumably J. Membr. Sci. 2014.

if ~exist('substance.m'), addpath('../../program'); end

%Setup the flow problem
T1 = 293.15;
%s = substance('butane');
s = substance('isobutane');
%s = substance('nitrogen');
f = fmodel('plug');

pu1 = membrane(10e-9,0.6,36,'tube',3,8.1,20e-6);
pu2 = membrane(100e-9,0.6,36,'tube',3,8.1,150e-6);
pu3 = membrane(6e-6,0.6,36,'tube',3,8.1,2e-3);

psmems = {{pu1 pu2 pu3}}; pms = mstackstruct(0,psmems,f); pmsorig = pms;
prmems = {{pu3 pu2 pu1}}; pmr = mstackstruct(0,prmems,f); pmrorig = pmr;

psat1 = s.ps(T1);

% Plot with p1 - p2 = 0.1 bar.

poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

deltap = 0.1e5;
for i = 2:length(poben)
  mf(i) = mnumadiabat(T1,poben(i),poben(i)-deltap,s,pms);
  mr(i) = mnumadiabat(T1,poben(i),poben(i)-deltap,s,pmr);
end
[mf(1),pms] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pms);
[mr(1),pmr] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pmr);

fprintf('\nDirection separation - support layer');
printpcondensation(pms);
fprintf('Direction support - separation layer');
printpcondensation(pmr);

twoplots(1);

% Plot with p1 - p2 = 0.5 bar.

poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

deltap = 0.5e5;
for i = 2:length(poben)
  mf(i) = mnumadiabat(T1,poben(i),poben(i)-deltap,s,pms);
  mr(i) = mnumadiabat(T1,poben(i),poben(i)-deltap,s,pmr);
end
[mf(1),pms] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pms);
[mr(1),pmr] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pmr);

fprintf('Direction separation - support layer');
printpcondensation(pms);
fprintf('Direction support - separation layer');
printpcondensation(pmr);

twoplots(3);

poben = [2.1:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.1:0.2:psat1/1e5 psat1/1e5]*1e5;
mf = poben; mr = poben;

deltap = 1e5;
for i = 2:length(poben)
  mf(i) = mnumadiabat(T1,poben(i),poben(i)-1e5,s,pms);
  mr(i) = mnumadiabat(T1,poben(i),poben(i)-1e5,s,pmr);
end
[mf(1),pms] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pms);
[mr(1),pmr] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pmr);

fprintf('\nDirection separation - support layer');
printpcondensation(pms);
fprintf('Direction support - separation layer');
printpcondensation(pmr);

twoplots(5);

% Plot with p1 - p2 = 1.5 bar.

%%poben = [1.6:0.01:psat1/1e5 psat1/1e5]*1e5;
%poben = [1.6:0.2:psat1/1e5 psat1/1e5]*1e5;
%mf = poben; mr = poben;
%
%deltap = 1.5e5;
%for i = 2:length(poben)
%  mf(i) = mnumadiabat(T1,poben(i),poben(i)-1.5e5,s,pms);
%  mr(i) = mnumadiabat(T1,poben(i),poben(i)-1.5e5,s,pmr);
%end
%[mf(1),pms] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pms);
%[mr(1),pmr] = mnumadiabat(T1,poben(1),poben(1)-deltap,s,pmr);
%
%fprintf('\nDirection separation - support layer');
%printpcondensation(pms);
%fprintf('Direction support - separation layer');
%printpcondensation(pmr);
%
%twoplots(5);


% Default PaperUnits are centimeters
fprintf(['In the eps-files, change\n'...
	'0 cap  to 1 cap (approx. line 156),\n'...
	'and approx. at line 68, change\n'...
	'/DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef  to\n'...
	'/DO { [.1 dpi2point mul 2 dpi2point mul] 0 setdash } bdef\n']);

command = ['sed -i ''20,220 {/^0 cap/s//1 cap/;/DO { \[\.5 dpi2point mul 4 dpi2point/'...
	      's//DO { [.1 dpi2point mul 2 dpi2point/}'' ' mfilename '?.eps'];
status = unix(command);
if status
  fprintf(['sed-command failed\n' command '\n']);
else
  fprintf('Changed the eps-files!\n');
end

fprintf('# xelatex %1$s.tex && rm %1$s.log %1$s.aux\n',mfilename);

%- NESTED FUNCTIONS --------------------------------------- NESTED FUNCTIONS ---

function twoplots(i) %------------------------------------------------- twoplots
% the figure is always less wide than given! add 6 mm, instead of
% PaperPosition [0 0 8.9 5]; Perfectly 16x9 would be: 8.8 4.95
% Plot with p1 - p2 = 1 bar.
figure('PaperPosition',[0 0 9.9 5],'PaperSize',[9.9 5]);
hl = plot(poben/psat1,mf*1e3,'k-',poben/psat1,mr*1e3,'k:');
set(gca,'FontName','Times','FontSize',7,'OuterPosition',[0 0 1 1]);
set(hl(1),'LineWidth',0.3);
set(hl(2),'LineWidth',0.5);
xlabel('{\it p}_{\fontsize{6}1}/{\it p}_{\fontsize{6}sat}');
ylabel('mass flux [gm^{\fontsize{6}-2}s^{\fontsize{6}-1}]');
xlim([0.7 1]);
%ylim([0.004 0.016]);
legend('separation layer upstream (flow direction A)',...
       'separation layer downstream (flow direction B)')
legend('Location','NorthWest'); legend('boxoff');
print('-deps2',sprintf('%s%u.eps',mfilename,i));

figure('PaperPosition',[0 0 9.9 5],'PaperSize',[9.9 5]);
hl = plot(poben/psat1,mf./mr,'k-');
set(gca,'FontName','Times','FontSize',7);
set(hl(:),'LineWidth',0.3);
xlabel('{\it p}_{\fontsize{6}1}/{\it p}_{\fontsize{6}sat}');
ylabel('mass flux ratio flow direction A/B');
xlim([0.7 1]);
print('-deps2',sprintf('%s%u.eps',mfilename,i+1));
end %-------------------------------------------------------------- end twoplots

function printpcondensation(ms) %---------------------------- printpcondensation
% Print relative pressure values from linear theory, when condensation starts
[dpk,pk] = ms.membrane(1).layer(1).flsetup.dpkdT(T1);
dpfl = ms.membrane(1).layer(1).flow(end).p(end) ...
       - ms.membrane(1).layer(1).flow(1).p(1);
n = ms.substance.jt(T1,poben(1)) * dpk;
fprintf(', delta p = %.1f bar.\n',deltap/1e5);
fprintf(['  pk/psat = %.2f, (dT/dp)_h(dpk/dT) = %.2f, (p1-p2)/psat = %.2f,\n'...
  '  pdiff(first layer)/psat = %.2f. pcrel(p1-p2) = %.2f, pcrel(fl) = %.2f.\n'],...
  pk/psat1,n,deltap/psat1,dpfl/psat1,(pk-n*deltap)/psat1,(pk-n*dpfl)/psat1);
end %---------------------------------------------------- end printpcondensation

end %%% END2.1GURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FIGURES %%%
