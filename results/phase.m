function phase()
%PHASE      Plot vapor volume fraction in membrane layers

if ~exist('substance.m'), addpath('../program'); end

%Setup the flow problem
T1 = 293.15;
%s = substance('butane');
s = substance('isobutane');
%s = substance('nitrogen');
f = fmodel('plug');

pu1 = membrane(10e-9,0.5,24,'tube',2,8.1,20e-6);
pu2 = membrane(100e-9,0.4,24,'tube',2.5,8.1,150e-6);
pu3 = membrane(6e-6,0.3,24,'tube',3,8.1,2e-3);

psmems = {{pu1 pu2 pu3}}; pms = mstackstruct(0,psmems,f); pmsorig = pms;
prmems = {{pu3 pu2 pu1}}; pmr = mstackstruct(0,prmems,f); pmrorig = pmr;

psat1 = s.ps(T1);

%produce(0.1e5,0.836);
%produce(0.1e5,0.85);
%produce(0.1e5,0.9);
%produce(0.1e5,0.99);
%produce(0.1e5,1);

produce(0.5e5,0.84);
produce(0.5e5,0.9);
produce(0.5e5,0.95);
produce(0.5e5,0.97);
produce(0.5e5,1);

%produce(1e5,0.78);
%produce(1e5,0.83);
%produce(1e5,0.88);
%produce(1e5,0.93);
%produce(1e5,1);


%-- NESTED FUNCTIONS -------------------------------------- NESTED FUNCTIONS ---

function produce(deltap,prel) %----------------------------------------- produce

% Construct a file name
if prel == 1
  id = '1';
else
  id = sprintf('%.0f',100*prel);
end
id = [sprintf('%02.0f-',deltap/1e4) id];

p1 = prel*psat1;
pms = pmsorig;
[mf,pms] = mnumadiabat(T1,p1,p1-deltap,s,pms);
plotphase(pms,[id 'A']);
pms.printsolution(pms,'phase');
fprintf('\n');
pmr = pmrorig;
[mr,pmr] = mnumadiabat(T1,p1,p1-deltap,s,pmr);
plotphase(pmr,[id 'B']);
pmr.printsolution(pmr,'phase');
fprintf('\n');
end %--------------------------------------------------------------- end produce

end %%% END PHASE.M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PHASE.M %%%


%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%

function plotphase(ms,id) %------------------------------------------- plotphase
%PLOTPHASE  Plot and export the phase distribution within the membrane.

hp = figure('Name','Phase','PaperPosition',[0 0 9.9 3],'PaperSize',[9.9 3]);
hppos = get(hp,'OuterPosition');
hppos(3) = 2*hppos(3);
set(hp,'OuterPosition',hppos);

% Oh, use a comma-separated list to assign plot style
lprops = {'Color','k','LineStyle','-','LineWidth',0.3};
aprops = {'Ylim',[-0.1 1.1],'FontName','Times','FontSize',7,'Box','on',...
	  'YTick',[0 1],'XTick',[0 1],'TickLength',[0.025 0.025]};

% We assume, only one membrane
nlayers = length(ms.membrane.layer); % nlayers = 3
sp(nlayers) = 0;
sp(1) = axes('Position',[0.11 0.23 0.27 0.76],'XLim',[-0.1 1],'YTickLabel',[0 1],aprops{:});
line([-0.1 0],[1 1],lprops{:});
ylabel('\it\alpha','Rotation',-0);
xlabel('{\it x/L}_{\fontsize{6}1}','VerticalAlignment','baseline')
sp(2) = axes('Position',[0.4  0.23 0.24 0.76],'XLim',[0 1],'YTickLabel',[],aprops{:});
xlabel('{\it x/L}_{\fontsize{6}2}','VerticalAlignment','baseline')
sp(3) = axes('Position',[0.66 0.23 0.27 0.76],'XLim',[0 1.1],'YTickLabel',[],aprops{:});
xlabel('{\it x/L}_{\fontsize{6}3}','VerticalAlignment','baseline')
line([1 1.1],[1 1],lprops{:});
for i = 1:nlayers
  %subplot(1,nlayers,i);
  jl = ms.membrane.layer(i);
  L = jl.matrix.L;
  nflow = length(jl.flow);
  axes(sp(i));
  for j = nflow:-1:1
    line(jl.flow(j).z/L,jl.flow(j).a,lprops{:})
  end
end

print('-dpdf',sprintf('%s%s.pdf',mfilename,id));
fprintf('Export to file %s%s.pdf\n',mfilename,id);
end %------------------------------------------------------------- end plotphase
