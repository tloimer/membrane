function map(T)
%MAP        Print the flow map.
%  MAP(T) prints the flow map for the temperature T.
%
%  Calls K, KL, NUG, NUL, PS, R.

% some variables
%Ccc,p1pk_2s,jtdpkdT

[R M]=molm;
pspk_2s = ps(T)*M/(R*rho(T)*T);
jtdpk = dpkdT(T)*jt(T,ps(T));
% the simpler variant, dpk/dT = dps/dT (which is very nearly true!)
% jtdpk = dpkdT(T).*jt(T,ps(T));

% Ccaps are caculated as the inverse, since they will be plotted on
% inverted axis
Ccc = (1-jtdpk)/pspk_2s;
kfkc = jtdpk/(1+Ccc);
% Now curves 1 and 2
% the x-coordinates of the plot points. logarithmic scaling!
kkc1 = log(kfkc) : -log(kfkc)/10 : 0;
kkc1 = [kfkc exp(kkc1(2:10)) 1];
kkc2 = 0 : log(4)/4 : log(4);
kkc2 = [1 exp(kkc2(2:end))];
C1 = (1-kkc1)./(kkc1+pspk_2s);
C2 = (1+Ccc)./(1 + pspk_2s*(1+jtdpk./kkc2));

% plot
% for a logarithmic plot, insert a small number, e.g.: [0 0.001 kkc2(end)]
% version without scaling
%plot([0 kkc2(end)],[0 0],'k',[0 kkc2(end)],[-1 -1],'k',[1 1],[-1 C2(1)],'k',...
%  [0 kkc2(end)],[Ccc Ccc],'k',[kfkc kfkc],[Ccc 5*Ccc],'k',...
%  kkc1,C1,'r',kkc2,C2,'k');
%ylim([-2 1.5*Ccc]);
% the wetted region is scaled
plot([0 kkc2(end)],[0 0],'k',[0 kkc2(end)],[-1 -1],'k',...
  [1 1],[-1 C2(1)/Ccc],'k',...
  [0 kkc2(end)],[1 1],'k',[kfkc kfkc],[1 5],'k',...
  kkc1,C1/Ccc,'k',kkc2,C2/Ccc,'k');
xlim([0 1.5]);
ylim([-2 1.5]);
set(gca,'YTick',[-1 0 1/Ccc 1]);
set(gca,'YTickLabel',{'-1','0','1','Ccc'})


