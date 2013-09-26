% the two membranes
mem(2) = membrane(50e-9,0.5,1.38,'tube',1,8.1,1e-3);
mem(1) = membrane(10e-9,0.5,1.38,'tube',1,8.1,0.2e-3);

% the fluid properties
theta = 0;
s = substance('isobutane');
T1 = 22 + 273.15;
p1red = [0.2:0.1:0.8 0.81:0.01:1];
p12 = 0.5e5;
p1 = p1red.*s.ps(T1);

% the two-phase flow model
f = fmodel('homogeneous');

mfolio = p1; % the forward configuration
%mverso = p1; % the backward configuration
p3f=p1; %p3v=p1;
mone=p1; mtwo=p1;
len = size(p1,2);
fid = fopen('table10-50nm','w');
fprintf(fid,['Data for a stack p1 - (10 nm, 1 mm) - (50 nm, 0.2 mm) - p2' ...
 ' and backward ( p1 - 50 nm - 10 nm - p2 ).\n\n']);
fprintf(fid,['p1/psat  mforward mback    mf/mb    p1 [bar] p2 [bar] p1-pint  '...
  'flow-forward flow-backward\n']);
  
for i = 1:len
  [mfolio(i) mfgas flf] = mstack(T1,p1(i),p1(i)-p12,theta,s,mem,f);
  [mverso(i) mvgas flv] = mstack(T1,p1(i),p1(i)-p12,theta,s,[mem(2) mem(1)],f);
%  mfolio(i) = mfolio(i)/mfgas;
%  mverso(i) = mverso(i)/mvgas;
  p3f(i) = flf(2).sol.p1;
%  Tint(i) = flf(2).sol.T1;
%  p3v(i) = flv(2).sol.p1;
%  mone(i) = mnum(T1,p1(i),p3f(i),theta,s,mem(1),f,'m',mfolio(i));
%  mtwo(i) = mnum(Tint,p3f(i),p1(i)-p12,theta,s,mem(2),f,'m',mfolio(i));
  fprintf(fid,'%-9.2f%-9.3g%-9.3g%-9.3f%-9.3f%-9.3f%-9.3f%-13s%s\n',p1red(i),...
  mfolio(i),mverso(i),mfolio(i)/mverso(i),p1(i)/1e5,(p1(i)-p12)/1e5,...
  (p1(i)-p3f(i))/1e5,[flf(1).sol.states flf(2).sol.states],...
  [flv(1).sol.states flv(2).sol.states]);
end

fclose(fid);

plot(p1red,mfolio,'k*',p1red,mverso,'bo');
%plot(p1red,T1-Tint,'k*',p1red,(p1-p3f)*2e-5,'ro');
xlabel('p_1/p_{sat}');
%ylabel('T_1 - T_{int} [K], 2*(p_1-p_{int}) [bar]');
ylabel('mass flux [kg/m^2s]');
%legend('10 nm, 0.2 mm and 50 nm, 2mm membrane','10 nm membrane, p_1 - p_{int}','50 nm membrane, p_{int} - p_2');
%title('Mass flux through two membranes, and separately');
title('10 nm and 50 nm membrane, theta = 90');
legend('forward: p1 - 10 nm - 50 nm - p2','backward: p1 - 50 nm - 10 nm - p2');

%figure;
%plot(p1,p3f,'k*',p1,p3v,'bo');
