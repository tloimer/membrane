function f2ph1ph = flow2ph1ph(m,T0,x0,q0,zrange)
%FLOW2PH1PH(M,T0,X0,Q0,ZRANGE) returns the solution for an 2ph flow with
%  evaporation front within the membrane.
%  The result is a structure F2PH1PH(I) that always contains these
%  fields: .z, .t, .p, .q., .x, .color.

% first get the 2ph-flow
f2ph=int2ph(m,T0,x0,q0,zrange);

% the remaining heat flux, if all evaporates
% q2ph=qevap(m,f2ph.y(1,:),f2ph.y(2,:));

% now integrate the 1ph-region backward, starting from Tend=T0
% fgb = intg(m,T0,ps(T0),0,[zrange(2) zrange(1)]);
% ygb = deval(fgb,f2ph.x)
% ygb = deval(intg(m,T0,ps(T0),0,[zrange(2) zrange(1)]),f2ph.x);

% and find the point, where the q's are equal
% zevap =interp1q((q2ph-ygb(3,:))',f2ph.x',0)

dz = abs(zrange(2)-zrange(1))/1000;
origopts = optimset('fzero');
options = optimset(origopts,'TolX',dz);
zevap = fzero(@qend,[zrange(1) zrange(2)-dz],options,m,f2ph,zrange(2));

% now integrate 1-ph to the end
yevap = deval(f2ph,zevap);
fg = intg(m,yevap(1),ps(yevap(1)),qevap(m,yevap(1),yevap(2)),...
       [zevap zrange(2)]);

% concatenate for the final result
i2ph = find(f2ph.x<zevap);
z2ph = f2ph.x(i2ph);
T2ph = f2ph.y(1,i2ph);
x2ph = f2ph.y(2,i2ph);
p2ph = ps(f2ph.y(1,i2ph));
q2ph = m*q_m(f2ph.y(1,i2ph),f2ph.y(2,i2ph));

% construct the result struct
f2ph1ph = struct('z',{[z2ph zevap],fg.x},...
  'T',{[T2ph yevap(1)],fg.y(1,:)},...
  'p',{[p2ph ps(yevap(1))],fg.y(2,:)},...
  'q',{[q2ph m*q_m(yevap(1),yevap(2))],fg.y(3,:)},...
  'x',{[x2ph yevap(2)],ones(size(fg.x))},...
  'color',{'g','r'});

%-----------------------------------------------------------------------
function qend = qend(zevap,m,f2ph,zend)
y2ph=deval(f2ph,zevap);
fg = intg(m,y2ph(1),ps(y2ph(1)),qevap(m,y2ph(1),y2ph(2)),[zevap zend]);
qend = fg.y(3,end);
