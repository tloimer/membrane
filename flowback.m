function flowstruct = flowback(m,Te,pe,L)
%FLOWBACK   Backward integration of flow field.
%  FLOWSTRUCT = FLOWBACK(M,TE,PE,L) integrates the flow field backwards.
%  The temperature distribution for the liquid flow within the membrane
%  and the liquid film upstream of the membrane are solved analytically,
%  where constant material properties are assumed. Returns a structure
%  FLOWSTRUCT that consists of three parts, FLOWSTRUCT.INFO, .SOL and
%  FLOWSTRUCT.FLOW. Writes very incomplete .INFO.
%  For 2-ph flow, FLOW(1).X contains the solution for the vapor region,
%  FLOW(2).X the 2-ph  region. For a flow with a liquid film, FLOW(1).X
%  is for the vapor region, FLOW(2).X for the liquid flow within the
%  membrane and FLOW(3).X for the liquid film in front of the membrane.
%
%  Calls CPL, INT2PH, INTG, INTP, K, PCAP, PS, Q_M, R, TEMP, TS, X.
%  Called from MBACK.
%  See also FLOWSTRUCT.

% disp(sprintf('m=%.9g',m)); %debug
info = struct('kap',kappa,'m',m,'Ca',[],'kapc',[],'kapf',[],'L',L,...
  'T0',[],'p0',[],'dp',[],'ph',fmodel);
sol = struct('len',{},'a3',{},'q3',{},'T0',{},'Te',{},'p0',{},'pe',{},...
  'de',{},'df',{});

% flow struct count
i=1; 

sol(1).Te=Te; % only the first struct-assignment needs an index
sol.pe=pe;
zrange = [L 0];

% start to integrate from the end
%----------------------------- vapor flow
fg = intg(m,Te,pe,0,zrange);
%fg.xe %debug
%fg.ye %debug

if isempty(fg.ie)
  %disp('got through') %debug
  %------------------- vapor flow all through
  % got through, write this part of the flow solution
  % construct the result struct
  flow = struct('z',{fg.x},...
    'T',{fg.y(1,:)},...
    'p',{fg.y(2,:)},...
    'q',{fg.y(3,:)},...
    'a',{ones(size(fg.x))},...
    'x',{ones(size(fg.x))},...
    'color',{'r'});

  % 2 - z=0
  T2 = fg.y(1,end);
  p2 = fg.y(2,end);
  ps2 = ps(T2);

  if p2-curv*sig(T2)<ps2
    %----------------- vapor flow in front of the membrane
    %disp('Vapor flow in front of the membrane not supported.')
    % but write the solution
    sol.len = -i;
    sol.a3 = 1;
    sol.T0 = fg.y(1,end);
    sol.p0 = fg.y(2,end);
    sol.q3 = fg.y(3,end);
    sol.df = 0;
    flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return
  end
  % p2-curv*sig(T2)>=ps2
  %--------------- vapor and liquid film
  % calculate the film thickness by solving heat conduction eq, see TEMP
  % it is assumed z2=0
  % evaporation enthalpy is increased!
  q2 = fg.y(3,end)+m*rkelv(T2); % + dhdp(Te)*(ps2-p2);
  c12 = cpl(T2);
  k12 = kl(T2);
  T1 = T2;
  zold = 0;
  z1 = L;
  while 1-zold/z1>0.001
    zold = z1;
    c12 = 0.5*(cpl(T1)+cpl(T2));
    k12 = 0.5*(kl(T1)+kl(T2));
    z1 = k12/(m*c12)*log(m*r(T1)/q2);
    T1 = temp(z1,m,T2,q2,0,c12,k12);
  end
  % and the result is
  z01 = [0 z1/2 z1];
  [T01,q01] = temp(z01,m,T2,q2,0,c12,k12);

  % add liquid film to the result struct
  i = i + 1;
  flow(i).z = z01;
  flow(i).T = T01;
  flow(i).p = ps(T1)*ones(1,3);
  flow(i).q = q01;
  flow(i).a = zeros(1,3);
  flow(i).x = zeros(1,3);
  flow(i).color = 'b';

  % write the solution
  sol.len = -i;
  sol.a3 = 0;
  sol.T0 = T01(3);
  sol.p0 = flow(i).p(3);
  sol.q3 = q01(3);
  sol.de = 0;
  sol.df = z1;
  flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return 

else % ~isempty(fg.ie)
  %--------- vapor flow to evaporation front
  % position of the evaporation front
  zvap = fg.xe(end);
  yvap = fg.ye(:,end);

  % construct the result struct
  flow = struct('z',{[fg.x(1:end-1) zvap]},...
    'T',{[fg.y(1,1:end-1) yvap(1)]},...
    'p',{[fg.y(2,1:end-1) yvap(2)]},...
    'q',{[fg.y(3,1:end-1) yvap(3)]},...
    'a',{ones(size(fg.x))},...
    'x',{ones(size(fg.x))},...
    'color',{'r'});

end %vapor flow

%----------------------------- 2ph flow
try
  %--------- 2ph flow all through
  % continue to integrate 2-ph (if possible!)
  % (check of crit. heat flux in int2ph - therefore 'try'
  f2ph=int2ph(m,yvap(1),1,yvap(3),[zvap 0]);

  % add to the result struct
  i = i + 1;
  flow(i).z = f2ph.x;
  flow(i).T = f2ph.y(1,:);
  flow(i).p = kelv(f2ph.y(1,:)).*ps(f2ph.y(1,:));
  flow(i).q = m*q_m(f2ph.y(1,:),f2ph.y(2,:));
  flow(i).a = f2ph.y(2,:);
  flow(i).x = x(f2ph.y(1,:),f2ph.y(2,:));
  flow(i).color = 'g';

  % and add a liquid film for a nonwetting fluid, costheta<=0
  if ps(flow(i).T(end))<flow(i).p(end)
    %disp('liq. film!')
    % not necessary
    if costheta<=0

  p2 = flow(i).p(end);
  T2 = flow(i).T(end);
  q2 = flow(i).q(end)+xdot(T2,flow(i).a(end))*m*rkelv(T2);
      
  %--------- liquid film

  T1 = Ts(p2);
  % mean material properties
  c12 = 0.5*(cpl(T1)+cpl(T2)); % mean
  k12 = 0.5*(kl(T1)+kl(T2));
  % calculate the film thickness (for z2=0)
  z1 = k12/(m*c12)*log(1-(T1-T2)*c12*m/q2);
  z01 = [0 z1/2 z1];
  [T01,q01] = temp(z01,m,T2,q2,0,c12,k12);
  % strange test
  if ( (T01(end)-T1)/T1>1e-6 ) | z1>0 
    warning(sprintf(['calculating film thickness and/or T1\n'...
    '  z1/L = %g, (T-T1)/T1 = %g'],z1/L,(T01(end)-T1)/T1));
  end
  % film part of the result struct
  %z01=z1/2;
  %[t01,q01] = temp(z01,m,T2,q2,0,c12,k12);

  % add liquid film to the result struct
  i = i + 1;
  flow(i).z = z01;
  flow(i).T = T01;
  flow(i).p = [p2 p2 p2];
  flow(i).q = q01;
  flow(i).a = [0 0 0];
  flow(i).x = [0 0 0];
  flow(i).color = 'b';

  % write the solution
  sol.len = -i;
  sol.a3 = 0;
  sol.T0 = T01(end);
  sol.p0 = p2;
  sol.q3 = q01(end);
  sol.de = zvap;
  sol.df = z1;
  flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return
    else % costheta<0
      %during iteration, that may happen
      disp(sprintf('strange, ps(T3)<p3: ps(T3)-p3 = %g',...
        ps(flow(i).T(end))-flow(i).p(end)));
    end
  end
  
  % write the solution
  sol.len = -i;
  sol.a3 = f2ph.y(2,end);
  sol.T0 = f2ph.y(1,end);
  sol.p0 = flow(i).p(end); % not ps(f2ph.y(1,end));
  sol.q3 = m*q_m(f2ph.y(1,end),f2ph.y(2,end));
  sol.de = zvap;
  sol.df = 0;
  flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return
%catch

end %2ph-flow

%----------------------------- liquid flow
% defs: 3 - evaporation front , 2 - z=0, 1 - film surface

% 3 -- 2
T3 = yvap(1);
q3 = yvap(3)+m*rkelv(yvap(1)); % jump condition!
p3 = yvap(2)-curv*sig(T3); % subtract cap. pressure
z3 = zvap;

% calculate the temperature distribution
% initially, properties at z3
T2 = temp(0,m,T3,q3,z3,cpl(T3),k(T3,0));
% then, use mean material properties
c23 = 0.5*(cpl(T3)+cpl(T2)); % mean
k23 = 0.5*(k(T3,0)+k(T2,0));

% compute pressure distribution
fp = intp(m,T3,p3,q3,[z3 0],c23,k23);

if isempty(fp.ie)
  %--------- liquid flow all through
  % temperature distribution in the liquid flow
  [tp,qp] = temp(fp.x,m,T3,q3,z3,c23,k23);
  % add liquid flow to the result struct
  i = i + 1;
  flow(i).z = fp.x;
  flow(i).T = tp;
  flow(i).p = fp.y(1,:);
  flow(i).q = qp;
  flow(i).a = zeros(size(fp.x));
  flow(i).x = zeros(size(fp.x));
  flow(i).color = 'b';

  p2 = fp.y(1,end);
  q2 = qp(end);

  if p2<=ps(T2)
    %--------- just a meniscus

    % write the solution
    sol.len = -i;
    sol.a3 = 0;
    sol.T0 = T2;
    sol.p0 = ps(T2);
    sol.q3 = qp(end);
    sol.de = zvap;
    sol.df = 0;
    flowstruct = struct('info',info,'sol',sol,'flow',flow);
    return

  end
  % else % p2>ps(T2)
  %--------- liquid film

  T1 = Ts(p2);
  % mean material properties
  c12 = 0.5*(cpl(T1)+cpl(T2)); % mean
  k12 = 0.5*(kl(T1)+kl(T2));
  % calculate the film thickness (for z2=0)
  z1 = k12/(m*c12)*log(1-(T1-T2)*c12*m/q2);
  z01 = [0 z1/2 z1];
  [T01,q01] = temp(z01,m,T2,q2,0,c12,k12);
  % strange test
  if ( (T01(end)-T1)/T1>1e-6 ) | z1>0 
    warning(sprintf(['calculating film thickness and/or T1\n'...
    '  z1/L = %g, (T-T1)/T1 = %g'],z1/L,(T01(end)-T1)/T1));
  end
  % film part of the result struct
  %z01=z1/2;
  %[t01,q01] = temp(z01,m,T2,q2,0,c12,k12);

  % add liquid film to the result struct
  i = i + 1;
  flow(i).z = z01;
  flow(i).T = T01;
  flow(i).p = [p2 p2 p2];
  flow(i).q = q01;
  flow(i).a = [0 0 0];
  flow(i).x = [0 0 0];
  flow(i).color = 'b';

  % write the solution
  sol.len = -i;
  sol.a3 = 0;
  sol.T0 = T01(end);
  sol.p0 = p2;
  sol.q3 = q01(end);
  sol.de = zvap;
  sol.df = z1;
  flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return

else
  %--------- liquid flow with 2ph flow
  % condensation front within the membrane
  % def: 2 - condensation front
  z2 = fp.xe(end);
  % 2nd pass
  % recalculate temperature distribution
  T2 = temp(z2,m,T3,q3,z3,c23,k23);
  c23 = 0.5*(cpl(T3)+cpl(T2));
  k23 = 0.5*(k(T3,0)+k(T2,0));
  fp = intp(m,T3,p3,q3,[z3 0],c23,k23);
  if isempty(fp.ie)
    error(['This error should not occur!\n First there is a'...
      'condensation front, then it vanishes.']);
  end
  z2 = fp.xe(end);
  p2 = fp.ye(1,end); % do not add cap. pressure here; done in INTP
  ind = find(fp.x>fp.xe(end));
  zp = [fp.x(ind) z2];
  [tp,qp] = temp(zp,m,T3,q3,z3,c23,k23);
  T2 = tp(end);
  q2 = qp(end);

  % add liquid flow to the result struct
  i = i + 1;
  %zp
  %zp(1) - zp(2)
  flow(i).z = zp;
  flow(i).T = tp;
  flow(i).p = [fp.y(1,ind) p2];
  flow(i).q = qp;
  flow(i).a = zeros(size(zp));
  flow(i).x = zeros(size(zp));
  flow(i).color = 'b';

  % now 2ph-flow should be possible
  try
    %--------- 2ph flow after liq. flow
    f2ph=int2ph(m,T2,0,q2,[z2 0]);
    if ~isempty(f2ph.ie)
      error('Now this is funny!');
    end
  catch
    disp('I increase the heat flux by 6%.');
    f2ph=int2ph(m,T2,0,1.01*q2,[z2 0]);
    %error('Now what?');
  end

  % add to the result struct
  i = i + 1;
  flow(i).z = f2ph.x;
  flow(i).T = f2ph.y(1,:);
  flow(i).p = kelv(f2ph.y(1,:)).*ps(f2ph.y(1,:));
  flow(i).q = m*q_m(f2ph.y(1,:),f2ph.y(2,:));
  flow(i).a = f2ph.y(2,:);
  flow(i).x = x(f2ph.y(1,:),f2ph.y(2,:));
  flow(i).color = 'g';

  % and add a liquid film for a nonwetting fluid, costheta<=0
  if ps(flow(i).T(end))<flow(i).p(end)
    %disp('liq. film!')
    % not necessary
    if costheta<=0

  p2 = flow(i).p(end);
  T2 = flow(i).T(end);
  q2 = flow(i).q(end)+xdot(T2,flow(i).a(end))*m*rkelv(T2);
      
  %--------- liquid film

  T1 = Ts(p2);
  % mean material properties
  c12 = 0.5*(cpl(T1)+cpl(T2)); % mean
  k12 = 0.5*(kl(T1)+kl(T2));
  % calculate the film thickness (for z2=0)
  z1 = k12/(m*c12)*log(1-(T1-T2)*c12*m/q2);
  z01 = [0 z1/2 z1];
  [T01,q01] = temp(z01,m,T2,q2,0,c12,k12);
  % strange test
  if ( (T01(end)-T1)/T1>1e-6 ) | z1>0 
    warning(sprintf(['calculating film thickness and/or T1\n'...
    '  z1/L = %g, (T-T1)/T1 = %g'],z1/L,(T01(end)-T1)/T1));
  end
  % film part of the result struct
  %z01=z1/2;
  %[t01,q01] = temp(z01,m,T2,q2,0,c12,k12);

  % add liquid film to the result struct
  i = i + 1;
  flow(i).z = z01;
  flow(i).T = T01;
  flow(i).p = [p2 p2 p2];
  flow(i).q = q01;
  flow(i).a = [0 0 0];
  flow(i).x = [0 0 0];
  flow(i).color = 'b';

  % write the solution
  sol.len = -i;
  sol.a3 = 0;
  sol.T0 = T01(end);
  sol.p0 = p2;
  sol.q3 = q01(end);
  sol.de = zvap;
  sol.df = z1;
  flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return
    else % costheta<0
      %during iteration, that may happen
      disp(sprintf('strange, ps(T3)<p3: ps(T3)-p3 = %g',...
        ps(flow(i).T(end))-flow(i).p(end)));
    end
  end
  

  % write the solution
  sol.len = -i;
  sol.a3 = f2ph.y(2,end);
  sol.T0 = f2ph.y(1,end);
  sol.p0 = flow(i).p(end); % not ps(f2ph.y(1,end));
  sol.q3 = m*q_m(f2ph.y(1,end),f2ph.y(2,end));
  sol.de = zvap;
  sol.df = 0;
  flowstruct = struct('info',info,'sol',sol,'flow',flow);
  return

end %~isempty(fp.ie)
