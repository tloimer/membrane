function flowstruct = mback(deltap,T0,L)
%MBACK      Backward integration.
%  MBACK(DELTAP,T0,L) calculates the solution for the flow through a
%  membrane for given pressure difference DELTAP, initial temperature T0
%  and thickness of the membrane L. Returns the solution in a struct
%  FLOWSTRUCT. Non-unique solutions are returned in an FLOWSTRUCT array.
%
%  See FLOWSTRUCT.
%
%  Calls INTG, JT, MLIN, PS, FLOWBACK.

% input check
p0=ps(T0);
if (deltap>p0)
  error( sprintf('deltap must be smaller than ps(T0) = %g',p0) );
  return
end

n = 1; % flowstruct counter, for multiple solutions
tolm = 1e-4; % tolerance in m
pe = p0 - deltap;
Te = T0 - jt(T0,p0,deltap);
kap = kappa;
kapc = kappac(T0);

if costheta<0
  % ------------------------------------------------- non-wetting liquid
  % liquid film possible?
  if -curv*sig(Te)>p0-ps(Te)
    % most probably a liquid film; now check exactly
    mb = kap*(ps(Te)-pe)*2/((nug(T0,p0)+nug(Te,pe))*L);
    zrange = [L 0];
    fg = intg(mb,Te,pe,0,zrange);
    while ~isempty(fg.ie)
      mb = mb*0.9;
      fg = intg(mb,Te,pe,0,zrange);
    end
    % now mb is smaller than m
    fac = (ps(fg.y(1,end))-pe)/(fg.y(2,end)-pe);
    while fac>1+tolm
      mb = mb*fac;
      fg = intg(mb,Te,pe,0,zrange);
      %fac = (ps(T2)-pe)/(p2-pe); % here, 2 - z=0.
      fac = (ps(fg.y(1,end))-pe)/(fg.y(2,end)-pe);
      % really, fac never becomes < 1?
      % Comment that out, after suitable testing.
      if fac<1
        error('The approach to p(z=0) = ps(T(z=0)) failed.')
      end
    end
    m = mb;
    % now we have calculated p(z=0) and m

    T2 = fg.y(1,end);
    q2 = fg.y(3,end)+m*r(T2);
    % liquid film really possible?
    if -curv*sig(T2)>p0-ps(T2)
      % ---------------------------------- liq.film - vapor, non-wetting
      % calculate the thickness of the liquid film
      % mean material properties
      c02 = 0.5*(cpl(T0)+cpl(T2)); % mean
      k02 = 0.5*(kl(T0)+kl(T2));
      % calculate the film thickness (for z2=0)
      z0 = k02/(m*c02)*log(1-(T0-T2)*c02*m/q2);
      z02 = [0 z0/2 z0];
      [T02,q02] = temp(z02,m,T2,q2,0,c02,k02);

      % and write the result
      flowstruct = struct('info',struct('kap',kap,'kapc',kapc,...
        'Ca',jt(T0,p0)*dpsdT(T0)*deltap/(curv*sig(T0)),'kapf',[],'m',m,...
	'L',L,'T0',T0,'p0',p0,'dp',deltap,'ph',fmodel),...
	'sol',struct('len',-2,'a3',0,'q3',q02(end),'T0',T0,'Te',Te,'p0',p0,...
	'pe',pe,'de',0,'df',z0),'flow',struct('z',{fg.x,z02},...
	'T',{fg.y(1,:),T02},'p',{fg.y(2,:),p0*ones(1,3)},'q',{fg.y(3,:),q02},...
	'a',{ones(size(fg.x)),zeros(1,3)},'x',{ones(size(fg.x)),zeros(1,3)},...
	'color',{'r','b'}));

      % since there is a liquid film, the range of m to look for a
      % second solution should be limited m larger than m(liq. film).
      % overwrite mguess
      n = n+1;
      mguess = 1.2*mb;
      % in case only the film solution exists, it will be returned a
      % second time, approximately identical
      disp('film solution exists');
    end
  end
end
      
% the linear theory provides us with an estimate for the mass flux
flin = mlin(deltap,T0,L);
% which is used only if mguess was not given above
if ~exist('mguess','var')
  mguess = flin(1).lin.m;
end

% find an interval for m
% initially, mguess should be too small
warning off;
pres=m2ph(mguess,Te,pe,L,p0);
if pres>0
  fac=1.2; % i carefully extend the range!
  mult=1;
else
  fac=1/1.2;
  mult=-1;
end
mb = mguess;
while mult*pres>0
  mb = fac*mb;
  pres = m2ph(mb,Te,pe,L,p0);
end

%disp(' - - start iteration') %debug
% Display: none iter final notify
options =optimset(optimset('fzero'),'TolX',mb*tolm,'Display','none');
m = fzero(@m2ph,[mb/fac mb],options,Te,pe,L,p0);
% warning backtrace; % flowback is called again anyway
%disp(' - - last call to FLOWBACK') %debug
warning on;

% only write the flowstruct, if the solution is different from before
if n==1 | abs(1-m/flowstruct(1).info.m)>tolm
  flowstruct(n)=flowback(m,Te,pe,L);
  flowstruct(n).info.T0 = T0;
  flowstruct(n).info.p0 = p0;
  flowstruct(n).info.dp = deltap;
  flowstruct(n).info.kapc = kapc;
  if costheta~=0
    flowstruct(n).info.Ca = dpsdT(T0)*jt(T0,p0)*deltap/(curv*sig(T0));
    if costheta>0
      flowstruct(n).info.kapf = kapc + kap*(1-sqrt(1+4*...
        flowstruct(n).info.Ca^2*kapc/kap))/(2*flowstruct(n).info.Ca^2);
    end
  end
end

% add linear theory values
if length(flin)~=n % length(flowstruct)
  warning(...
  'linear and numerical calculations give different number of solutions.');
  return
end

[flowstruct(:).lin] = deal(flin(:).lin);

%-----------------------------------------------------------------------
function pres = m2ph(m,Te,pe,L,p0)
flowstruct=flowback(m,Te,pe,L);
pres = p0 - flowstruct.sol.p0;
