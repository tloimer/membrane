function fi = addint(fs)
%ADDINT     Add numerically integrated solution.
%  ADDLIN(FS) Add the solution from linear theory to the flowstruct FS.
%  ADDINT returns a new flowstruct and leaves the old one unchanged.

%we write a new flowstruct
fi = fs;

% only continue if we really have a bit of liquid flow
if fi.flow(2).a~=0
  return
end

% values given
m = fi.info.m;
L = fi.info.L;
T0 = fi.info.T0;
p0 = fi.info.p0;
%deltap = fi.info.dp;

% for function calls that depend on the 2ph-model
if ~strcmp(fmodel,fi.info.ph)
  oldmodel = fmodel;
  phmodel(fi.info.ph);
end

i = length(fi.flow) + 1;

% integrate the liquid flow within the membrane
%z0 = fi.sol.de;
Ti = fi.flow(1).T(end);
pi = fi.flow(1).p(end);
qi = fi.flow(1).q(end) + m*r(Ti);

%[fi.sol.de 0]

fl = intl(m,Ti,pi,qi,[fi.sol.de 0]);

% add this part of the result
fi.flow(i).z = fl.x;
fi.flow(i).T = fl.y(1,:);
fi.flow(i).p = fl.y(2,:);
fi.flow(i).q = fl.y(3,:);
fi.flow(i).a = zeros(size(fl.x));
fi.flow(i).x = zeros(size(fl.x));
fi.flow(i).color = 'k:';

i = i+1;

%if ~isempty(fl.ie)
%  % we have a liquid film of zero thickness
%  warning(sprintf(['integration of liquid flow stopped at '...
%    'z/L = %.3g,\n  T-T0 = %.3gK, ps(T)-p = %.3gPa, p-p0 = %.3gPa.'],...
%    fl.x(end)/L, fl.y(1,end)-T0, ps(fl.y(1,end))-fl.y(2,end), fl.y(2,end)-p0));
%else
Ti = fl.y(1,end);
pi = fl.y(2,end);
qi = fl.y(3,end);
if ps(Ti)>pi
  % end
  return
else
  %disp(sprintf('liquid flow, z=%.3g: T = %.3gK, ps(T)-p = %.3gPa',...
    %fl.x(end)/L, Ti, ps(Ti)-pi));
  %return
  % how far to integrate?
  if fi.sol.df<0
    zf = 2*fi.sol.df;
  else
    zf = -2*(Ts(pi)-Ti)*kl(Ti)/qi;
  end
  fliq=struct('ie',[]);
  n = 0;
  while isempty(fliq.ie)
    fliq = intliq(m,Ti,pi,qi,[0 zf]);
    n = n+1;
    zf=2*zf;
    if n>9
      error('too long while-loop');
    end
  end %while

  % add this part of the result
  fi.flow(i).z = fliq.x;
  fi.flow(i).T = fliq.y(1,:);
  fi.flow(i).p = fliq.y(2,:);
  fi.flow(i).q = fliq.y(3,:);
  fi.flow(i).a = zeros(size(fliq.x));
  fi.flow(i).x = zeros(size(fliq.x));
  fi.flow(i).color = 'k:';

end %if ps>p

if ( exist('oldmodel','var') )
  phmodel(oldmodel);
end
