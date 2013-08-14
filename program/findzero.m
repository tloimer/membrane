function b = findzero(presiduum,mguess,p1tol,solver)
%FINDZERO  Find a solution to PRESIDUUM(M,SOLVER) = 0.

%  INITIALIZE
% Information is plotted for trace > 1. Iteration results for trace > 2.
trace = 1;
global VERBOSE; if VERBOSE > 0, trace = trace + VERBOSE; end
fcount = 0; procedure = ' ';
savesolver = solver;

%  ASSIGN CRUDE SOLVER TOLERANCES
% During iteration, use crude tolerances.
solver = solverstruct('crude');
% defaults: solver.writesolution = false;
%	    solver.fullsolution = false; solver.partialsolution = true;

a(1) = mguess; fa(1) = presiduum(a(1),solver); fcount = 1;

%  SEARCH AN INTERVAL
% Now look, in which direction we have to search for pres =  0.
if fa(1) > 0, fac = 0.5; else fac = 1.4; end

for i = 2:6
  old = i - 1;
  a(i) = a(old)*fac;
  fa(i) = presiduum(a(i),solver);
  fcount = fcount + 1;
  if (fa(old) < 0) == (fa(i) > 0) || fa(i) == 0
    % found an interval.
    a = [a(i) a(old)];
    fa = [fa(i) fa(old)];
%DEBUG disp(sprintf('Intervall: %8f < m < %8f, %8f < pres < %8f',...
%DEBUG   a(1),a(2),fa(1),fa(2)));
    break
  end
end

%  PLOT SEARCH, IF NO INTERVAL FOUND
if size(a,2) ~= 2 %DEBUG || true(); %DEBUG
  % No interval found.
  % also plot the result for mvap and mliq
  %mvap = p12*mem.kappa/(flsetup.nuapp(T2,(p1+p2)/2)*mem.L);
  p%vap = presiduum(mvap,solver);
  % No, really no interval found.
  disp(sprintf('MNUM: No interval found! mguess = %g',mguess));
  %disp('  Call MNUM with a guess for the mass flow,');
  %disp('  mnum(T1,p1,p2,theta,s,mem,f,''m'',mguess)');
  disp('  Mass flow:'); disp(a);
  disp('  p_res :'); disp(fa);
  % reuse fac as figure handle
  fac = figure('Name','Search for Interval','NumberTitle','off');
  plot(a,fa/p1,'ko-'); %,mvap,pvap/p1,'r*'); %,mliq,pliq/p1,'b*');
  xlabel('mflux'); ylabel('p_{1,calc}/p_1 - 1');
  % AND ASK FOR AN INTERVAL
  while size(a,2)~= 2
    b = input(['  Provide an interval [mflux mflux] or try\n' ...
      '  several points [mflux mflux mflux ...]: ']);
    lenb = size(b,2); fb = b;
    for i = 1:lenb,  fb(i) = presiduum(b(i),solver);  end
    % reuse old as figure handle; save a possible different current figure.
    old = get(0,'CurrentFigure'); figure(fac);
    line(b,fb/p1,'Marker','+','LineStyle','none');
    figure(old);
    if lenb == 2 && (fb(1) < 0) == (fb(2) > 0)
      a = b; fa = fb;
    end
  end
  close(fac);
end

%  INTERVAL FOUND OR ASSIGNED
% check for zero at boundaries
b = a(2); a = a(1); fb = fa(2); fa = fa(1);
lowprecision = true;

% if a or b == 0, the while-loop is run once, lowprecision set to false, and the
% result re-calculated.

tol = eps; % we do not care for accuracy in m
fc = fb;
procedure = 'initial';
header2 = ' Func-count    x          f(x)             Procedure';
if trace > 2
    disp(header2)
end
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end

    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(fb) < p1tol) || (abs(m) <= toler) % || (fb == 0.0)
        if lowprecision
	  if trace > 1
	    disp(sprintf(['p1calc - p1 = %g Pa, m = %g kg/m2s. We jumped '...
			  'over 100xp1tol.'], fb, b));
	  end
	else
          break
	end
    end
    if trace > 2
        disp(sprintf('%5.0f   %13.6g %13.6g        %s',fcount, b, fb, procedure));
    end

    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
        procedure='bisection';
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, q = -q; else p = -p; end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
            procedure='interpolation';
        else
            d = m;  e = m;
            procedure='bisection';
        end;
    end % Interpolation

    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, b = b + d;
    elseif b > c, b = b - toler;
    else b = b + toler;
    end
    if lowprecision && abs(fb) < 100*p1tol
      % Continue with the demanded solver tolerance
      solver = savesolver;
      lowprecision = false;
      if trace > 1
	disp(sprintf(['p1calc - p1 = %g Pa, m = %g kg/m2s. Switched to '...
		      'required precision.'], fb, a));
      end
    end
    fb = presiduum(b,solver);
    fcount = fcount + 1;
end % Main loop

% Output last chosen b
if trace > 2
    disp(sprintf('%5.0f   %13.6g %13.6g        %s',fcount, b, fb, procedure));
end

%%% END FINDZERO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END FINDZERO %%%
