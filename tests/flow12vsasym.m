function flow12vsasym(T1,p1,p2,theta,s,mem,f,accuracy)
%FLOW12VSASYM Plot results from flow12.m and asym.m, for comparison.
%  FLOW12VSASYM(T1,P1,P2,THETA,SUBSTANCE,MEMBRANE,FMODEL,ACCURACY) ACCURACY is
%  optional, default 'accurate', can be set to 'crude'.

if nargin < 8
  accuracy = 'accurate';
end

% First run mnum
% Although one could provide the solution from asym to mnum as mguess, asym did
% produces errors, which does not help at all.
[mflow,fl] = mnum(T1,p1,p2,theta,s,mem,f); %'m',masym);
try % 2-phase flow does not have q
  pTzplots(fl,'flow12',false,true); %,535 - is the default vertical position
catch
end

% Then run asym
ms = mstackstruct(theta,mem,f);
% a too precise guess probably wrecked havoc,
%ms.mguess = mflow; % Changed mnumadiabat, can deal with ms.mguess.
[masym,ms] = mnumadiabat(T1,p1,p2,s,ms,accuracy);

% Print the solution
ms.printsolution(ms);

fprintf('mnum %.3g g/m2s, mnumadiabat %.3g g/m2s; mnumadiabat %+.0f ppm\n',...
	mflow*1e3, masym*1e3, (masym/mflow-1)*1e6);

% Plot
flasym = ms.singlemstofl(ms);
pTzplots(flasym,'asym',false,true,25);
