function [mflow,mcrude,maccurate,pflow,pcrude,pacc] = flow12andasym(T1,p1,p2,theta,s,mem,f)
%FLOW12ANDASYM Compute results for flow12.m and asym.m with crude and accurate tolerances.

% First run mnum
% Although one could provide the solution from asym to mnum as mguess, asym did
% produces errors, which does not help at all.
[mflow,fl] = mnum(T1,p1,p2,theta,s,mem,f); %'m',masym);
pflow = fl.sol.p1 - fl.info.p1;

% Then run asym
ms = mstackstruct(theta,mem,f);
ms.mguess = mflow;
[mcrude,ms] = mnumadiabat(T1,p1,p2,s,ms,'crude');
pcrude = ms.p1sol - ms.p1in;

ms = mstackstruct(theta,mem,f);
ms.mguess = mflow;
[maccurate,ms] = mnumadiabat(T1,p1,p2,s,ms);
pacc = ms.p1sol - ms.p1in;
