function solver = solverstruct(solvertype) %----------------------- solverstruct
%SOLVERSTRUCT Construct a solver struct.
%  SOLVERSETUP returns a struct SOLVER set up for crude tolerances
%
%  SOLVERSETUP('accurate') sets accurate solver tolerances.
%
%  By default, SOLVER.writesolution = false and SOLVER.fullsolution = false.
%
%  SOLVER.partialsolution is obsolete, but still supported.
%  SOLVER.partialsolution = ~SOLVER.fullsolution.
%
%  See also ASYM, MSTACK>SOLVERSETUP.

solver = struct('rtol',1e-3,'atol',1e-5,'tola',1e-6,'maxTperstep',2,...
  'maxpperstep',1e5,'odemaxstep',[],...
  'writesolution',false,'fullsolution',false,'partialsolution',true); %%CHANGED!
%  'maxpperstep',1e5,'writesolution',false,'partialsolution',true); %%CHANGED!
if strcmp(solvertype,'accurate')
  solver.rtol = 1e-4;  solver.atol = 1e-6;
  solver.maxTperstep = .4;  solver.maxpperstep = 4e4;
end

% Step size calculation for ode45, ode23t.
% minsteps (= 2): minimum number of integration steps;
% odemaxstep: number of (power of 2: 2, 4, 8, ...) steps such that T or p does
% not change more than maxTperstep (maxpperstep) within one step.
solver.odemaxstep = ... % 1/max(minsteps,...
  @(range,maxstep) 1/max(2,pow2(ceil(log2(abs(range)/maxstep))));


end %---------------------------------------------------------- end solverstruct
