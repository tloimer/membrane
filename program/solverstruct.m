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
  'maxpperstep',1e5,...
  'writesolution',false,'fullsolution',false,'partialsolution',true); %%CHANGED!
%  'maxpperstep',1e5,'writesolution',false,'partialsolution',true); %%CHANGED!
if strcmp(solvertype,'accurate')
  solver.rtol = 1e-4;  solver.atol = 1e-6;
  solver.maxTperstep = .4;  solver.maxpperstep = 4e4;
end

end %---------------------------------------------------------- end solverstruct
