function pmout = phmodel(model)
%PHMODEL    Set or query the 2-ph flow model.
%  PHMODEL returns a string denoting the 2ph-model currently in use.
%  
%  PHMODEL(MODEL) sets the 2ph-model. MODEL can be either 'plug',
%  'bundle' or 'annular'.

if (nargin==0)
  pmout = fmodel;
else
  arg = [model '/nu.m ' model '/xdot.m ' model '/fmodel.m .'];
  stat = unix(['ln -sf ' arg]);
end
