function fj = join(fs,fl)
%JOIN       Join two flowstructs.
%  JOIN(FS,FL) joins the two flowstructs or flowstruct arrays FS and FL.
%  Returns a flowstruct where the .INFO and .SOL field is filled from FS
%  and the .LIN field is filled from FL.
%
%  See also FLOWSTRUCT.

fj = fl;
n = length(fl);
if n ~= length(fs)
  error('The two flowstructs are of unequal dimension.');
end
for i = 1:n
  m = length(fs(i).flow);
  [fj(i).info fj(i).sol fj(i).flow(end+1:end+m)]...
    = deal(fs(i).info,fs(i).sol,fs(i).flow(:));
end
