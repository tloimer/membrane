function h = memplot(flowstruct,fieldname)
%MEMPLOT    Plot one variable.
%  MEMPLOT(FLOWSTRUCT,VAR) plots the the value of a variable VAR
%  contained in the structure FLOWSTRUCT across the membrane.
%
%  Called from RESPLOT.

com = 'plot(';
for i=1:length(flowstruct)
  ffi=sprintf('flowstruct(%d).',i);
  com = [com ffi 'z,' ffi fieldname ',' ffi 'color,'];
end
com(end)=')';
eval(com);
