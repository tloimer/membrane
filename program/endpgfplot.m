function endpgfplot(pid,legend)
%ENDPGFPLOT  Close a pgfplot.
%  ENDPGFPLOT(PID) closes the pgfplot with file identifier PID.
%  ENDPGFPLOT(PID,LEGEND) adds the legend LEGEND, a string where commas separate
%  the legend identifiers: ENDPGFPLOT(PID,'data 1, data2').

if nargin==2
  fprintf(pid,['\\legend{' legend '}\n']);
else
  fprintf(pid,'%%\\legend{data1, data2}\n');
end
fprintf(pid,'\\end{axis}\\end{tikzpicture}\n');
fclose(pid);
