function addcoords(pid,X,Y,options)

fprintf(pid,...
  ['\\addplot[black,' options '] coordinates {']);
fprintf(pid,...
 '\n (%.4g,%.4g) (%.4g,%.4g) (%.4g,%.4g) (%.4g,%.4g) (%.4g,%.4g)',...
 [X';Y']);
% Die letzten zwei Zeichen, ' (', löschen, falls nötig.
if mod(size(X,1),5) ~=0, fseek(pid,-2,0); end
fprintf(pid,'\n};\n');
