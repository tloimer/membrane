function addcoords(pid,X,Y,options)

if nargin==3
  options = [];
end

if size(X,1) == 1
  X = X';
end
if size(Y,1) == 1
  Y = Y';
end

% the accuracy; think, e.g., of  a range between 292.1 and 292.5 K
%xmin = min(X);
%ymin = min(Y);
%xrange = max(X) - xmin;
%yrange = max(Y) - ymin;

%sX = int2str(floor(log10(xmin))-floor(log10(xrange))+3);
%sY = int2str(floor(log10(xmin))-floor(log10(xrange))+3);
%sformat = sprintf(' (%%.%sg,%%.%sg)',sX,sY);
%sformat = ['\n' sformat sformat sformat sformat sformat];

fprintf(pid, ['\\addplot[black,' options '] coordinates {']);
%fprintf(pid, sformat, [X';Y']);
% the format string, if the above is too obscure
fprintf(pid, ...
 '\n (%.5g,%.5g) (%.5g,%.5g) (%.5g,%.5g) (%.5g,%.5g) (%.5g,%.5g)',...
 [X';Y']);
% Die letzten zwei Zeichen, ' (', löschen, falls nötig.
if mod(size(X,1),5) ~=0, fseek(pid,-2,0); end
fprintf(pid,'\n};\n');
