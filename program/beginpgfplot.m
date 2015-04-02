function pid = beginpgfplot(pname,options)
%BEGINPGFPLOT Write the begin of a pgfplot.
%  BEGINPGFPLOT(FILENAME) writes the pgfplot FILENAME.PGFPLOT.
%  BEGINPGFPLOT(FILENAME,OPTIONS) adds the string OPTIONS to the axis.
%  OPTIONS is parsed by fprint, hence \\ gives \, \n is a newline,
%  %% yields %.
%
%  See also ADDCOORDS, ENDPGFPLOT.

if nargin==1, options=''; end
pid = fopen([pname '.pgfplot'],'w');
fprintf(pid, ['%%\\usepackage{pgfplots} \\pgfplotsset{compat=newest}\n'...
 '%%\\pgfplotsset{/tikz/font=\\small, %% kleinere Beschriftungen\n'...
 '%%  /pgfplots/every axis y label/.append style={rotate=-90,xshift=10pt}\n'...
 '%%  \\addplot[mark=*,mark options={scale=0.6}], legend style={draw=white}\n'...
 '\\begin{tikzpicture}[baseline]\n\\begin{axis}[' options ']\n']);
% '%% xlabel= , xmax= , xmin=, ylabel=,\n'],'w');
