%RUNALLTESTS A simple reminder to run all tests. Stay in the directory tests/.

[path,name] = fileparts(cd);
if ~strcmp(name,'tests'), error('Run from directory tests/.'); end
clear('path','name');

if ~exist('substance.m','file'), addpath('../program'); end

global VERBOSE;  VERBOSE = 1;
if ~exist('INTERACTIVE','var')
  INTERACTIVE = true;
end

entropy
%fprintf('To remove figures, type "close all". To resume, type "return".\n');
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf(['Try VERBOSE = 0; to suppress messages.\n'...
	   'Clear INTERACTIVE or try INTERACTIVE = false to batch run all '...
	   'tests.\nTo resume, type "return" (literally, six keys!).\n']);
  keyboard;
end
close all; clear all

% loads ../matlab3/results1105,
% read1103 - reads ../data/membranes11.tsv, ../data/data0905.tsv
frvmCcc
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf('To resume, type "return".\n');
  keyboard;
end

% needs `load' and `read' from before
frvmvskap
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf('To resume, type "return".\n');
  keyboard;
end
close all; clear all

testmstack
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf('To resume, type "return".\n');
  keyboard;
end
close all; clear all

% this is set to black in frvmCcc.m and frvmvskap.m
set(0,'DefaultLineMarkerEdgeColor','auto');
testmnumadiabat;

fprintf(['Display the difference between the correlations '...
	 'for ''propane'' and ''propaneperry''.\n']);
twosubstances;
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf('To remove figures, type "close all". To resume, type "return".\n');
  keyboard;
else
  close all;
end

qtestsubstance;

for sname = {'propane', 'propaneperry', 'isobutane', 'butane', 'nitrogen',...
             'R142b'}
  testsubstance(sname{1});
  if exist('INTERACTIVE','var') && INTERACTIVE
    fprintf('To remove figures, type "close all". To resume, type "return".\n');
    keyboard;
  else
    close all;
  end
end

testgaseous;
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf('To remove figures, type "close all". To resume, type "return".\n');
  keyboard;
else
  close all;
end

testcalctaubeta;
if exist('INTERACTIVE','var') && INTERACTIVE
  fprintf('To remove figures, type "close all". To resume, type "return".\n');
  keyboard;
else
  close all;
end
clear all;

test2ph;
fprintf('\n');

testcond;

fprintf('All tests passed.\n');
