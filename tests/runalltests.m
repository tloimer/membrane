%RUNALLTESTS A simple reminder to run all tests. Stay in the directory tests/.

[path,name] = fileparts(cd);
if name ~= 'tests', error('Run from directory tests/.'); end
clear('path','name');

if ~exist('substance.m'), addpath('../program'); end

%global VERBOSE;  VERBOSE = 1;
entropy
%fprintf('To remove figures, type "close all". To resume, type "return".\n');
fprintf('To resume, type "return".\n');
keyboard
close all; clear all

frvmCcc
fprintf('To resume, type "return".\n');
keyboard

frvmvskap
fprintf('To resume, type "return".\n');
keyboard
close all; clear all

testmstack
fprintf('To resume, type "return".\n');
keyboard
close all; clear all

% this is set to black in frvmCcc.m and frvmvskap.m
set(0,'DefaultLineMarkerEdgeColor','auto');
testmnumadiabat;
