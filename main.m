clear;
clc;
close all;

a=dir('*mA');

for i=1:length(a)
% how to get voltage and amperage from file names?
% - options include strsplit, regex, etc.
% ultimately, we need to use the format of each file name
% 'material'_'volts'V_'amps'mA
b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3
% {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
v = strsplit(b{2},'V'); % volts are always in the second portion
ampval= strsplit(b{3},'mA'); % amps are always in the third portion
volts(i) = str2num(v{1}); % convert string to number (vector)
amps(i) = str2num(ampval{1});
if i <= 2
    data.aluminum.(strcat('v', v{1})) = readmatrix(a(i).name);
elseif i == 3 || i == 4
    data.brass.(strcat('v', v{1})) = readmatrix(a(i).name);
else
    data.steel.(strcat('v', v{1})) = readmatrix(a(i).name);
end
end

% properties(1) = density; properties(2) = cp; properties(3) = k

data.aluminum.properties(1) = 2810;
data.aluminum.properties(2) = 960;
data.aluminum.properties(3) = 130;

data.brass.properties(1) = 8500;
data.brass.properties(2) = 380;
data.brass.properties(3) = 115;

data.steel.properties(1) = 8000;
data.steel.properties(2) = 500;
data.steel.properties(3) = 16.2;
