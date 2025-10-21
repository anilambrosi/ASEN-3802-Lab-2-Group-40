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
% data table format: [time (s)] [temp at each thermocouple (C)]

% properties(1) = density; (2) = cp; (3) = k; (4) = alpha
data.aluminum.properties(1) = 2810;
data.aluminum.properties(2) = 960;
data.aluminum.properties(3) = 130;
data.aluminum.properties(4) = 130/(2810*960);

data.brass.properties(1) = 8500;
data.brass.properties(2) = 380;
data.brass.properties(3) = 115;

data.steel.properties(1) = 8000;
data.steel.properties(2) = 500;
data.steel.properties(3) = 16.2;

%% part 2 task 1
T_0 = 17.065;
H_an = 91.0858;
x = 0.1238;
L = .1338;
b = @(H,L,n) ((-1)^n)*(8*H*L)/(((2*n-1)*pi)^2);
lambda = @(L,n) (2*n-1)*pi/(2*L);
b_n = zeros(10,1);
lambda_n = zeros(10,1);
fourierTerms1 = zeros(10,1);
fourierTerms1000 = zeros(10,1);
i = 0:10;
u1 = zeros(11,1);
u1(1) = T_0 + H_an*x;
u1000 = zeros(11,1);
u1000(1) = T_0 + H_an*x;


for j=1:10
    b_n(j) = b(H_an,L,j);
    lambda_n(j)=lambda(L,j);
    fourierTerms1(j) = b_n(j)*sin(lambda_n(j)*x)*exp(-data.aluminum.properties(4)*(lambda_n(j)^2));
    fourierTerms1000(j)= b_n(j)*sin(lambda_n(j)*x)*exp(-data.aluminum.properties(4)*(lambda_n(j)^2)*1000);
end

for j = 2:11
    u1(j) = u1(1) + sum(fourierTerms1(1:(j-1)));
    u1000(j) = u1000(1) + sum(fourierTerms1000(1:(j-1)));
end

figure()
plot(i,u1)
hold on
plot(i,u1000)
xlabel('Terms in fourier series')
ylabel(['Temperature at final thermocouple (' char(176) 'C)'])
legend('t=1','t=1000')
ylim([10 35])

fouriernumber(1) = data.aluminum.properties(4)/(L*L);
fouriernumber(2) = data.aluminum.properties(4)*1000/(L*L);