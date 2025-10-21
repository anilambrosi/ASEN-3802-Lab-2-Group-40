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


%{
[Hex,T0,Han] = SteadyStateSlope(data,volts,amps);
% Function call for Part 2 task 2
[u] = slopederivation(data.aluminum.v25(:,1),data.aluminum.properties,T0(1),Han(1,1));  
%Generic function call to get all of the data for each thermocouple
%Change slot 1 for data, slot 2 for material properties, 3 for t0 (change i
%value for each material), 4th slot for H values 

%Plotting to compare data
figure()
plot(data.aluminum.v25(1:339,1),u(:,1))
hold on
plot(data.aluminum.v25(1:339,1),data.aluminum.v25(1:339,2))
%}

figure()
T_0 = 17.065;
H = 54.931;
u = slopederivation(data.aluminum.v25(:,1),data.aluminum.properties,T_0,H);
for i = 1:size(u,2)
hold on
plot(data.aluminum.v25(1:339,1),data.aluminum.v25(1:339,i+1),'green');
hold on
plot(data.aluminum.v25(1:339,1),u(:,i),'blue');
hold on
text(data.aluminum.v25(end-5,1),u(end,i)+0.2,['Thermocouple ' num2str(i)],"HorizontalAlignment","right")
end
ylim([16 25])
xlabel('Time (s)')
ylabel(['Temperature at each thermocouple location(' char(176) 'C)'])
legend('Experimental','Analytical','','','','','','','','','','','','','','','Location','southeast')
print('p2t3_aluminum25','-dpng')

figure()
T_0 = 17.275;
H = 78.272;
u = slopederivation(data.aluminum.v30(:,1),data.aluminum.properties,T_0,H);
for i = 1:size(u,2)
hold on
plot(data.aluminum.v30(1:(end-5),1),data.aluminum.v30(1:(end-5),i+1),'green');
hold on
plot(data.aluminum.v30(1:(end-5),1),u(:,i),'blue');
hold on
text(data.aluminum.v30(end-5,1),u(end,i)+0.2,['Thermocouple ' num2str(i)],"HorizontalAlignment","right")
end
ylim([16 28])
xlim([0 3650])
xlabel('Time (s)')
ylabel(['Temperature at each thermocouple location(' char(176) 'C)'])
legend('Experimental','Analytical','','','','','','','','','','','','','','','Location','southeast')
print('p2t3_aluminum30','-dpng')

figure()
T_0 = 16.602;
H = 104.708;
u = slopederivation(data.brass.v25(:,1),data.brass.properties,T_0,H);
for i = 1:size(u,2)
hold on
plot(data.brass.v25(1:(end-5),1),data.brass.v25(1:(end-5),i+1),'green');
hold on
plot(data.brass.v25(1:(end-5),1),u(:,i),'blue');
hold on
text(data.brass.v25(end-5,1),u(end,i)+0.2,['Thermocouple ' num2str(i)],"HorizontalAlignment","right")
end
xlabel('Time (s)')
ylabel(['Temperature at each thermocouple location(' char(176) 'C)'])
legend('Experimental','Analytical','','','','','','','','','','','','','','','Location','southeast')
print('p2t3_brass25','-dpng')


figure()
T_0 = 16.78;
H = 150.169;
u = slopederivation(data.brass.v30(:,1),data.brass.properties,T_0,H);
for i = 1:size(u,2)
hold on
plot(data.brass.v30(1:(end-5),1),data.brass.v30(1:(end-5),i+1),'green');
hold on
plot(data.brass.v30(1:(end-5),1),u(:,i),'blue');
hold on
text(data.brass.v30(end-5,1),u(end,i)+0.2,['Thermocouple ' num2str(i)],"HorizontalAlignment","right")
end
xlim([0,5571])
xlabel('Time (s)')
ylabel(['Temperature at each thermocouple location(' char(176) 'C)'])
legend('Experimental','Analytical','','','','','','','','','','','','','','','Location','southeast')
print('p2t3_brass30','-dpng')

figure()
T_0 = 15.107;
H = 287.308;
u = slopederivation(data.steel.v22(:,1),data.steel.properties,T_0,H);
for i = 1:size(u,2)
hold on
plot(data.steel.v22(1:(end-5),1),data.steel.v22(1:(end-5),i+1),'green');
hold on
plot(data.steel.v22(1:(end-5),1),u(:,i),'blue');
hold on
text(data.steel.v22(end-5,1),u(end,i)+0.2,['Thermocouple ' num2str(i)],"HorizontalAlignment","right")
end
ylim([15 52])
xlim([0 11540])
xlabel('Time (s)')
ylabel(['Temperature at each thermocouple location(' char(176) 'C)'])
legend('Experimental','Analytical','','','','','','','','','','','','','','','Location','southeast')
print('p2t3_steel22','-dpng')



%% Functions 
% Ignore this
function [slopeAnalytical,T0,Han]  = SteadyStateSlope(data,V,I)
% Define Thermocouple x
xPos(1) = 0;
for i=1:8
        xPos(i+1) = .0127*(i-1)+.034925;
end

% Indexing definitions
namebrac = ["aluminum","brass","steel"];
voltbrac = ["v25","v30","v22"];
k=0;

% Nested for loops to calculate data and plot for each
for i=1:3
    for j=1:2
        % Corrects for Steel V variance
        if i==3
            j=3;
        end
% Experimental Data pull and analysis
desiredData(2:9) = data.(namebrac(i)).(voltbrac(j))(length(data.(namebrac(i)).(voltbrac(j)))-10,[2:9]);
[P,S] = polyfit(xPos(2:9),desiredData(2:9),1);
slopeAnalytical(i,j) = P(1); %Save Data for Table
T0(i,j) = P(2); %Save Data for Table
desiredData(1) = T0(i,j);
[polyY,err] = polyval(P,xPos,S);


% Analytical Model Equations
k=k+1; % Indexing for V & I
Qdot = (V(k) * I(k))/1000;                     
A = pi * (.0254/2)^2;                 
Han(i,j) = Qdot / (data.(namebrac(i)).properties(3) * A);             
T_an = T0(i,j) + Han(i,j) * xPos;

% Plots
L=xPos*100/(xPos(9)+.0254); % Prep to make X-axis factor of L
figure()
plot(L,T_an,'g',LineWidth=1) % Analy
hold on
plot(L,polyY,'b',LineWidth=1)
scatter(L,desiredData,'c',LineWidth=1) % Data pts
plot(L,polyY+2*err,'r--',L,polyY-2*err,'r--') % error
xlim([0 100])
xtickformat('%.0f%%') %make X-axis Factor of L
xlabel('Distance from X_{0} (% Length of rod)')
ylabel("Degrees (Celcius)")
legend("Analytical Slope","Experimental Slope","Experimental Data","95% Error", Location="northwest")
title("Steady State Temperature Distribution Slopes for " +namebrac(i) + " at " + voltbrac(j))
print("Steady Dist Slope " + namebrac(i)+ "_" + voltbrac(j),'-r300','-dpng')
% Break Loop after 1 iteration of steel
        if i==3
            break
        end
    end
end
end

% This function has the error
function [u] = slopederivation(time,props,T0,H)
% Define Constants
l = .01+7*.0127+.034925;
alpha = props(3)/(props(1)*props(2));
for i=1:8
    x(i) = .0127*(i-1)+.034925;
    for j=1:length(time)-5
        summa=0;
        for n=1:15
lambda = (2*n-1)*pi/(2*l);
b = (8*H*l/((2*n-1)^2*pi^2))*(-1)^n;
summa = summa + b*sin(lambda*x(i))*exp(-alpha*time(j)*lambda^2);
        end
        u(j,i) = T0+H*x(i) +summa;
    end
end
end

