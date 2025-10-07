clear;
clc;
close all;
% Take a row then 

%function []  = SteadyStateSlope(data,V,I)
xPos(1) = 0;
% Analytical Data
% equations 
Qdot = V * I;                     
A = pi * (D/2)^2;                 
Han = Qdot / (k * A);             
T_an = T0 + Han * x;

% Experimental Data
for i=1:8
        xPos(i+1) = .0127*(i-1)+.01375;
end

namebrac = ["aluminum","brass","steel"];
voltbrac = ["v25","v30"];
V =[25,30,25,30,22];
I = [240 290	237	285	203];
t0=16.5;
    %t0 = data.(namebrac(i)).properties(2);
    desiredData(1) = t0;


for i=1:2


    for j=1:2

desiredData(2:9) = data.(namebrac(i)).(voltbrac(j))(length(data.(namebrac(i)).(voltbrac(j)))-10,[2:9]);
P = polyfit(xPos,desiredData,1);
polyY = polyval(P,xPos);
figure()
plot(xPos,T_an)
hold on
plot(xPos,polyY)

    end
end




