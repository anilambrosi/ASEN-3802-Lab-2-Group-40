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
