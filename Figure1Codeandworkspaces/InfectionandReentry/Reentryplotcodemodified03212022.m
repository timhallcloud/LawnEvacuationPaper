cd('C:\Users\R123a\Desktop\Figure1\Figure1addition02282023\1E')
load('C:\Users\R123a\Desktop\Figure1\Figure1addition02282023\1E\Reentryfiguresmatlabdata.mat')

%%

[H,L] = size(Sick);
Distance = {};
for n = 1:L/3
    tempx = Sick(:,(3*(n-1)+2));
    tempx = tempx(~isnan(tempx));
    tempy = Sick(:,(3*(n-1)+3));
    tempy = tempy(~isnan(tempy));
    tempDistances = sqrt((tempx - 657).^2 + (tempy - 260).^2) - (84.017/2);
    endpoint = min(find(tempDistances<0));
    if isempty(endpoint)
    else
        tempDistances = tempDistances(1:endpoint);
    end
    Distance{n} = tempDistances;
end

%interpolate to find 0 cross
Distancenew = {};
Timenew = {};
index = 1;
for j = 1:length(Distance)
    if j~=6
        tempslope = (Distance{j}(end-1) - Distance{j}(end));
        tempadditionaltime = abs(Distance{j}(end-1)/tempslope);
        tempdistance = [Distance{j}(1:end-1);0];
        timetemp = [1:(length(Distance{j})-1) length(Distance{j})+tempadditionaltime];
    end
    Distancenew{index} = tempdistance;
    Timenew{index} = timetemp;
    index = index+1;
end

figure; hold on;
for j = 1:length(Distancenew)
    plot(Timenew{j}/20,Distancenew{j}*(26/390),'r')
    scatter(Timenew{j}(end)/20,Distancenew{j}(end)*(26/390),'rx')
end

xlim([0 5])
ylim([-0.5 1.2])
box on;
x0=100;
y0=100;
width=450;
height=400;
set(gcf,'position',[x0,y0,width,height])

%%
[H,L] = size(Healthy);
Distance = {};
for n = 1:L/3
    tempx = Healthy(:,(3*(n-1)+2));
    tempx = tempx(~isnan(tempx));
    tempy = Healthy(:,(3*(n-1)+3));
    tempy = tempy(~isnan(tempy));
    tempDistances = sqrt((tempx - 681).^2 + (tempy - 754).^2) - (86/2);
    endpoint = min(find(tempDistances<0));
    if isempty(endpoint)
    else
        tempDistances = tempDistances(1:endpoint);
    end
    Distance{n} = tempDistances;
end

%interpolate to find 0 cross
hold on
Distancenew = {};
Timenew = {};
index = 1;
for j = 1:length(Distance)
    tempslope = (Distance{j}(end-1) - Distance{j}(end));
    tempadditionaltime = abs(Distance{j}(end-1)/tempslope);
    tempdistance = [Distance{j}(1:end-1);0];
    timetemp = [1:(length(Distance{j})-1) length(Distance{j})+tempadditionaltime];
    Distancenew{index} = tempdistance;
    Timenew{index} = timetemp;
    index = index+1;
end

for j = 1:length(Distancenew)
    plot(Timenew{j}/20,Distancenew{j}*(26/390),'k')
    scatter(Timenew{j}(end)/20,Distancenew{j}(end)*(26/390),'kx')
end

xlim([0 5])
ylim([-0.5 1.2])
box on;
x0=100;
y0=100;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])
