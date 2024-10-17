figure;
x = TimeExit;
y = AverageExitRate;
yerror = ErrorExitRate;
yyaxis right
errorbar(x,y,yerror,'r');

hold on;

x = [1:1:20];
y = 1- mean(LawnEvacuationHourly');
yerror = std(LawnEvacuationHourly')/sqrt(5);
yyaxis left
errorbar(x,y(1:20),yerror(1:20),'k');
ylim([0 1.1])

xlim([0 21])

ax = gca;
ax.YAxis(2).Color = 'r';
ax.YAxis(1).Color = 'k';

%%
figure;
x = [5:1:20];
y = yTmean(3:18)/10;
yerror = yTerror(3:18)/10;
yyaxis right
errorbar(x,y,yerror,'b');
ylim([0 45])
hold on;

x = [1:1:20];
y = 1- mean(LawnEvacuationHourly');
yerror = std(LawnEvacuationHourly')/sqrt(5);
yyaxis left
errorbar(x,y(1:20),yerror(1:20),'k');
ylim([0 1.1])

xlim([0 21])

ax = gca;
ax.YAxis(2).Color = 'b';
ax.YAxis(1).Color = 'k';

