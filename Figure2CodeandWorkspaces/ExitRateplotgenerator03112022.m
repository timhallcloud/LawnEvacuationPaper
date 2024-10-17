Inputs = {};
Inputs{5} = EEX280;
Inputs{6} = EEX352;
Inputs{4} = EEX486;
Inputs{3} = EEX495;
Inputs{2} = EEX496;
Inputs{1} = EEX853;

meanleavingrateATR = [];
errorleavingrateATR = [];
meanleavingrateNoATR = [];
errorleavingrateNoATR = [];

for i = 1:length(Inputs)
    meanleavingrateNoATR = [meanleavingrateNoATR mean(Inputs{i}(:,1))];
    errorleavingrateNoATR = [errorleavingrateNoATR std(Inputs{i}(:,1))/sqrt(5)];
    meanleavingrateATR = [meanleavingrateATR mean(Inputs{i}(:,2))];
    errorleavingrateATR = [errorleavingrateATR std(Inputs{i}(:,2))/sqrt(5)];
end

figure;
bar([1 2 3 4 5 6]-0.15,meanleavingrateNoATR,0.25);
hold on;
bar([1 2 3 4 5 6]+0.15,meanleavingrateATR,0.25);
errorbar([1 2 3 4 5 6]+0.15,meanleavingrateATR,errorleavingrateATR,'.k');
errorbar([1 2 3 4 5 6]-0.15,meanleavingrateNoATR,errorleavingrateNoATR,'.k');

for i = 1:length(Inputs)
    scatter(ones(1,5)*(i-0.15),Inputs{i}(:,1),'.k');
    scatter(ones(1,5)*(i+0.15),Inputs{i}(:,2),'.k');
end
xticks([1:1:6])
xticklabels({'dop-2','flp-4','mpz-1','sams-5','ttx-3','npr-4'})
ylabel('Exit events per worm')
legend({'-ATR','+ATR'})
box on;

print('-painters','-dpdf','ExitEventPhenotypeBargraph03112022v2')



DifferenceLeaving = meanleavingrateATR - meanleavingrateNoATR;
ErrorDIfference = (errorleavingrateNoATR.^2 + errorleavingrateATR.^2).^0.5;