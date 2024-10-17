load('ttx3NPR4SIABfiles.mat')

Lines = {};
Lines{1} = R280;


MeanATR = [];
MeanNoATR = [];
errorATR = [];
errorNoATR = [];
for i = 1:length(Lines)
    Composite = Lines{i};
    NoATR = Composite(:,1);
    NoATR = NoATR(~isnan(NoATR));
    ATR = Composite(:,2);
    ATR = ATR(~isnan(ATR));
    MeanATR = [MeanATR mean(ATR)];
    errorATR = [errorATR std(ATR)/sqrt(length(ATR))];
    MeanNoATR = [MeanNoATR mean(NoATR)];
    errorNoATR = [errorNoATR std(NoATR)/sqrt(length(NoATR))];
end

 

figure;
hold on;
bar([1]-0.15,MeanNoATR/20,0.2,'b')
bar([1]+0.15,MeanATR/20,0.2,'r')
errorbar([1]-0.15,MeanNoATR/20,errorNoATR/20,'.k')
scatter(([1]-0.15)*ones(length(NoATR),1),NoATR/20,'.k');
errorbar([1]+0.15,MeanATR/20,errorATR/20,'.k')
scatter(([1]+0.15)*ones(length(ATR),1),ATR/20,'.k');
xticks([1])
ylabel('Re-entry timescale (min)')
legend({'-ATR','+ATR'})
xticklabels({'Ttx-3 Arch Inhibition'})
xlim([0.5 1.5])
print(gcf, '-dpdf', '-painters', 'Ttx3inhibition');

ylim([0 60])