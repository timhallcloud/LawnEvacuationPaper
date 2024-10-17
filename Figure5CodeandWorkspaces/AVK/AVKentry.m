MeanATR = [];
MeanNoATR = [];
errorATR = [];
errorNoATR = [];
for i = 1
    MeanATR = [MeanATR mean(ATR)];
    errorATR = [errorATR std(ATR)/sqrt(length(ATR))];
    MeanNoATR = [MeanNoATR mean(NoATR)];
    errorNoATR = [errorNoATR std(NoATR)/sqrt(length(NoATR))];
end
    

figure;
hold on;
bar([1]-0.25,MeanNoATR/20,0.2,'b')
bar([1]+0.25,MeanATR/20,0.2,'r')

errorbar([1]-0.25,MeanNoATR/20,errorNoATR/20,'.k')
scatter(([1]-0.25)*ones(length(NoATR),1),NoATR/20,'.k');
errorbar([1]+0.25,MeanATR/20,errorATR/20,'.k')
scatter(([1]+0.25)*ones(length(ATR),1),ATR/20,'.k');
xticks([1])
xticklabels({'flp-1'})
ylabel('Re-entry timescale (min)')
legend({'-ATR','+ATR'})


ylim([0 60])