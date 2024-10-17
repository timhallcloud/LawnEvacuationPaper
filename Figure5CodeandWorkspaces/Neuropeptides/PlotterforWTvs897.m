figure; hold on;
bar([1 2],[mean(WT) mean(Line897)]/5)
errorbar([1 2],[mean(WT) mean(Line897)]/5,[std(WT)/sqrt(length(WT)) std(Line897)/sqrt(length(Line897))]/5,'.k')
scatter(ones(length(WT),1),WT/4,'.k')
scatter(2*ones(length(Line897),1),Line897/4,'.k')
xticks([1 2])
xticklabels({'WT','pdf-2'})
ylabel('Entry Timescale (hrs)')