%%%
figure;
bar([1 2 3 4 5 6]-0.15,[mean(NoNPR4) mean(NoAIY) mean(NoSAMS5) mean(NoDOP2) mean(NoFLP4) mean(NoMPZ1)],0.25,'b')
hold on
bar([1 2 3 4 5 6]+0.15,[mean(ATRNPR4) mean(ATRAIY) mean(ATRSAMS5) mean(ATRDOP2) mean(ATRFLP4) mean(ATRMPZ1)],0.25,'r')
errorbar([1 2 3 4 5 6]-0.15,[mean(NoNPR4) mean(NoAIY) mean(NoSAMS5) mean(NoDOP2) mean(NoFLP4) mean(NoMPZ1)],[std(NoNPR4)/sqrt(length(NoNPR4)) std(NoAIY)/sqrt(length(NoAIY)) std(NoSAMS5)/sqrt(length(NoSAMS5)) std(NoDOP2)/sqrt(length(NoDOP2)) std(NoFLP4)/sqrt(length(NoFLP4)) std(NoMPZ1)/sqrt(length(NoMPZ1))],'.k')
errorbar([1 2 3 4 5 6]+0.15,[mean(ATRNPR4) mean(ATRAIY) mean(ATRSAMS5) mean(ATRDOP2) mean(ATRFLP4) mean(ATRMPZ1)],[std(ATRNPR4)/sqrt(length(ATRNPR4)) std(ATRAIY)/sqrt(length(ATRAIY)) std(ATRSAMS5)/sqrt(length(ATRSAMS5)) std(ATRDOP2)/sqrt(length(ATRDOP2)) std(ATRFLP4)/sqrt(length(ATRFLP4)) std(ATRMPZ1)/sqrt(length(ATRMPZ1))],'.k')
xticks([1:6])
xticklabels({'npr-4','ttx-3','sams-5','dop-2','flp4','mpz1'})
scatter((1-0.15)*ones(1,length(NoNPR4)),NoNPR4,'.k')
scatter((1+0.15)*ones(1,length(ATRNPR4)),ATRNPR4,'.k')
scatter((2-0.15)*ones(1,length(NoAIY)),NoAIY,'.k')
scatter((2+0.15)*ones(1,length(ATRAIY)),ATRAIY,'.k')
scatter((3-0.15)*ones(1,length(NoSAMS5)),NoSAMS5,'.k')
scatter((3+0.15)*ones(1,length(ATRSAMS5)),ATRSAMS5,'.k')
scatter((4-0.15)*ones(1,length(NoDOP2)),NoDOP2,'.k')
scatter((4+0.15)*ones(1,length(ATRDOP2)),ATRDOP2,'.k')
scatter((5-0.15)*ones(1,length(NoFLP4)),NoFLP4,'.k')
scatter((5+0.15)*ones(1,length(ATRFLP4)),ATRFLP4,'.k')
scatter((6-0.15)*ones(1,length(NoMPZ1)),NoMPZ1,'.k')
scatter((6+0.15)*ones(1,length(ATRMPZ1)),ATRMPZ1,'.k')
ylim([0 120]);
legend({'-ATR','+ATR'})
ylabel('Entry Timescale (min)');
box on
%%
[h,p] = ttest2(NoNPR4,NPR4)