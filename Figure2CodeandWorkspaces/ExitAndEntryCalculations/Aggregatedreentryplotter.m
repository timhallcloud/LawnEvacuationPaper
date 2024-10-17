%%
figure;
bar([6 4 5 1 3 2]-0.15,[mean(NoNPR4) mean(NoSAMS5) mean(NoAIY) mean(NoDOP2) mean(NoMPZ1) mean(NoFLP4)],0.2,'b')
hold on
bar([6 4 5 1 3 2]+0.15,[mean(ATRNPR4) mean(ATRSAMS5) mean(ATRAIY) mean(ATRDOP2) mean(ATRMPZ1) mean(ATRFLP4)],0.2,'r')
errorbar([6 4 5 1 3 2]-0.15,[mean(NoNPR4) mean(NoSAMS5) mean(NoAIY) mean(NoDOP2) mean(NoMPZ1) mean(NoFLP4)],[std(NoNPR4)/sqrt(length(NoNPR4)) std(NoSAMS5)/sqrt(length(NoSAMS5)) std(NoAIY)/sqrt(length(NoAIY)) std(NoDOP2)/sqrt(length(NoDOP2)) std(NoMPZ1)/sqrt(length(NoMPZ1)) std(NoFLP4)/sqrt(length(NoFLP4))],'.k')
errorbar([6 4 5 1 3 2]+0.15,[mean(ATRNPR4) mean(ATRSAMS5) mean(ATRAIY) mean(ATRDOP2) mean(ATRMPZ1) mean(ATRFLP4)],[std(ATRNPR4)/sqrt(length(ATRNPR4)) std(ATRSAMS5)/sqrt(length(ATRSAMS5)) std(ATRAIY)/sqrt(length(ATRAIY)) std(ATRDOP2)/sqrt(length(ATRDOP2)) std(ATRMPZ1)/sqrt(length(ATRMPZ1)) std(ATRFLP4)/sqrt(length(ATRFLP4))],'.k')
xticks([1:6])
xticklabels({'dop-2','flp4','mpz1','sams-5','ttx-3','npr-4'})
scatter((6-0.15)*ones(1,length(NoNPR4)),NoNPR4,'.k')
scatter((6+0.15)*ones(1,length(ATRNPR4)),ATRNPR4,'.k')
scatter((4-0.15)*ones(1,length(NoSAMS5)),NoSAMS5,'.k')
scatter((4+0.15)*ones(1,length(ATRSAMS5)),ATRSAMS5,'.k')
scatter((5-0.15)*ones(1,length(NoAIY)),NoAIY,'.k')
scatter((5+0.15)*ones(1,length(ATRAIY)),ATRAIY,'.k')
scatter((1-0.15)*ones(1,length(NoDOP2)),NoDOP2,'.k')
scatter((1+0.15)*ones(1,length(ATRDOP2)),ATRDOP2,'.k')
scatter((3-0.15)*ones(1,length(NoMPZ1)),NoMPZ1,'.k')
scatter((3+0.15)*ones(1,length(ATRMPZ1)),ATRMPZ1,'.k')
scatter((2-0.15)*ones(1,length(NoFLP4)),NoFLP4,'.k')
scatter((2+0.15)*ones(1,length(ATRFLP4)),ATRFLP4,'.k')
ylim([0 120]);
legend({'-ATR','+ATR'})
ylabel('Entry Timescale (min)');
box on
%%
[h,p] = ttest2(NoNPR4,NPR4)

%%

ylabel('Re-entry time (min)')
title('Re-entry timescale of worms during neural inhibition')
legend({'No ATR','ATR'})



ATRdata = [mean(ATRNPR4) mean(ATRSAMS5) mean(ATRAIY) mean(ATRDOP2) mean(ATRMPZ1) mean(ATRFLP4)];
errorATR = [std(ATRNPR4)/sqrt(length(ATRNPR4)) std(ATRSAMS5)/sqrt(length(ATRSAMS5)) std(ATRAIY)/sqrt(length(ATRAIY)) std(ATRDOP2)/sqrt(length(ATRDOP2)) std(ATRMPZ1)/sqrt(length(ATRMPZ1)) std(ATRFLP4)/sqrt(length(ATRFLP4))];
NoATRdata = [mean(NoNPR4) mean(NoSAMS5) mean(NoAIY) mean(NoDOP2) mean(NoMPZ1) mean(NoFLP4)];
errorNoATR = [std(NoNPR4)/sqrt(length(NoNPR4)) std(NoSAMS5)/sqrt(length(NoSAMS5)) std(NoAIY)/sqrt(length(NoAIY)) std(NoDOP2)/sqrt(length(NoDOP2)) std(NoMPZ1)/sqrt(length(NoMPZ1)) std(NoFLP4)/sqrt(length(NoFLP4))];




DeltaATRandNo = ATRdata - NoATRdata;
ErrorATRandNo = sqrt(errorATR.^2 + errorNoATR.^2);


RatioATRandNo = ATRdata./NoATRdata;
ErrorATRandNoRatio = sqrt((errorATR./ATRdata).^2 + (errorNoATR./NoATRdata).^2).*RatioATRandNo;
