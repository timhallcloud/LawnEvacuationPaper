load('C:\Users\R123a\Documents\01212021newdatanalysis\CompiledAnalysis10182021\DatacompilationAUCandZthreshold\10182021analysisAUCdata.mat')
clearvars -except Sortedpromoterorder
load('C:\Users\R123a\Documents\01212021newdatanalysis\CompiledAnalysis10182021\DatacompilationAUCandZthreshold\updated01212021aggregatedlawnleavingdataATRnoATR.mat')

figure;
for i = 1:29
    N = Sortedpromoterorder(i);
    subplot(6,5,i);
    if i == 1
        title('dop-2')
    else
        title(Name{N})
    end
    hold on;
    Y = mean(ATR{N})/10;Y = Y(1:18);
    Ye = std(ATR{N})/(10*sqrt(5));Ye = Ye(1:18);
    errorbar([1:18],Y,Ye,'r-')
    Y = mean(NoATR{N})/10;Y = Y(1:18);
    Ye = std(NoATR{N})/(10*sqrt(5));Ye = Ye(1:18);
    errorbar([1:18],Y,Ye,'b-')
    ylim([0 1.1])
    xlim([0 18])
    xticks([0 6 12 18])
end

