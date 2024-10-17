% clear all
cd('C:\Users\Bob2.0\Desktop\Figure2CodeandWorkspaces')
load('CorrectlyLabeled01212021aggregatedlawnleavingdataATRnoATR.mat')
ATR = {}; NoATR = {}; Names = {};
ATR{1} = ATR06152019;  NoATR{1} = NoATR06152019; Name{1} = 'Dop-2';
ATR{2} = ATR06282019;  NoATR{2} = NoATR06282019; Name{2} = 'flp-21';
ATR{3} = ATR06302019;  NoATR{3} = NoATR06302019; Name{3} = 'ttx-3';
ATR{4} = ATR11082020; NoATR{4} = NoATR11082020; Name{4} = 'str-2';
ATR{5} = ATR08012019; NoATR{5} = NoATR08012019; Name{5} = 'flp-3';
ATR{6} = ATR08072019; NoATR{6} = NoATR08072019; Name{6} = 'flp-22';
ATR{7} = ATR08182019; NoATR{7} = NoATR08182019; Name{7} = 'flp-19';
ATR{8} = ATR08122019; NoATR{8} = NoATR08122019; Name{8} = 'mpz-prom1';
ATR{9} = ATR09022019; NoATR{9} = NoATR09022019; Name{9} = 'flp-11';
ATR{10} = ATR09092019; NoATR{10} = NoATR09092019; Name{10} = 'flp-12';
ATR{11} = ATR09232019; NoATR{11} = NoATR09232019; Name{11} = 'flp-4';
ATR{12} = ATR09242019; NoATR{12} = NoATR09242019; Name{12} = 'inx-4';
ATR{13} = ATR10052019; NoATR{13} = NoATR10052019; Name{13} = 'flp-7';
ATR{14} = ATR10072019; NoATR{14} = NoATR10072019; Name{14} = 'mod-1';
ATR{15} = ATR10222019; NoATR{15} = NoATR10222019; Name{15} = 'mbr-1';
ATR{16} = ATR10272019; NoATR{16} = NoATR10272019; Name{16} = 'ser2p2';
ATR{17} = ATR11252019; NoATR{17} = NoATR11252019; Name{17} = 'lin-11';
ATR{18} = ATR11292019; NoATR{18} = NoATR11292019; Name{18} = 'sra-11';
ATR{19} = ATR12012019; NoATR{19} = NoATR12012019; Name{19} = 'mgl-1';
ATR{20} = ATR12042019; NoATR{20} = NoATR12042019; Name{20} = 'opt-3';
ATR{21} = ATR12192019; NoATR{21} = NoATR12192019; Name{21} = 'nmr-1';
ATR{22} = ATR01142020; NoATR{22} = NoATR01142020; Name{22} = 'rig-5';
ATR{23} = ATR01272020; NoATR{23} = NoATR01272020; Name{23} = 'npr-4';
ATR{24} = ATR01282020; NoATR{24} = NoATR01282020; Name{24} = 'ser2-p3';
ATR{25} = ATR12172019; NoATR{25} = NoATR12172019; Name{25} = 'odr-2(16)';
ATR{26} = ATR02252020; NoATR{26} = NoATR02252020; Name{26} = 'sams-5';
ATR{27} = ATR02292020; NoATR{27} = NoATR02292020; Name{27} = 'gpa-14';
ATR{28} = ATR03062020; NoATR{28} = NoATR03062020; Name{28} = 'odr-2(18)';
ATR{29} = ATR07172020; NoATR{29} = NoATR07172020; Name{29} = 'mgl-3';
wildtype73_10col = LL03102020;
%%
clearvars -except wildtype73_10col NoATR ATR Name
STP = 1;
ETP = 18;
%%
%Area under curve calculations
%Area between curve initialized
ABC = {};
ABCmean = [];
ABCstdev = [];
for i = 1:length(ATR)
    clear AUCATR; AUCATR = sum(ATR{i}(:,STP:ETP)');
    clear AUCnoATR; AUCnoATR = sum(NoATR{i}(:,STP:ETP)');
    clear ABCtemp; ABCtemp = [];
    for n = 1:5
        for m = 1:5
            ABCtemp = [ABCtemp AUCATR(m) - AUCnoATR(n)];
        end
    end   
    ABC{i} = ABCtemp;
    ABCmean = [ABCmean mean(ABCtemp)];
    ABCstdev = [ABCstdev std(ABCtemp)/sqrt(length(ABCtemp))];
end

%%
%Way 1 of calculating area between curves for control
%split 10 colonies into 2 sets of 5 and compare.
allcombinationsATR = nchoosek([1:10],5);
fivecontrolABC = [];
for i = 1:length(allcombinationsATR)
    allcombinationsNoATR = setdiff([1:10],allcombinationsATR(i,:));
    fivecontrolABC = [fivecontrolABC sum(mean(wildtype73_10col(allcombinationsATR(i,:),[STP:ETP])) - mean(wildtype73_10col(allcombinationsNoATR,[STP:ETP])))];
end

ABCcontrol1 = fivecontrolABC
ABCcontrolmean1 = mean(ABCcontrol1);
ABCcontrolstd1 = std(ABCcontrol1);

%%
%Generate the phenotype graphs
[sortedphenotypes,Sortedpromoterorder] = sort(ABCmean);
figure
hold on
bar([1:29],sortedphenotypes);
errorbar([1:29],sortedphenotypes,ABCstdev(Sortedpromoterorder),'.','Color', 'k')
xticks([1:29])
xticklabels(Name(Sortedpromoterorder));
xtickangle(90)
h=fill([0 0 30 30],[-ABCcontrolstd1 ABCcontrolstd1 ABCcontrolstd1 -ABCcontrolstd1],'k');
h.FaceAlpha=0.3;

%%
%Organize all the neurons
[numbers, TEXT, linenamesfromfile] = xlsread('donelines03172020.xlsx','A2:A30');
[numbers, TEXT,neuronlistpromoters] = xlsread('donelines03172020.xlsx','B2:B30');
AllneuronsNotunique = {};
Stringsplitneuronnames = {};
for i =1:length(neuronlistpromoters)
AllneuronsNotunique = [AllneuronsNotunique strsplit(neuronlistpromoters{i},{';' ' '})];
Stringsplitneuronnames{i} = strsplit(neuronlistpromoters{i},{';' ' '});
end
AllneuronsNotunique
Uniqueneurons = unique(AllneuronsNotunique);

Measurementmatrix = zeros(length(linenamesfromfile),length(Uniqueneurons));
Measurementanotation = {};
for i = 1:length(Stringsplitneuronnames)
    Rowindex = i;
    for j = 1:length(Stringsplitneuronnames{i})
        Colindex = find(strcmp(Uniqueneurons,Stringsplitneuronnames{i}(j)));
        Measurementmatrix(Rowindex,Colindex) = 1;
    end
end

figure;
Measurementmatrixfig = heatmap(Uniqueneurons,linenamesfromfile,Measurementmatrix,'GridVisible','off','ColorbarVisible','off')

%%
fileIDnew = fopen('Uniqueneuronnames.txt','w');
for i = 1:length(Uniqueneurons)-1
    fprintf(fileIDnew,'%s', Uniqueneurons{i})
    fprintf(fileIDnew,',')
end
fprintf(fileIDnew,'%s', Uniqueneurons{end})
fclose(fileIDnew);

%%
fileIDnew = fopen('Uniqueneuronnames.txt','w');
for i = 1:length(Uniqueneurons)
    fprintf(fileIDnew,'%s ', Uniqueneurons{i})
end
fclose(fileIDnew);


%%
%print out data to be used by python script
fileID = fopen('Measurementmatrix10102020.txt','w');
for i = 1:29
    temp = Measurementmatrix(i,:);
    fprintf(fileID,'%d,',round(temp(1:end-1)));
    fprintf(fileID,'%d',round(temp(end)));
    fprintf(fileID,'\n')
end
fclose(fileID);

%%
%do they come from same distribution?
probabilityforeachline = [];
for i = 1:29
    [h,p] = ttest2(ABC{i},ABCcontrol2)
    probabilityforeachline = [probabilityforeachline p];
end
%%
%Generate phenotypes
ZstaticABC = [];
probabilityttest = [];
propagatederror = [];
for i = 1:length(ABC)
    [h,p,ci,stats] = ttest2(ABCcontrol2,ABC{i});
    errorZtemp = sqrt(((std(ABCcontrol2)^2)/(length(ABCcontrol2))) + ((std(ABC{i})^2)/(length(ABC{i}))));
    ZstaticABC = [ZstaticABC ABCmean(i)/errorZtemp];
    probabilityttest = [probabilityttest p];
end
[sortedphenotypes,Sortedpromoterorder] = sort(ZstaticABC);
figure
hold on
bar([1:29],sortedphenotypes);xticks([1:29])
xticklabels(Name(Sortedpromoterorder));
xtickangle(90)



%%
probabilityZ = (2*(1-normcdf(abs(ZstaticABC))));
ZstaticABC1percent = ZstaticABC.*(probabilityZ < 0.01/length(ABC));

[sortedphenotypes,Sortedpromoterorder] = sort(ABCmean);
figure
hold on
bar([1:29],sortedphenotypes);
errorbar([1:29],sortedphenotypes,ABCstdev(Sortedpromoterorder),'.','Color', 'k')
xticks([1:29])
xticklabels(Name(Sortedpromoterorder));
xtickangle(90)
h=fill([0 0 30 30],[-ABCcontrolstd1 ABCcontrolstd1 ABCcontrolstd1 -ABCcontrolstd1],'k');
h.FaceAlpha=0.3;

%%
figure; 
y1 = mean(ATR{1}(:,1:18),1)/10;
y1error = std(ATR{1}(:,1:18),1)/(10*sqrt(5));
y2 = mean(NoATR{1}(:,1:18),1)/10;
y2error = std(NoATR{1}(:,1:18),1)/(10*sqrt(5));
x = [STP:1:ETP];
patch([x fliplr(x)], [y1 fliplr(y2)], 'k','FaceAlpha',.3)
hold on;
errorbar(x,y1,y1error,'r');
errorbar(x,y2,y2error,'b');
ylim([0 1.2])
%legend({'ATR','No ATR'})
plot([1 3],[1.1 1.1],'g','LineWidth',3)
ylabel('Fractional Lawn Occupancy')
xlabel('Time (hours)')
