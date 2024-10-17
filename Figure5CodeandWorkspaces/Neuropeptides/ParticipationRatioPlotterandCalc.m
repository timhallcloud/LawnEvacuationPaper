%%
load('WorkspaceforheatmapsNPandNPR.mat')
%Calculate entropy of all NPs
%Convert gene expression to probability
TotalExpression = sum(NParray');
ProbabilityNP = [];
for i = 1:length(Neuronsunique)
    ProbabilityNP = [ProbabilityNP NParray(:,i)./TotalExpression'];
end
test = ProbabilityNP.*log(ProbabilityNP);
test(isnan(test))=0
[ValueEntropy,OrderEntropy] = sort(-sum(test'));


[ValueIPR,OrderIPR] = sort(1./sum(ProbabilityNP'.^2));
IPR=(1./sum(ProbabilityNP'.^2));
Participationratio = 1./IPR;


scaledNParray = [];
for i = 1:length(NParray)
    scaledNParray = [scaledNParray; NParray(i,:)./max(NParray(i,:))];
end

temp = scaledNParray(:,find(Neuronsunique == 'AIY'));
figure; hold on; scatter(temp(temp>=0.01),Participationratio(temp>=0.01),'x')
text(temp(temp>=0.01),Participationratio(temp>=0.01),genes(temp>=0.01),'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('Relative Expression')
ylabel('Participation Ratio')
xlim([0 1])
ylim([0 1])
box on
title('AIY Neuropeptide Candidates')

temp = scaledNParray(:,find(Neuronsunique == 'AVK'));
figure; hold on; scatter(temp(temp>=0.01),Participationratio(temp>=0.01),'x')
text(temp(temp>=0.01),Participationratio(temp>=0.01),genes(temp>=0.01),'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('Relative Expression')
ylabel('Participation Ratio')
xlim([0 1])
ylim([0 1])
box on
title('AVK Neuropeptide Candidates')

temp = scaledNParray(:,find(Neuronsunique == 'SIA'));
figure; hold on; scatter(temp(temp>=0.01),Participationratio(temp>=0.01),'x')
text(temp(temp>=0.01),Participationratio(temp>=0.01),genes(temp>=0.01),'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('Relative Expression')
ylabel('Participation Ratio')
xlim([0 1])
ylim([0 1])
box on
title('SIA Neuropeptide Candidates')


temp = scaledNParray(:,find(Neuronsunique == 'MI'));
figure; hold on; scatter(temp(temp>=0.01),Participationratio(temp>=0.01),'x')
text(temp(temp>=0.01),Participationratio(temp>=0.01),genes(temp>=0.01),'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('Relative Expression')
ylabel('Participation Ratio')
xlim([0 1])
ylim([0 1])
box on
title('MI Neuropeptide Candidates')

temp = max([scaledNParray(:,find(Neuronsunique == 'AIY')) scaledNParray(:,find(Neuronsunique == 'SIA')) scaledNParray(:,find(Neuronsunique == 'AVK'))]');
figure; hold on; scatter(temp(temp>=0.01),Participationratio(temp>=0.01),'x')
text(temp(temp>=0.01),Participationratio(temp>=0.01),genes(temp>=0.01),'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('Relative Expression')
ylabel('Participation Ratio')
xlim([0 1])
ylim([0 1])
box on
title('SIA/AIY/AVK Neuropeptide Candidates')

temp = max([scaledNParray(:,find(Neuronsunique == 'AIY')) scaledNParray(:,find(Neuronsunique == 'SIA')) scaledNParray(:,find(Neuronsunique == 'AVK')) scaledNParray(:,find(Neuronsunique == 'MI'))]');
figure; hold on; scatter(temp(temp>=0.01),Participationratio(temp>=0.01),'x')
text(temp(temp>=0.01),Participationratio(temp>=0.01),genes(temp>=0.01),'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('Relative Expression')
ylabel('Participation Ratio')
xlim([0 1])
ylim([0 1])
box on
title('SIA/AIY/AVK/MI Neuropeptide Candidates')
