ABCerrorcompilation = [];
ABCfullcompilation = [];

for i = 1:length(ABC)
    ABCerrorcompilation = [ABCerrorcompilation std(ABC{i})/sqrt(25)];
    ABCfullcompilation = [ABCfullcompilation mean(ABC{i})];
end

ABC1percent = ABCfullcompilation.*(probabilityZ < 0.01/length(ABC));
errorABC1percent = ABCerrorcompilation.*(probabilityZ < 0.01/length(ABC));
fileID1 = fopen('1pABC03162021.txt','w');
fprintf(fileID1,'%6.6f\r\n',ABC1percent);
fclose(fileID1);
fileID2 = fopen('1pABC03162021error.txt','w');
fprintf(fileID2,'%6.6f\r\n',errorABC1percent);
fclose(fileID2);
temp1 = probabilityZ < (0.01/length(ABC));
temp2 = probabilityZ < (0.05/length(ABC)) - temp1;
temp3 = ones(1,29) - (temp2+temp1);
[Zvaluessorted,sortedZ] = sort(ZstaticABC);
