Col5exitrateTotal = []
Col5exitevents = Col5(find(~isnan(Col5(:,2))),1)
Col5exitrate = Col5e(find(~isnan(Col5(:,2))),1)

for i = 1:14
    [C,ia,ib] = intersect(Col5exitevents,[60*i:(60*i + 59)]);
    Col5exitrateTotal = [Col5exitrateTotal sum(Col5exitrate(ia))];
end
    

figure; errorbar(Time, op50,op50e,'-b')
hold on
errorbar(Time, PA14,PA14e,'-r')
xlabel('Time (hours)')
ylabel('Fractional Occupancy')
legend({'OP50','PA14'})


