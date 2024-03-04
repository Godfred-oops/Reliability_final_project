function plotSDR = SdrIMplot(intensity_Levels,PSDA_existingData,PSDA_retrofittedData)

%plot expected SDRmax vs IM for B1&B3 retrofiited and existing buildings
plotSDR = figure; 
for i  = 1:size(PSDA_existingData,2)
    %plot of the existing B1 buildings
    subplot(2,2,i)
    plot(intensity_Levels, PSDA_existingData(:,i))
    xlabel('Ground Motion Intensity Measure (IM)')
    ylabel('Expected SDR_{max}')
    if i == 1
        title('Expected SDR_{max} vs IM for B1 existing buildings at the cripple level')
    else
        title('Expected SDR_{max} vs IM for B1 existing buildings at the 1st story level')
    end
 
    %plot of the retrofitting B1 buildings
    subplot(2,2,i+2)
    plot(intensity_Levels, PSDA_retrofittedData(:,i))
    xlabel('Ground Motion Intensity Measure (IM)')
    ylabel('Expected SDR_{max}')
    if i+2 == 3
        title('Expected SDR_{max} vs IM for B1 retrofitted buildings at the cripple level')
    else
        title('Expected SDR_{max} vs IM for B1 retrofitted buildings at the 1st story level')
    end

end

end
