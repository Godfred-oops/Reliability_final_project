clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CEE 244 Structural Reliability Final Project
%% 2/28/24
%% Anthony Yen & Ababio Godfred Opoku
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use given code to load input data
LoadInputData
cd(BaseDirectory)

%% Part a
% plot the expected SDRmax vs expected ground motion SaT1, at the cripple
% wall level and 1st story for B1 and B3 buildings (Existing and
% retrofitted)

% reading expected SDRmax data for B1 & B3 existing buildings
PSDA_existingB1 = PSDADataExisting{1,1}.medianSDR;

PSDA_existingB3 = PSDADataExisting{3,1}.medianSDR;

% reading expected SDRmax data for B1 & B3 existing buildings
PSDA_retrofittedB1 = PSDADataRetrofitted{1,1}.medianSDR;

PSDA_retrofittedB3 = PSDADataRetrofitted{3,1}.medianSDR;

%plotting dimensions
height = 1000;
width = 2000; 

%plot for expected SDR vs IM for B1 existing and retrofitted buildings
figure('Position',[100, 100, width, height]); 
for i  = 1:size(PSDA_existingB1,2)
    %plot of the existing B1 buildings
    subplot(2,2,i)
    plot(intensityLevels, PSDA_existingB1(:,i))
    xlabel('Ground Motion Intensity Measure (IM)')
    ylabel('Expected SDR_{max}')
    if i == 1
        title('Expected SDR_{max} vs IM for B1 existing buildings at the cripple level')
    else
        title('Expected SDR_{max} vs IM for B1 existing buildings at the 1st story level')
    end

    %plot of the retrofitting B1 buildings
    subplot(2,2,i+2)
    plot(intensityLevels, PSDA_retrofittedB1(:,i))
    xlabel('Ground Motion Intensity Measure (IM)')
    ylabel('Expected SDR_{max}')
    if i+2 == 3
        title('Expected SDR_{max} vs IM for B1 retrofitted buildings at the cripple level')
    else
        title('Expected SDR_{max} vs IM for B1 retrofitted buildings at the 1st story level')
    end

end
saveas(gcf, 'b1sdr.png')

%plot for expected SDR vs IM for B3 existing and retrofitted buildings
figure('Position',[100, 100, width, height]); 
for i  = 1:size(PSDA_existingB3,2)-1
    %plot of the existing B3 buildings
    subplot(2,2,i)
    plot(intensityLevels, PSDA_existingB3(:,i))
    xlabel('Ground Motion Intensity Measure (IM)')
    ylabel('Expected SDR_{max}')
    if i == 1
        title('Expected SDR_{max} vs IM for B3 existing buildings at the cripple level')
    else
        title('Expected SDR_{max} vs IM for B3 existing buildings at the 1st story level')
    end

    %plot of the retrofitting B3 buildings
    subplot(2,2,i+2)
    plot(intensityLevels, PSDA_retrofittedB3(:,i))
    xlabel('Ground Motion Intensity Measure (IM)')
    ylabel('Expected SDR_{max}')
    if i+2 == 3
        title('Expected SDR_{max} vs IM for B3 retrofitted buildings at the cripple level')
    else
        title('Expected SDR_{max} vs IM for B3 retrofitted buildings at the 1st story level')
    end

end
saveas(gcf, 'b3sdr.png')
%%  Part b
% plot the expected PFA vs expected ground motion SaT1, at the cripple
% wall level and 1st story for B1 and B3 buildings (Existing and
% retrofitted)

% reading expected SDRmax data for B1 & B3 existing buildings
PFA_existingB1 = PSDADataExisting{1,1}.medianPFA;

PFA_existingB3 = PSDADataExisting{3,1}.medianPFA;

% reading expected SDRmax data for B1 & B3 existing buildings
PFA_retrofittedB1 = PSDADataRetrofitted{1,1}.medianPFA;

PFA_retrofittedB3 = PSDADataRetrofitted{3,1}.medianPFA;

%plot for expected PFA vs IM for B1 existing and retrofitted buildings for
%1st story level
figure('Position',[100, 100, width, height]);
subplot(2,1,1)

plot(intensityLevels, PFA_existingB1(:,1))
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected PFA')
title('Expected PFA vs IM for B1 existing buildings at the 1st story level')

subplot(2,1,2)
plot(intensityLevels, PFA_retrofittedB1(:,1))
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected PFA')
title('Expected PFA vs IM for B1 retrofitted buildings at the 1st story level')
saveas(gcf, 'b1PFA.png')

%plot for expected PFA vs IM for B3 existing and retrofitted buildings for
%1st story level
figure('Position',[100, 100, width, height]);
subplot(2,1,1)

plot(intensityLevels, PFA_existingB3(:,2))
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected PFA')
title('Expected PFA vs IM for B3 existing buildings at the 1st story level')

subplot(2,1,2)
plot(intensityLevels, PFA_retrofittedB3(:,2))
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected PFA')
title('Expected PFA vs IM for B3 retrofitted buildings at the 1st story level')


saveas(gcf, 'b3PFA.png')

%% part c
sdr = 1e-6:0.0003:0.3;
%logStdsdr for existing building 
logSdr_existingB1 = PSDADataExisting{1,1}.logSTDSDR;
logSdr_existingB3 = PSDADataExisting{3,1}.logSTDSDR;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Sdr and logSTDsdr for cripple and 1st story
%level for existing building 
stdSdrB1 = logSdr_existingB1(5, 1:2);
stdSdrB3 = logSdr_existingB3(5, 1:2);

muB1 = PSDA_existingB1(5, 1:2);
muB3 = PSDA_existingB3(5, 1:2);

%logStdsdr for retrofitted building 
logSdr_retrofittedB1 = PSDADataRetrofitted{1,1}.logSTDSDR;
logSdr_retrofittedB3 = PSDADataRetrofitted{3,1}.logSTDSDR;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Sdr and logSTDsdr for cripple and 1st story
%level for retrofitted building 
stdSdrRftB1 = logSdr_retrofittedB1(5, 1:2);
stdSdrRftB3 = logSdr_retrofittedB3(5, 1:2);

muRftB1 = PSDA_retrofittedB1(5, 1:2);
muRftB3 = PSDA_retrofittedB3(5, 1:2);

%probability of exceedance P(SDR > sdr) for B1 for cripple level and 1st
%story level
for i = 1:size(stdSdrB1,2)
    mu1 = muB1(:,i);
    stdsdr1 = stdSdrB1(:,i);
    mu2 = muRftB1(:,i);
    stdsdr2 = stdSdrRftB1(:,i);
    for j = 1:size(sdr, 2)
        prob_sdr_existing(j,i) = 1 - logncdf(sdr(j), log(mu1), stdsdr1);
        prob_sdr_retrofitted(j,i) = 1 - logncdf(sdr(j), log(mu2), stdsdr2);
    end
end

figure('Position',[100, 100, width, height]);
for i  = 1:size(prob_sdr_existing,2)
    %plot of the existing B1 buildings
    subplot(2,2,i)
    loglog(sdr', prob_sdr_existing(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i == 1
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 existing buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 existing buildings at the 1st story level')
    end

    %plot of the retrofitting B1 buildings
    subplot(2,2,i+2)
    loglog(sdr', prob_sdr_retrofitted(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i+2 == 3
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 retrofitted buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 retrofitted buildings at the 1st story level')
    end

end
saveas(gcf, 'prosdrB1.png')

%probability of exceedance P(SDR > sdr) for B3 for cripple level and 1st
%story level
for i = 1:size(stdSdrB3,2)
    mu3 = muB1(:,i);
    stdsdr3 = stdSdrB3(:,i);
    mu4 = muRftB3(:,i);
    stdsdr4 = stdSdrRftB3(:,i);
    for j = 1:size(sdr, 2)
        prob_sdr_existingB3(j,i) = 1 - logncdf(sdr(j), log(mu3), stdsdr3);
        prob_sdr_retrofittedB3(j,i) = 1 - logncdf(sdr(j), log(mu4), stdsdr4);
    end
end

figure('Position',[100, 100, width, height]);
for i  = 1:size(prob_sdr_existingB3,2)
    %plot of the existing B3 buildings
    subplot(2,2,i)
    loglog(sdr', prob_sdr_existingB3(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i == 1
        title('Probability of exceedance P(SDR > sdr) vs sdr for B3 existing buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B3 existing buildings at the 1st story level')
    end

    %plot of the retrofitting B3 buildings
    subplot(2,2,i+2)
    loglog(sdr', prob_sdr_retrofittedB3(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i+2 == 3
        title('Probability of exceedance P(SDR > sdr) vs sdr for B3 retrofitted buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B3 retrofitted buildings at the 1st story level')
    end

end
saveas(gcf, 'prosdrB3.png')

%% part d
pfa = 0.01:0.01:6;

%logStdPfa for existing building 
logPfa_existingB1 = PSDADataExisting{1,1}.logSTDPFA;
logPfa_existingB3 = PSDADataExisting{3,1}.logSTDPFA;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Pfa and logSTDPfa for 1st story
%level for existing building 
stdPfaB1 = logPfa_existingB1(5, 1);
stdPfaB3 = logPfa_existingB3(5, 2);

muPfa_B1 = PFA_existingB1(5, 1);
muPfa_B3 = PFA_existingB3(5, 2);

%logStdPfa for retrofitted building 
logPfa_retrofittedB1 = PSDADataRetrofitted{1,1}.logSTDPFA;
logPfa_retrofittedB3 = PSDADataRetrofitted{3,1}.logSTDPFA;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Pfa and logSTDPfa for cripple and 1st story
%level for retrofitted building 
stdPfaRftB1 = logPfa_retrofittedB1(5, 1);
stdPfaRftB3 = logPfa_retrofittedB3(5, 2);

muPfa_RftB1 = PFA_retrofittedB1(5, 1);
muPfa_RftB3 = PFA_retrofittedB3(5, 2);

%probability of exceedance P(PFA > pfa) for B1 & B3 buildings for 1st
%story level
for i = 1:size(pfa,2)
    %probability of exceedance P(PFA > pfa) for B1 buildings
    prob_pfa_existing(i,1) = 1 - logncdf(pfa(i), log(muPfa_B1), stdPfaB1);
    prob_pfa_retrofitted(i,1) = 1 - logncdf(pfa(i), log(muPfa_RftB1), stdPfaRftB1);

    %probability of exceedance P(PFA > pfa) for B3 buildings
    prob_pfa_existingB3(i,1) = 1 - logncdf(pfa(i), log(muPfa_B3), stdPfaB3);
    prob_pfa_retrofittedB3(i,1) = 1 - logncdf(pfa(i), log(muPfa_RftB3), stdPfaRftB3);

end

%plot for probability of exceedance for B1
prob_pfa_B1(:,1) = prob_pfa_existing;
prob_pfa_B1(:,2) = prob_pfa_retrofitted;
figure('Position',[100, 100, width, height]);
list = {'existing', 'retrofitted'};
for i = 1:size(prob_pfa_B1,2)
    subplot(2,1,i)
    loglog(pfa', prob_pfa_B1(:,i))
    xlabel('Peak Floor Acceleration (pfa)')
    ylabel('Probability of exceedance P(PFA > pfa)')  
    title(sprintf('Probability of exceedance P(PFA > pfa) vs pfa for B1 %s buildings at the 1st floor level', list{i}))
end
saveas(gcf, 'propfaB1.png')


%plot for probability of exceedance for B3
prob_pfa_B3(:,1) = prob_pfa_existingB3;
prob_pfa_B3(:,2) = prob_pfa_retrofittedB3;
figure('Position',[100, 100, width, height]);
list = {'existing', 'retrofitted'};
for i = 1:size(prob_pfa_B3,2)
    subplot(2,1,i)
    loglog(pfa', prob_pfa_B3(:,i))
    xlabel('Peak Floor Acceleration (pfa)')
    ylabel('Probability of exceedance P(PFA > pfa)')  
    title(sprintf('Probability of exceedance P(PFA > pfa) vs pfa for B3 %s buildings at the 1st floor level', list{i}))
end
saveas(gcf, 'propfaB3.png')

%% part e
%collapse fragility 
%values from the problem statement
collapse_existing(:,1) = [1.22; 2.38; 0.81; 1.41];
collapse_existing(:,2) = [0.3; 0.6; 0.4; 0.5]; 

collapse_retrofitted(:,1) = [3.15; 4.44; 2.73; 2.67]; 
collapse_retrofitted(:,2) = [0.3; 0.4; 0.5; 0.5];

log_mean_existing = log(collapse_existing(:,1)) - 0.5*(collapse_existing(:,2).^2);
log_mean_retrofitted = log(collapse_retrofitted(:,1)) - 0.5*(collapse_retrofitted(:,2).^2);

for i = 1:length(log_mean_existing)
    for j = 1:size(intensityLevels,1)
        fragility_exisiting(j,i) = logncdf(intensityLevels(j), log_mean_existing(i), collapse_existing(i,2));
        fragility_retrofitted(j,i) = logncdf(intensityLevels(j), log_mean_retrofitted(i), collapse_retrofitted(i,2));

    end
end
%developing collapse fragility functions 
buildinglist = {'B1', 'B2', 'B3', 'B4'};
figure('Position',[100, 100, width, height])
for i = 1:length(log_mean_existing)
    subplot(2,2, i)
    plot(intensityLevels, fragility_exisiting(:,i))
    hold on
    plot(intensityLevels, fragility_retrofitted(:,i))
    xlabel('Sa_{T} (g)')
    ylabel('Probability of collapse (P_{col})')  
    title(sprintf('Collapse Fragility function for %s existing and retrofitted buildings', buildinglist{i}))
    legend('Existing', 'Retrofitted')

end
saveas(gcf, 'collapsefragility.png')

%% part f
%damage fragility
gypsum_median = [0.021, 0.0071, 0.012];
gypsum_disp = 0.5; 
stucco_median = [0.002, 0.005, 0.015];
stucco_disp = 0.45; 
potableWater_median = [2.25, 4.1];
potableWater_disp = 0.4; 

%damage fragility for gpsum wallboard and exterior stucco
for i = 1:length(gypsum_median)
    fragility_gypsum(:,i) = logncdf(sdr, log(gypsum_median(i)) - 0.5*gypsum_disp^2,gypsum_disp);
    fragility_stucco(:,i) = logncdf(sdr, log(stucco_median(i)) - 0.5*stucco_disp^2, stucco_disp);
end

for i = 1:length(potableWater_median)
    fragility_potableWater(:,i) = logncdf(pfa, log(potableWater_median(i)) - 0.5*potableWater_disp^2,potableWater_disp);
end

%plotting damage fragility for gypsum wallboard
figure
for i = 1:size(fragility_gypsum,2)
    plot(sdr, fragility_gypsum(:,i))
    hold on
end
xlabel('SDR')
ylabel('P(DS|SDR)')
title('Damage fragility for Gypsum Wallboard')
legend('DS1', 'DS2', 'DS3')
hold off
saveas(gcf, 'gypsum_fragility.png')
%plotting damage fragility for exterior stucco 
figure
for i = 1:size(fragility_stucco,2)
    plot(sdr, fragility_stucco(:,i))
    hold on
end
xlabel('SDR')
ylabel('P(DS|SDR)')
title('Damage fragility for exterior stucco')
legend('DS1', 'DS2', 'DS3')
hold off
saveas(gcf, 'stucco_fragility.png')
%plotting damage fragility for potableWater piping
figure
for i = 1:size(fragility_potableWater,2)
    plot(pfa, fragility_potableWater(:,i))
    hold on
end
xlabel('PFA')
ylabel('P(DS|PFA)')
title('Damage fragility for potableWater piping')
legend('DS1', 'DS2')
hold off
saveas(gcf, 'potablewater_fragility.png')
%% part L

fault_start = [2;1];
fault_end = [19;9];

%coordinates for the 51 points along the fault
x_cord_fault = linspace(fault_start(1), fault_end(1), 51);
y_cord_fault = linspace(fault_start(2), fault_end(2), 51);

fault_coord(:,1) = x_cord_fault';
fault_coord(:,2) = y_cord_fault';

for i = 1:size(BuildingCoordinates,1)
    buildingcoord = BuildingCoordinates(i,:);
    for j = 1:size(fault_coord,1)
        %distance from the building site to each 51 epicenter locations
        site_dist(i,j) = norm(buildingcoord' - fault_coord(j,:)');
    end
end
%inverse cdf of magnitudes
Mmin = 6;
Mmax = 8; 

m = rand(100,1);

b =1;

for i = 1:size(m,1)
    y = (1- 10^(-b*(m(i) - Mmin)))/(1- 10^(-b*(Mmax - Mmin)));
    cdfM(i) = Mmin - log10(1 - y*(1 - 10^(-b*(Mmax - Mmin))))/b;
end

