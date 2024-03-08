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

PSDA_existingB2 = PSDADataExisting{2,1}.medianSDR;

PSDA_existingB3 = PSDADataExisting{3,1}.medianSDR;

PSDA_existingB4 = PSDADataExisting{4,1}.medianSDR;

% reading expected SDRmax data for B1 & B3 existing buildings
PSDA_retrofittedB1 = PSDADataRetrofitted{1,1}.medianSDR;

PSDA_retrofittedB2 = PSDADataRetrofitted{2,1}.medianSDR;

PSDA_retrofittedB3 = PSDADataRetrofitted{3,1}.medianSDR;

PSDA_retrofittedB4 = PSDADataRetrofitted{4,1}.medianSDR;

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

PFA_existingB2 = PSDADataExisting{2,1}.medianPFA;

PFA_existingB3 = PSDADataExisting{3,1}.medianPFA;

PFA_existingB4 = PSDADataExisting{4,1}.medianPFA;

% reading expected SDRmax data for B1 & B3 existing buildings
PFA_retrofittedB1 = PSDADataRetrofitted{1,1}.medianPFA;

PFA_retrofittedB2 = PSDADataRetrofitted{2,1}.medianPFA;

PFA_retrofittedB3 = PSDADataRetrofitted{3,1}.medianPFA;

PFA_retrofittedB4 = PSDADataRetrofitted{4,1}.medianPFA;

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
logSdr_existingB2 = PSDADataExisting{2,1}.logSTDSDR;
logSdr_existingB3 = PSDADataExisting{3,1}.logSTDSDR;
logSdr_existingB4 = PSDADataExisting{4,1}.logSTDSDR;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Sdr and logSTDsdr for cripple and 1st story
%level for existing building 
cell_stdsdr = cell(4,1);

cell_stdsdr{1} = logSdr_existingB1(1:16, 1:2);
cell_stdsdr{2} = logSdr_existingB2(1:16, 1:2);
cell_stdsdr{3} = logSdr_existingB3(1:16, 1:2);
cell_stdsdr{4} = logSdr_existingB4(1:16, 1:2);

cell_mu =  cell(4,1);
cell_mu{1} = PSDA_existingB1(1:16, 1:2);
cell_mu{2} = PSDA_existingB2(1:16, 1:2);
cell_mu{3} = PSDA_existingB3(1:16, 1:2);
cell_mu{4} = PSDA_existingB4(1:16, 1:2);

%logStdsdr for retrofitted building 
logSdr_retrofittedB1 = PSDADataRetrofitted{1,1}.logSTDSDR;
logSdr_retrofittedB2 = PSDADataRetrofitted{2,1}.logSTDSDR;
logSdr_retrofittedB3 = PSDADataRetrofitted{3,1}.logSTDSDR;
logSdr_retrofittedB4 = PSDADataRetrofitted{4,1}.logSTDSDR;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Sdr and logSTDsdr for cripple and 1st story
%level for retrofitted building 
cell_stdsdrRft = cell(4,1);

cell_stdsdrRft{1} = logSdr_retrofittedB1(1:16, 1:2);
cell_stdsdrRft{2} = logSdr_retrofittedB2(1:16, 1:2);
cell_stdsdrRft{3} = logSdr_retrofittedB3(1:16, 1:2);
cell_stdsdrRft{4} = logSdr_retrofittedB4(1:16, 1:2);


cell_muRft =  cell(4,1);

cell_muRft{1} = PSDA_retrofittedB1(1:16, 1:2);
cell_muRft{2} = PSDA_retrofittedB2(1:16, 1:2);
cell_muRft{3} = PSDA_retrofittedB3(1:16, 1:2);
cell_muRft{4} = PSDA_retrofittedB4(1:16, 1:2);

muRftB1 = PSDA_retrofittedB1(5, 1:2);
muRftB3 = PSDA_retrofittedB3(5, 1:2);

%probability of exceedance P(SDR > sdr) for B1 for cripple level and 1st
%story level

for i = 1:size(intensityLevels,1)

    %existing for cripple
    for j = 1:numberOfBuildings

        prob_sdr_existing{i,j}(:,1) = 1 - logncdf(sdr', log(cell_mu{j}(i,1)), cell_stdsdr{1}(i,1));

        %existing for 1st story
        prob_sdr_existing{i,j}(:,2) = 1 - logncdf(sdr', log(cell_mu{j}(i,2)), cell_stdsdr{1}(i,2));

        %retrofitted for cripple
        prob_sdr_retrofitted{i,j}(:,1) = 1 - logncdf(sdr', log(cell_muRft{j}(i,1)), cell_stdsdrRft{1}(i,1));

        %retrofitted for 1st story
        prob_sdr_retrofitted{i,j}(:,2) = 1 - logncdf(sdr', log(cell_muRft{j}(i,2)), cell_stdsdrRft{1}(i,2));
    end

end

figure('Position',[100, 100, width, height]);
for i  = 1:2
    %plot of the existing B1 buildings
    subplot(2,2,i)
    cell_pro_existingB1 = prob_sdr_existing{5,1};
    loglog(sdr', cell_pro_existingB1(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i == 1
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 existing buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 existing buildings at the 1st story level')
    end

    %plot of the retrofitting B1 buildings
    subplot(2,2,i+2)
    cell_pro_retrofittedB1 = prob_sdr_retrofitted{5,1};
    loglog(sdr', cell_pro_retrofittedB1(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i+2 == 3
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 retrofitted buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B1 retrofitted buildings at the 1st story level')
    end

end
saveas(gcf, 'prosdrB1.png')

figure('Position',[100, 100, width, height]);
for i  = 1:2 %change this later
    %plot of the existing B3 buildings
    subplot(2,2,i)
    cell_pro_existingB3 = prob_sdr_existing{5,3};
    loglog(sdr', cell_pro_existingB3(:,i))
    xlabel('Maximum Story Drfit Ratio (SDR)')
    ylabel('Probability of exceedance P(SDR > sdr)')
    if i == 1
        title('Probability of exceedance P(SDR > sdr) vs sdr for B3 existing buildings at the cripple level')
    else
        title('Probability of exceedance P(SDR > sdr) vs sdr for B3 existing buildings at the 1st story level')
    end

    %plot of the retrofitting B3 buildings
    subplot(2,2,i+2)
    cell_pro_retrofittedB3 = prob_sdr_retrofitted{5,3};
    loglog(sdr', cell_pro_retrofittedB3(:,i))
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
logPfa_existingB2 = PSDADataExisting{2,1}.logSTDPFA;
logPfa_existingB3 = PSDADataExisting{3,1}.logSTDPFA;
logPfa_existingB4 = PSDADataExisting{4,1}.logSTDPFA;

%IM (0.65g) is on the fifth row
%Extracting the 5th median Pfa and logSTDPfa for 1st story
%level for existing building 
stdPfaB1 = logPfa_existingB1(5, 1);
stdPfaB3 = logPfa_existingB3(5, 2);

cell_stdPfa = cell(4,1);
cell_stdPfa{1} = logPfa_existingB1;
cell_stdPfa{2} = logPfa_existingB2;
cell_stdPfa{3} = logPfa_existingB3;
cell_stdPfa{4} = logPfa_existingB4;

cell_muPfa = cell(4,1);

cell_muPfa{1} = PFA_existingB1;
cell_muPfa{2} = PFA_existingB2;
cell_muPfa{3} = PFA_existingB3;
cell_muPfa{4} = PFA_existingB4;

muPfa_B1 = PFA_existingB1(5, 1);
muPfa_B3 = PFA_existingB3(5, 2);

%logStdPfa for retrofitted building 
logPfa_retrofittedB1 = PSDADataRetrofitted{1,1}.logSTDPFA;
logPfa_retrofittedB2 = PSDADataRetrofitted{2,1}.logSTDPFA;
logPfa_retrofittedB3 = PSDADataRetrofitted{3,1}.logSTDPFA;
logPfa_retrofittedB4 = PSDADataRetrofitted{4,1}.logSTDPFA;

cell_stdPfaRft = cell(4,1);
cell_stdPfaRft{1} = logPfa_retrofittedB1;
cell_stdPfaRft{2} = logPfa_retrofittedB2;
cell_stdPfaRft{3} = logPfa_retrofittedB3;
cell_stdPfaRft{4} = logPfa_retrofittedB4;
%IM (0.65g) is on the fifth row
%Extracting the 5th median Pfa and logSTDPfa for cripple and 1st story
%level for retrofitted building 
stdPfaRftB1 = logPfa_retrofittedB1(5, 1);
stdPfaRftB3 = logPfa_retrofittedB3(5, 2);

cell_muPfa_Rft = cell(4,1);

cell_muPfa_Rft{1} = PFA_retrofittedB1;
cell_muPfa_Rft{2} = PFA_retrofittedB2;
cell_muPfa_Rft{3} = PFA_retrofittedB3;
cell_muPfa_Rft{4} = PFA_retrofittedB4;
muPfa_RftB1 = PFA_retrofittedB1(5, 1);
muPfa_RftB3 = PFA_retrofittedB3(5, 2);

%probability of exceedance P(PFA > pfa) for B1 & B3 buildings for 1st

for i = 1:size(intensityLevels,1)

    %existing for cripple
    for j = 1:numberOfBuildings

        %existing for 1st story
        prob_pfa_existing{i,j}(:,1) = 1 - logncdf(pfa', log(cell_muPfa{j}(i,1)), cell_stdPfa{1}(i,1));


        %retrofitted for 1st story
        prob_pfa_retrofitted{i,j}(:,1) = 1 - logncdf(pfa', log(cell_muPfa_Rft{j}(i,1)), cell_stdPfaRft{1}(i,1));
    end

end

%plot for probability of exceedance for B1
prob_pfa_B1(:,1) = prob_pfa_existing{5,1};
prob_pfa_B1(:,2) = prob_pfa_retrofitted{5,1};
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
prob_pfa_B3(:,1) = prob_pfa_existing{5,3};
prob_pfa_B3(:,2) = prob_pfa_retrofitted{5,3};
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

% part f
%damage fragility
gypsum_median = [0.021, 0.0071, 0.012];
gypsum_disp = 0.5; 
stucco_median = [0.002, 0.005, 0.015];
stucco_disp = 0.45; 
potableWater_median = [2.25, 4.1];
potableWater_disp = 0.4; 
cermic_median = 0.0021;
cermic_disp = 0.6; 
wallpaper_median = 0.0021; 
wallpaper_disp = 0.6;
pendant_median = 1.5;
pendant_disp = 0.4;
hvac_median = [1.5, 2.25];
hvac_disp = 0.4;
sprinkler_water_median = [1.5, 2.6];
sprinkler_water_disp = 0.4; 
sprinkler_drop_median = [1.5, 2.25];
sprinkler_drop_disp  = 0.4; 
cripple_median = [0.0053, 0.0132, 0.0396];
cripple_disp = 0.45; 


%damage fragility for gpsum wallboard and exterior stucco
for i = 1:length(gypsum_median)
    fragility_gypsum(:,i) = logncdf(sdr, log(gypsum_median(i)) - 0.5*gypsum_disp^2,gypsum_disp);
    fragility_stucco(:,i) = logncdf(sdr, log(stucco_median(i)) - 0.5*stucco_disp^2, stucco_disp);
    fragility_cripple(:,i) = logncdf(sdr, log(cripple_median(i)) - 0.5*cripple_disp^2, cripple_disp);

end

for i = 1:length(potableWater_median)
    fragility_potableWater(:,i) = logncdf(pfa, log(potableWater_median(i)) - 0.5*potableWater_disp^2,potableWater_disp);
    fragility_hvac(:,i) = logncdf(pfa, log(hvac_median(i)) - 0.5*hvac_disp^2,hvac_disp);
    fragility_sprinker_water(:,i) = logncdf(pfa, log(sprinkler_water_median(i)) - 0.5*sprinkler_water_disp^2,sprinkler_water_disp);
    fragility_sprinkler_drop(:,i) = logncdf(pfa, log(sprinkler_drop_median(i)) - 0.5*sprinkler_drop_disp^2,sprinkler_drop_disp);

end

fragility_cermic(:,1) = logncdf(sdr, log(cermic_median)- 0.5*cermic_disp^2,cermic_disp);
fragility_wallpaper(:,1) = logncdf(sdr, log(wallpaper_median) - 0.5*wallpaper_disp^2,wallpaper_disp);
fragility_pendant(:,1) = logncdf(pfa, log(pendant_median) - 0.5*pendant_disp^2,pendant_disp);

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

y = rand(100,1);

b =1;

for i = 1:size(y,1)
    m(i) = Mmin - log10(1 - y(i)*(1 - 10^(-b*(Mmax - Mmin))))/b;
end

%computing the log-mean for the intensity measure

%converting the 50x51 matrix to a column 
u = 1;
for i = 1:size(site_dist, 1)
    convert_site_dist(u:u+size(site_dist,1)) = site_dist(i,:);
    u = u+size(site_dist,2);
end
PGA = 0.01:0.01:6;

e1 = 0.4383; 
e2 = 0.106;
c1 = -0.5543;
c2 = 0.0195; 
c3 = -0.0075; 
gamma = 0.01; 

std_SA = 0.6;
mean_InSA = cell(size(m,2), 1);
%mean of the PGA (Sa)
for i = 1:size(m,2)
    row_meanInSA = zeros(1, size(convert_site_dist,2)); 
    for j = 1:size(convert_site_dist,2)
            row_meanInSA(1, j) = e1 + e2 *(m(i) - 6.75) + (c1 + c2 *(m(i) - 4.5))*log(convert_site_dist(j)) + c3 *(convert_site_dist(j)-1)+0.5646;

    end
    mean_InSA{i} = reshape(row_meanInSA, [size(site_dist,1),size(site_dist,2)]);
end

%generating shaking intensities per event 
n = 50; %number of samples for monte carlos simulation
monte_sample=rand(n,1);
cellmat_meanInSA = cell2mat(mean_InSA);

for i = 1:size(monte_sample,1)
    for j = 1:size(cellmat_meanInSA,1)
        row_shake_IM(j,:) = logninv(monte_sample(i), cellmat_meanInSA(j,:), std_SA);
    end
    shake_IM{i} = row_shake_IM;
end

% reshape the shaking intensities matrix (50 rows (sites) x 255000 intensities)
site_SaT1 = zeros(50,255000);
shakeIM_expanded = cell2mat(shake_IM);
index = 1;
for i = 1:50:4951
    site_SaT1(1:50,index:index+2549) = shakeIM_expanded(i:i+49,:);
    index = index + 2549;
end

% lambda(SaT1) computed as the fraction of realizations times lambda(M>m_min)
lambda_SaT1 = zeros(50,600);
PGA = 0.01:0.01:6;
for i = 1:50
    for j = 1:length(PGA)
        counter = 0;
        for k = 1:size(site_SaT1,2)
            if site_SaT1(i,k) > PGA(j)
                counter = counter + 1;
            end
        end
        lambda_SaT1(i,j) = gamma* (counter/255000);
    end
end

figure(1)
loglog(PGA,lambda_SaT1(1,:))