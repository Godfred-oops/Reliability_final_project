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
figure 
for i  = 1:size(PSDA_existingB1,2)
    %plot of the existing B1 buildings
    plot(intensityLevels, PSDA_existingB1(:,i))
    hold on


    %plot of the retrofitting B1 buildings

    plot(intensityLevels, PSDA_retrofittedB1(:,i))
    hold on 

end
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected SDR_{max}')
legend("CrippleWall_{existing}", "CrippleWall_{retrofitted}", "1stStory_{existing}", "1stStory_{retrofitted}")
title("Expected SDR vs IM for B1 existing and retrofitted buildings")
saveas(gcf, 'b1sdr.png')

%plot for expected SDR vs IM for B3 existing and retrofitted buildings
figure
for i  = 1:size(PSDA_existingB3,2)-1
    %plot of the existing B3 buildings

    plot(intensityLevels, PSDA_existingB3(:,i));
    hold on

    %plot of the retrofitting B3 buildings

    plot(intensityLevels, PSDA_retrofittedB3(:,i))
    hold on

end
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected SDR_{max}')
legend("CrippleWall_{existing}", "CrippleWall_{retrofitted}", "1stStory_{existing}", "1stStory_{retrofitted}")
title("Expected SDR vs IM for B3 existing and retrofitted buildings")
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
figure

plot(intensityLevels, PFA_existingB1(:,1))
hold on
plot(intensityLevels, PFA_retrofittedB1(:,1))
hold off
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected PFA')
title('PFA vs IM for B1 existing & retrofitted buildings at the 1st story level')
saveas(gcf, 'b1PFA.png')

%plot for expected PFA vs IM for B3 existing and retrofitted buildings for
%1st story level
figure

plot(intensityLevels, PFA_existingB3(:,2))
hold on
plot(intensityLevels, PFA_retrofittedB3(:,2))
hold off
xlabel('Ground Motion Intensity Measure (IM)')
ylabel('Expected PFA')
title('PFA vs IM for B3 existing & retrofitted buildings at the 1st story level')
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

figure
for i  = 1:2
    %plot of the existing B1 buildings

    cell_pro_existingB1 = prob_sdr_existing{5,1};
    loglog(sdr', cell_pro_existingB1(:,i))
    hold on

    %plot of the retrofitting B1 buildings

    cell_pro_retrofittedB1 = prob_sdr_retrofitted{5,1};
    loglog(sdr', cell_pro_retrofittedB1(:,i))
    hold on

end
title('P(SDR > sdr) vs sdr for B1  buildings')
xlabel('SDR')
ylabel('P(SDR > sdr)')
legend("CrippleWall_{existing}", "CrippleWall_{retrofitted}", "1stStory_{existing}", "1stStory_{retrofitted}",'Location','southwest')
saveas(gcf, 'prosdrB1.png')

figure
for i  = 1:2 %change this later
    %plot of the existing B3 buildings

    cell_pro_existingB3 = prob_sdr_existing{5,3};
    loglog(sdr', cell_pro_existingB3(:,i))
    hold on


    %plot of the retrofitting B3 buildings

    cell_pro_retrofittedB3 = prob_sdr_retrofitted{5,3};
    loglog(sdr', cell_pro_retrofittedB3(:,i))
    hold on

end
title('P(SDR > sdr) vs sdr for B3 buildings')
xlabel('SDR')
ylabel('P(SDR > sdr)')
legend("CrippleWall_{existing}", "CrippleWall_{retrofitted}", "1stStory_{existing}", "1stStory_{retrofitted}", 'Location','southwest')
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

figure
for i = 1:size(prob_pfa_B1,2)
    loglog(pfa', prob_pfa_B1(:,i))
    hold on

end
xlabel('Peak Floor Acceleration (pfa)')
ylabel('P(PFA > pfa)')  
title('P(PFA > pfa) vs pfa for B1 buildings at the 1st floor level')
legend('existing', 'retrofitted')
saveas(gcf, 'propfaB1.png')


%plot for probability of exceedance for B3
prob_pfa_B3(:,1) = prob_pfa_existing{5,3};
prob_pfa_B3(:,2) = prob_pfa_retrofitted{5,3};

figure
for i = 1:size(prob_pfa_B3,2)
    loglog(pfa', prob_pfa_B3(:,i))
    hold on

end
xlabel('Peak Floor Acceleration (pfa)')
ylabel('P(PFA > pfa)')  
title('P(PFA > pfa) vs pfa for B3 buildings at the 1st floor level')
legend('existing', 'retrofitted')
saveas(gcf, 'propfaB3.png')


%% part e
collapse fragility 
values from the problem statement
collapse_existing(:,1) = [1.22; 2.38; 0.81; 1.41];
collapse_existing(:,2) = [0.3; 0.6; 0.4; 0.5]; 

collapse_retrofitted(:,1) = [3.15; 4.44; 2.73; 2.67]; 
collapse_retrofitted(:,2) = [0.3; 0.4; 0.5; 0.5];

log_mean_existing = log(collapse_existing(:,1));
log_mean_retrofitted = log(collapse_retrofitted(:,1));

for i = 1:length(log_mean_existing)
    for j = 1:size(intensityLevels,1)
        fragility_exisiting(j,i) = logncdf(intensityLevels(j), log_mean_existing(i), collapse_existing(i,2));
        fragility_retrofitted(j,i) = logncdf(intensityLevels(j), log_mean_retrofitted(i), collapse_retrofitted(i,2));

    end
end
%developing collapse fragility functions 
buildinglist = {'B1', 'B2', 'B3', 'B4'};

for i = 1:length(log_mean_existing)
    figure
    plot(intensityLevels, fragility_exisiting(:,i))
    hold on
    plot(intensityLevels, fragility_retrofitted(:,i))
    hold on 
    xlabel('Sa_{T} (g)')
    ylabel('Probability of collapse (P_{col})')  
    title(sprintf('Collapse Fragility function for %s buildings', buildinglist{i}))
    legend('Existing', 'Retrofitted')
    filename = sprintf('collapsefragility_%s.png', buildinglist{i});
    saveas(gcf, filename)

end


%% part f
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
    fragility_gypsum(:,i) = logncdf(sdr, log(gypsum_median(i)),gypsum_disp);
    fragility_stucco(:,i) = logncdf(sdr, log(stucco_median(i)), stucco_disp);
    fragility_cripple(:,i) = logncdf(sdr, log(cripple_median(i)), cripple_disp);

end

for i = 1:length(potableWater_median)
    fragility_potableWater(:,i) = logncdf(pfa, log(potableWater_median(i)),potableWater_disp);
    fragility_hvac(:,i) = logncdf(pfa, log(hvac_median(i)) ,hvac_disp);
    fragility_sprinker_water(:,i) = logncdf(pfa, log(sprinkler_water_median(i)) ,sprinkler_water_disp);
    fragility_sprinkler_drop(:,i) = logncdf(pfa, log(sprinkler_drop_median(i)) ,sprinkler_drop_disp);

end

fragility_cermic(:,1) = logncdf(sdr, log(cermic_median),cermic_disp);
fragility_wallpaper(:,1) = logncdf(sdr, log(wallpaper_median) ,wallpaper_disp);
fragility_pendant(:,1) = logncdf(pfa, log(pendant_median) ,pendant_disp);

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

%% part g
%storing 9 components of the building
comp_database.CompName = ["gypsum","stucco", "cermic", "wallpaper", "pendant", "piping", "hvac", "sprinkler_piping", "sprinkler_drop"];
comp_database.EDP = ["SDR","SDR","SDR","SDR","PFA","PFA", "PFA","PFA","PFA"];
comp_database.CompFrag = {{[0.0021, 0.5], [0.0071, 0.5], [0.012, 0.5]},{[0.0020, 0.45], [0.005, 0.45], [0.015, 0.45]},...
    {[0.0021, 0.6]}, {[0.0021, 0.6]}, ...
    {[1.5, 0.4]}, {[2.25, 0.4], [4.1, 0.4]}, {[1.5, 0.4], [2.25, 0.4]}, ...
    {[1.5, 0.4], [2.6, 0.4]}, {[1.5, 0.4], [2.25, 0.4]}};

comp_database.RepairCost = {{3937,8658,27824},{1195,1696,5225},{34992},{3240},...
    {79},{56,536},{364,3556},{72,543},{168,168}};
comp_database.Qty = [1.5,12,0.25,0.5,4,0.5,0.25,0.25,0.1];

%storing cripple info
comp_database_cripple = {[0.0053, 0.45], [0.0132, 0.45], [0.0396, 0.45]};
comp_database_repairCripple = {{1195,1696,5225}};
comp_cripple_qty = 2;

%non-collapse loss for existing & retrofitted buildings
%first colum is existing and second column is retrofitted

total_loss = zeros(10, length(intensityLevels));
buildingNonCollapse = cell(4,2);

combinedPSDA=[PSDADataExisting, PSDADataRetrofitted];

for a = 1:size(combinedPSDA,2)
    %extracting the building condition (existing or retrofitted)
    existing = combinedPSDA(:,a);
    for b = 1:numberOfBuildings
            for i = 1:length(intensityLevels)

                for k = 1:size(existing{b,1}.medianSDR,2)
                    if k ~=1
                        %non-collapse for the other story other than cripple wall
                        for j = 1:length(comp_database.CompName)
                            comp_name = comp_database.CompName(j);
                            comp_EDP = comp_database.EDP(j);
                            comp_frag = comp_database.CompFrag{j};
                            comp_consq = comp_database.RepairCost{j};
                            comp_qty = comp_database.Qty(j);

                            if comp_EDP == "SDR"

                                edp_values = sdr; 
                                edp_median = existing{b,1}.medianSDR(:,k);
                                edp_disp = existing{b,1}.logSTDSDR(:,k);
                                edp_median_idx = edp_median(i);
                                edp_disp_idx = edp_disp(i);
                            elseif comp_EDP == "PFA" 

                                edp_values = pfa; 
                                edp_median = existing{b,1}.medianPFA(:,k-1);
                                edp_disp = existing{b,1}.logSTDPFA(:,k-1);
                                edp_median_idx = edp_median(i);
                                edp_disp_idx = edp_disp(i);

                            end

                            % Loop over all damage states for the current component 
                            E_Lc_IMNC1 = 0;
                            comp_DSnum = length(comp_frag);


                            for m = 1:comp_DSnum
                                comp_frag_DS = comp_frag{m}; % the fragility parameters under the current DS for the current component
                                 comp_consq_DS = comp_consq{m}; % the mean loss under the current DS for the current component

                                % Determine the probability of attaining (versus exceedance) the current damage state for the current component (via Eq. 5 in the tutorial)
                                clear Pr_DS_EDP
                                if m == comp_DSnum
                                    Pr_DS_EDP = normcdf(log(edp_values/comp_frag_DS(1))/comp_frag_DS(2));
                                else
                                    comp_frag_DS_next = comp_frag{m+1};
                                    Pr_DS_EDP = normcdf(log(edp_values/comp_frag_DS(1))/comp_frag_DS(2))- normcdf(log(edp_values/comp_frag_DS_next(1))/comp_frag_DS_next(2));
                                end

                                % calculate Pr(EDP=edp)*d(edp)
                                Pr_EDP = [abs(diff(normcdf(log(edp_values/edp_median_idx)/edp_disp_idx))), 1-normcdf(log(edp_values(end)/edp_median_idx)/edp_disp_idx)];

                                % Multiply Pr_DS_EDP by Pr_EDP to obtain Pr(DS=ds|IM,NC) (eq.4)
                                Pr_DS_IM = Pr_DS_EDP * Pr_EDP';

                                %Append the loss under the current DS(Eq.3)
                                E_Lc_IMNC_total = E_Lc_IMNC1 + Pr_DS_IM*comp_qty*comp_consq_DS;

                            end % end DS_idx

                            Loss_comp_IM(j+1,k-1) = E_Lc_IMNC_total;  

                        end


                    else
                           %cripple wall calculations
                           %values for the cripple wall
                            comp_frag = {[0.0053, 0.45], [0.0132, 0.45], [0.0396, 0.45]};
                            comp_consq = {1195,1696,5225};
                            comp_qty = 2;

                            edp_values = sdr; 
                            edp_median = existing{b,1}.medianSDR(:,1);
                            edp_disp = existing{b,1}.logSTDSDR(:,1);
                            edp_median_idx = edp_median(i);
                            edp_disp_idx = edp_disp(i);

                            % Loop over all damage states for the current component 
                            E_Lc_IMNC = 0; % initiate the loss for the current component under the current IM. It will later be accumulated during the loop (via Eq.3 in the tutorial)
                            comp_DSnum = length(comp_frag);


                            for j = 1:comp_DSnum
                            comp_frag_DS = comp_frag{j}; % the fragility parameters under the current DS for the current component
                            comp_consq_DS = comp_consq{j}; % the mean loss under the current DS for the current component

                            % Determine the probability of attaining (versus exceedance) the current damage state for the current component (via Eq. 5 in the tutorial)
                            clear Pr_DS_EDP
                            if j == comp_DSnum
                                Pr_DS_EDP = normcdf(log(edp_values/comp_frag_DS(1))/comp_frag_DS(2));
                            else
                                comp_frag_DS_next = comp_frag{j+1};
                                Pr_DS_EDP = normcdf(log(edp_values/comp_frag_DS(1))/comp_frag_DS(2))- normcdf(log(edp_values/comp_frag_DS_next(1))/comp_frag_DS_next(2));
                            end

                            % calculate Pr(EDP=edp)*d(edp)
                            Pr_EDP = [abs(diff(normcdf(log(edp_values/edp_median_idx)/edp_disp_idx))), 1-normcdf(log(edp_values(end)/edp_median_idx)/edp_disp_idx)];

                            % Multiply Pr_DS_EDP by Pr_EDP to obtain Pr(DS=ds|IM,NC) (eq.4)
                            Pr_DS_IM = Pr_DS_EDP * Pr_EDP';

                            %Append the loss under the current DS(Eq.3)
                            E_Lc_IMNC = E_Lc_IMNC + Pr_DS_IM*comp_qty*comp_consq_DS;
                            end


                            Loss_comp_IM(1,k) = E_Lc_IMNC;

                    end


                end
                if size(Loss_comp_IM,2) == 1
                    total_loss(:,i) = Loss_comp_IM;
                else
                    total_loss(:,i) = Loss_comp_IM(:,1) + Loss_comp_IM(:,2);
                end
                clear Loss_comp_IM
            end
            buildingNonCollapse{b,a} = total_loss; 
    end
end

%probability of non-collapse existing and retrofitted buildings
for i = 1:size(fragility_exisiting,2)
    NonCollapse_fragility_existing(:,i) = 1-fragility_exisiting(:,i);
    NonCollapse_fragility_retrofitted(:,i) = 1 - fragility_retrofitted(:,i);
end

%existing & retrofitted non-collapse loss
for i = 1:size(buildingNonCollapse,1)
    non_collapse_index = buildingNonCollapse(i,:);
    Non_collapse_existing(:,i) = NonCollapse_fragility_existing(:,i) .* sum(cell2mat(non_collapse_index(:,1)),1)';
    Non_collapse_retrofitted(:, i) = NonCollapse_fragility_retrofitted(:,i) .* sum(cell2mat(non_collapse_index(:,2)),1)';

end


%total for both existing and retrofit for B1
exist_B1 = Non_collapse_existing(:,1);
retrofit_B1 = Non_collapse_retrofitted(:,1);

%total for both existing and retrofit for B4
exist_B4 = Non_collapse_existing(:,4);
retrofit_B4 = Non_collapse_retrofitted(:,4);

figure
plot(intensityLevels, exist_B1)
hold on 
plot(intensityLevels, retrofit_B1)
hold off
xlabel("IM")
ylabel("E(L_T|NC)")
title("Non-collapse loss for B1 buildings")
legend("B1_{Existing}", "B1_{Retrofitted}", 'Location','northwest')
saveas(gcf, 'nonCollpaseB1.png')

figure
plot(intensityLevels, exist_B4)
hold on 
plot(intensityLevels, retrofit_B4)
hold off
xlabel("IM")
ylabel("E(L_T|NC)")
title("Non-collapse loss for B4 buildings")
legend("B4_{Existing}", "B4_{Retrofitted}", 'Location','northwest')
saveas(gcf, 'nonCollpaseB4.png')
%% part h) 

collapse_existing(:,1) = [1.22; 2.38; 0.81; 1.41];
collapse_existing(:,2) = [0.3; 0.6; 0.4; 0.5]; 

collapse_retrofitted(:,1) = [3.15; 4.44; 2.73; 2.67]; 
collapse_retrofitted(:,2) = [0.3; 0.4; 0.5; 0.5];

% compute the total expected collapse (direct) losses at each intensity 
% level building types B1 and B4 (existing and retrofitted)
% collapse_existing = [median,std] collapse_retrofitted = [median,std]
collapse_loss_existing = zeros(length(intensityLevels),numberOfBuildings);
collapse_loss_retrofitted = zeros(length(intensityLevels),numberOfBuildings);
% estimated replacement cost per building type 
replacement_cost = zeros(numberOfBuildings,1);
% building types B1 and B2 are 1 story buildings, B3 and B4 are 2 story 
% cost = construction costs (not including temporary housing)
replacement_cost(1,1) = (150 * 1400); % B1
replacement_cost(2,1) = (150 * 1400); % B2
replacement_cost(3,1) = (150 * 2800); % B3
replacement_cost(4,1) = (150 * 2800); % B4

% cost = construction costs + cripple wall retrofit
replacement_cost_retrofitted(1,1) = (150 + 4) * 1400; % B1
replacement_cost_retrofitted(2,1) = (150 + 4) * 1400; % B2
replacement_cost_retrofitted(3,1) = (150 + 4) * 2800; % B3
replacement_cost_retrofitted(4,1) = (150 + 4) * 2800; % B4

% solve for existing buildings
for i = 1:numberOfBuildings % for 4 building types
    CollapseMed = collapse_existing(i,1);
    CollapseDisp = collapse_existing(i,2);
    reconstruction_cost = replacement_cost(i);

    collapse_loss_existing(:,i) = reconstruction_cost*normcdf(log(intensityLevels/CollapseMed)/CollapseDisp);
end 
% solve for retrofitted buildings
for i = 1:numberOfBuildings % for 4 building types
    CollapseMed = collapse_retrofitted(i,1);
    CollapseDisp = collapse_retrofitted(i,2);
    reconstruction_cost = replacement_cost_retrofitted(i);

    collapse_loss_retrofitted(:,i) = reconstruction_cost*normcdf(log(intensityLevels/CollapseMed)/CollapseDisp);
end 

% plot for building types B1 (existing and retrofitted)
figure

plot(intensityLevels,collapse_loss_existing(:,1))
hold on
plot(intensityLevels,collapse_loss_retrofitted(:,1))
hold off
title("Expected Collapse Loss for B1 buildings")
xlabel("Intensity Level")
ylabel("Expected Collapse Loss")
legend("B1_{existing}", "B1_{retrofitted}", 'Location','northwest')
saveas(gcf, 'collapseB1.png')

figure
plot(intensityLevels,collapse_loss_existing(:,4))
hold on 
plot(intensityLevels,collapse_loss_retrofitted(:,4))
hold off
title("Expected Collapse Loss for B4 buildings")
xlabel("IM")
ylabel("Expected Collapse Loss")
legend("B4_{existing}", "B4_{retrofitted}",'Location','northwest')
saveas(gcf, 'collapseB4.png')

%% part i
% use total probability theorem to compute the expected direct loss at each intensity
% measure and report the results in a plot for B1 and B4 (existing and retrofitted)

%sum of the loss for the collapse and non-collapse E(LT|SA, NC)P(NC|SA) +
%E(LT|SA,C)P(C|SA) = E(LT|SA)
total_directLoss_existing = Non_collapse_existing + collapse_loss_existing;
total_directLoss_retrofitted = Non_collapse_retrofitted + collapse_loss_retrofitted;

%direct loss for B1 buildings (existing and retrofitted)
figure
plot(intensityLevels,collapse_loss_existing(:,1))
hold on 
plot(intensityLevels,collapse_loss_retrofitted(:,1))
hold on
plot(intensityLevels, total_directLoss_existing(:,1), '--')
hold on
plot(intensityLevels, total_directLoss_retrofitted(:,1), '--')
title("Direct Loss for existing and retrofitted B1 Buildings")
xlabel("IM")
ylabel("Direct Collapse Loss")
legend("B1_{existing}", "B1_{retrofitted}", "B1_{DirectLossExisting}", "B1_{DirectLossRetrofitted}",'Location','southeast')
saveas(gcf, 'DirectB1.png')

figure
plot(intensityLevels,collapse_loss_existing(:,4))
hold on 
plot(intensityLevels,collapse_loss_retrofitted(:,4))
hold on
plot(intensityLevels, total_directLoss_existing(:,4), '--')
hold on
plot(intensityLevels, total_directLoss_retrofitted(:,4), '--')
title("Direct Loss for existing and retrofitted B4 buildings")
xlabel("IM")
ylabel("Direct Collapse Loss")
legend("B4_{existing}", "B4_{retrofitted}", "B4_{DirectLossExisting}", "B4_{DirectLossRetrofitted}",'Location','southeast')
saveas(gcf, 'DirectB4.png')
%% part j
% Use the replacement times provided to plot the  indirect loss (due to temporary housing costs) 
% versus SaT1 for building types B1 and B4 (existing and retrofitted).

% read txt file to matrix [SaT1 Existing Retrofitted]
SaT1VsRecoveryTime = readmatrix("SaT1VsRecoveryTime.txt");
% replacement times for each building type in days ([B1 B2 B3 B4])
replacement_time = 30*[9 9 14 14];
% temporary housing cost for each building type
temp_housing_cost = [50 50 75 75];
% calculate indirect loss and plot values
indirect_loss_existing = zeros(length(SaT1VsRecoveryTime),numberOfBuildings);
indirect_loss_retrofitted = zeros(length(SaT1VsRecoveryTime),numberOfBuildings);
% existing
for i = 1:numberOfBuildings
    indirect_loss_existing(:,i) = SaT1VsRecoveryTime(:,2) * replacement_time(i) * temp_housing_cost(i);
end 
% retrofitted
for i = 1:numberOfBuildings
    indirect_loss_retrofitted(:,i) = SaT1VsRecoveryTime(:,3) * replacement_time(i) * temp_housing_cost(i);
end 
%plot of the indirect losses for existing and retrofitted B1 

figure
plot(SaT1VsRecoveryTime(:,1),indirect_loss_existing(:,1))
hold on
plot(SaT1VsRecoveryTime(:,1),indirect_loss_retrofitted(:,1))
hold off
title("Indirect Loss for existing and retrofitted B1 buildings")
xlabel("Sa_{T1}")
ylabel("Indirect Loss")
legend("B1_{existing}", "B1_{retrofitted}",'Location','northwest')
saveas(gcf, 'saT1recoveryB1.png')

figure
plot(SaT1VsRecoveryTime(:,1),indirect_loss_existing(:,4))
hold on
plot(SaT1VsRecoveryTime(:,1),indirect_loss_retrofitted(:,4))
hold off
title("Indirect Loss for existing and retrofitted B4 buildings")
xlabel("Sa_{T1}")
ylabel("Indirect Loss")
legend("B4_{existing}", "B4_{retrofitted}",'Location','northwest')
saveas(gcf, 'saT1recoveryB4.png')

%% part k
% Compute the total (indirect + direct) loss versus SaT1 for B1 and B4 (existing and retrofitted)
% Overlay the plots from parts (i), (j) and (k) in the same figure (one for each building type).

for i = 1:length(intensityLevels)
    for j = 1:size(indirect_loss_retrofitted,2)
            inter_indirect_loss_existing(i,j) = interp1(SaT1VsRecoveryTime(:,1),indirect_loss_existing(:,j), intensityLevels(i), 'cubic');
            inter_indirect_loss_retrofitted(i,j) = interp1(SaT1VsRecoveryTime(:,1),indirect_loss_retrofitted(:,j), intensityLevels(i));
    end
end

%total loss (direct + indirect loss)
total_loss_retrofitted = inter_indirect_loss_retrofitted + total_directLoss_retrofitted;
total_loss_existing = inter_indirect_loss_existing + total_directLoss_existing;

%plot for B1
figure
plot(intensityLevels,total_loss_existing(:,1))
hold on
plot(intensityLevels,total_loss_retrofitted(:,1))
hold on
plot(intensityLevels,inter_indirect_loss_existing(:,1), '--')
hold on
plot(intensityLevels,inter_indirect_loss_retrofitted(:,1), '--')
hold on
plot(intensityLevels,total_directLoss_existing(:,1), '*-')
hold on
plot(intensityLevels,total_directLoss_retrofitted(:,1), '*-')


title("Total Loss for existing and retrofitted B1 buildings")
xlabel("IM")
ylabel("Direct Loss")
legend("B1_{TotalLossExisting}", "B1_{TotalLossRetrofitted}","B1_{IndirectLossExisting}", "B1_{IndirectLossRetrofitted}", "B1_{directLossExisting}","B1_{directLossRetrofitted}",'Location','southeast')
saveas(gcf, 'totalLossB1.png')

%plot for B4
figure
plot(intensityLevels,total_loss_existing(:,4))
hold on
plot(intensityLevels,total_loss_retrofitted(:,4))
hold on
plot(intensityLevels,inter_indirect_loss_existing(:,4), '--')
hold on
plot(intensityLevels,inter_indirect_loss_retrofitted(:,4), '--')
hold on
plot(intensityLevels,total_directLoss_existing(:,4), '*-')
hold on
plot(intensityLevels,total_directLoss_retrofitted(:,4), '*-')


title("Total Loss for existing and retrofitted B4 buildings")
xlabel("IM")
ylabel("Direct Loss")
legend("B4_{TotalLossExisting}", "B4_{TotalLossRetrofitted}","B4_{IndirectLossExisting}", "B4_{IndirectLossRetrofitted}", "B4_{directLossExisting}","B4_{directLossRetrofitted}",'Location','southeast')
saveas(gcf, 'totalLossB4.png')

part L

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

%reshaping the shaking intensities matrix
for i = 1:size(shake_IM,2)
    shakeIM_index = shake_IM{i};
    k = 1;
    for j = 1:size(m,2)
        reshape_shakeIM{j,i} = shakeIM_index(k:k+49, :);
        k = k+50;
    end
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

figure
%hazard curve for 1st building/site
loglog(PGA,lambda_SaT1(1,:))
ylabel('Annual Rate of Exceedance')
xlabel('Sa_{T1}')
title('Hazard Curve for building site 1')
saveas(gcf, 'hazard-curve1.png')
%hazard curve for 25th building/site
figure
loglog(PGA,lambda_SaT1(25,:))
ylabel('Annual Rate of Exceedance')
xlabel('Sa_{T1}')
title('Hazard Curve for building site 25')
saveas(gcf, 'hazard_curve25.png')
%hazard curve for 50th building/site
figure
loglog(PGA,lambda_SaT1(50,:))
ylabel('Annual Rate of Exceedance')
xlabel('Sa_{T1}')
title('Hazard Curve for building site 50')
saveas(gcf, 'hazard_curve50.png')

% %% part m
% % regional loss curve 
