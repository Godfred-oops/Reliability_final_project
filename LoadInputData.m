% This code will be used to load all of the input data needed for the
% project

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Extract Building Locations and Types                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load BuildingTypesAndCoordinates
BuildingCoordinates = BuildingTypesAndCoordinates(:,1:2);
BuildingNumbers = BuildingTypesAndCoordinates(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Extract Intensity Levels                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intensityLevels = readmatrix('intensityLevels.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Extract PSDA Data                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell array used to store PSDA data. The size of the array will
% be a 4 x 1. Each entry will contain a struct with the PSDA data (existing
% and retrofitted) for a given building.
PSDADataExisting = {};
PSDADataRetrofitted = {};

% Define path to base directory
BaseDirectory = strcat('C:\Users\ababi\Downloads\ClassProjectStarterCodeAndData');

numberOfBuildings = 4;

% Loop over the number of buildings
for i = 1:numberOfBuildings
    
    % Specify the number of stories
    if i == 1 || i == 2
        numberOfStories = 1;
    else
        numberOfStories = 2;
    end
         
    % Define building ID for current building
    BuildingID = sprintf('Building%d',i);

    % Define directory where PSDA data for current building is stored.
    RPSDADataDirectory = strcat(BaseDirectory,'\PSDAData\',...
        BuildingID);

    % Go to existing building directory
    cd(RPSDADataDirectory)
    cd Existing
    
      
    % Extract response data for existing building
    Existing.medianSDR = load('medianSDR.txt');
    Existing.medianPFA = load('medianPFA.txt');
    Existing.logSTDSDR = load('logSTDSDR.txt');
    Existing.logSTDPFA = load('logSTDPFA.txt');
    
    % Append to PSDA data directory 
    PSDADataExisting{i,1} = Existing;
    
    % Go to existing building directory
    cd ..
    cd Retrofitted    
      
    % Extract response data for retrofitted building
    Retrofitted.medianSDR = load('medianSDR.txt');
    Retrofitted.medianPFA = load('medianPFA.txt');
    Retrofitted.logSTDSDR = load('logSTDSDR.txt');
    Retrofitted.logSTDPFA = load('logSTDPFA.txt');

    % Append to PSDA data directory 
    PSDADataRetrofitted{i,1} = Retrofitted;



end