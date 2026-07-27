% RUN_INITIALIZE_ND2FILES
%
% Assumptions:
%   - All ND2 files are in one folder.
%   - Alphabetical order is the acquisition order.
%   - The first file is always hyb 0.
%   - Every ND2 series is one FOV and all series are extracted.
%   - Later files are MERFISH rounds plus optional territory rounds.
%   - Planes are interleaved by channel within each Z position:
%       Z1-C1, Z1-C2, ..., Z2-C1, Z2-C2, ...
%
% Output folders:
%   hyb0_channels   : DAPI, optional WFdapi, P-granule, optional lineage
%   fiducial_beads  : created only if a fiducial channel is specified
%   merfish_probe   : created only if MERFISH rounds exist
%   territory_probe : created only if territoryRounds > 0
%
% Files and variables are named by biological role rather than wavelength.
% Requires bfmatlab/Bio-Formats on the MATLAB path.

clear;
close all;
clc;
tic;

%% --------------------------- USER INPUTS -------------------------------

nd2Folder = "C:\Users\james\JamesSciNet\20260321_ANS9_RFPb_mNG_GAM488_R4_Denoised";

hyb0OutputFolder = "hyb0_channels";
fiducialOutputFolder = "fiducial_beads";
merfishOutputFolder = "merfish_probe";
territoryOutputFolder = "territory_probe";

% Territory rounds can occur before or after the MERFISH rounds.
territoryRounds = 0;
territoryPosition = "End";    % "Before" or "End"

% Typical hyb 0 channel order when both DAPI acquisitions are present:
%   1 = 640, P-granules
%   2 = 560, fiducial beads
%   3 = 488, lineage
%   4 = SD405, default DAPI
%   5 = WF405, optional widefield DAPI
hyb0ChannelCount = 4;
pGranuleChannel = 1;
hyb0FiducialChannel = 2;      % Use [] if absent
lineageChannel = 3;           % Use [] if absent
sdDapiChannel = [];            % Use [] if SD405 is absent
wfDapiChannel = 4;            % Use [] if WF405 is absent

% When only WF405 is present, for example:
%   hyb0ChannelCount = 4;
%   sdDapiChannel = [];
%   wfDapiChannel = 4;
%
% In that case WF405 is saved as the default DAPI output.

% Later-round channel order.
roundChannelCount = 2;
merfishChannel = 1;
territoryChannel = 1;
fiducialChannel = 2;          % Use [] if absent

%% ------------------------------ SETUP ----------------------------------

assert(exist("bfopen", "file") == 2, ...
    "bfopen was not found. Add bfmatlab to the MATLAB path.");
assert(isfolder(nd2Folder), "ND2 folder not found: %s", nd2Folder);

territoryPosition = lower(string(territoryPosition));
assert(ismember(territoryPosition, ["before", "end"]), ...
    'territoryPosition must be "Before" or "End".');
assert(territoryRounds >= 0 && mod(territoryRounds, 1) == 0, ...
    "territoryRounds must be a nonnegative integer.");

files = dir(fullfile(nd2Folder, "*.nd2"));
assert(~isempty(files), "No ND2 files found in: %s", nd2Folder);

[~, order] = sort(lower(string({files.name})));
files = files(order);

nLaterRounds = numel(files) - 1;
assert(territoryRounds <= nLaterRounds, ...
    "territoryRounds exceeds the number of files after hyb 0.");

nMerfishRounds = nLaterRounds - territoryRounds;
makeFolderIfNeeded(hyb0OutputFolder);

if territoryRounds > 0
    makeFolderIfNeeded(territoryOutputFolder);
end

if nMerfishRounds > 0
    makeFolderIfNeeded(merfishOutputFolder);
end

if ~isempty(hyb0FiducialChannel) || ~isempty(fiducialChannel)
    makeFolderIfNeeded(fiducialOutputFolder);
end

%% --------------------------- HYB ROUND 0 -------------------------------

hyb0File = fullfile(files(1).folder, files(1).name);
fprintf("Hyb 0: %s\n", hyb0File);

data = bfopen(hyb0File);

for series = 1:size(data, 1)

    planes = data{series, 1};

    % SD405 is preferred/default. If SD405 is absent, WF405 becomes DAPI.
    % When both exist, WF405 is also saved separately as WFdapi.
    if ~isempty(sdDapiChannel)
        saveChannel(planes, hyb0ChannelCount, sdDapiChannel, ...
            "DAPI_0", "DAPI", series, hyb0OutputFolder);

        if ~isempty(wfDapiChannel)
            saveChannel(planes, hyb0ChannelCount, wfDapiChannel, ...
                "WFdapi_0", "WFdapi", series, hyb0OutputFolder);
        end
    elseif ~isempty(wfDapiChannel)
        saveChannel(planes, hyb0ChannelCount, wfDapiChannel, ...
            "DAPI_0", "DAPI", series, hyb0OutputFolder);
    else
        error("At least one DAPI channel must be specified.");
    end

    saveChannel(planes, hyb0ChannelCount, pGranuleChannel, ...
        "PGranule_0", "PGranule", series, hyb0OutputFolder);

    if ~isempty(lineageChannel)
        saveChannel(planes, hyb0ChannelCount, lineageChannel, ...
            "Lineage_0", "Lineage", series, hyb0OutputFolder);
    end

    if ~isempty(hyb0FiducialChannel)
        saveChannel(planes, hyb0ChannelCount, hyb0FiducialChannel, ...
            "Fiducial_0", "Fiducial", series, fiducialOutputFolder);
    end
end

clear data;

%% ----------------------- CLASSIFY LATER ROUNDS -------------------------

isTerritory = false(1, nLaterRounds);

if territoryRounds > 0
    if territoryPosition == "before"
        isTerritory(1:territoryRounds) = true;
    else
        isTerritory(end-territoryRounds+1:end) = true;
    end
end

merfishRound = 0;
territoryRound = 0;

%% ----------------------- EXTRACT LATER ROUNDS --------------------------

for fileIndex = 1:nLaterRounds

    file = files(fileIndex + 1);
    filename = fullfile(file.folder, file.name);

    if isTerritory(fileIndex)
        territoryRound = territoryRound + 1;
        role = "Territory";
        roundNumber = territoryRound;
        signalChannel = territoryChannel;
        outputFolder = territoryOutputFolder;
    else
        merfishRound = merfishRound + 1;
        role = "MERFISH";
        roundNumber = merfishRound;
        signalChannel = merfishChannel;
        outputFolder = merfishOutputFolder;
    end

    fprintf("%s round %d: %s\n", role, roundNumber, filename);
    data = bfopen(filename);

    for series = 1:size(data, 1)

        planes = data{series, 1};

        saveChannel(planes, roundChannelCount, signalChannel, ...
            role + "_" + roundNumber, role, series, outputFolder);

        if ~isempty(fiducialChannel)
            saveChannel(planes, roundChannelCount, fiducialChannel, ...
                "Fiducial_" + role + "_" + roundNumber, ...
                "Fiducial", series, fiducialOutputFolder);
        end
    end

    clear data;
end

fprintf("\nFinished extracting ND2 data.\n");
toc;

%% --------------------------- LOCAL FUNCTIONS ---------------------------

function makeFolderIfNeeded(folder)
    if ~isfolder(folder)
        mkdir(folder);
    end
end

function saveChannel(planes, nChannels, channel, fileStem, role, series, folder)

    assert(channel >= 1 && channel <= nChannels, ...
        "Channel index is outside the specified channel count.");
    assert(mod(size(planes, 1), nChannels) == 0, ...
        "Plane count is not divisible by the specified channel count.");

    stack = cat(3, planes{channel:nChannels:end, 1});

    output = struct();
    output.(char("ImageStack" + role)) = stack;

    outputFile = fullfile(folder, fileStem + "_" + string(series) + ".mat");
    save(char(outputFile), "-struct", "output", "-v7.3");

    fprintf("  Saved %s\n", outputFile);
end
