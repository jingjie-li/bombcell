%% ~~ Example bombcell pipeline ~~
% Adjust the paths in the 'set paths' section and the parameters in bc_qualityParamValues
% This pipeline will:
%   (1) load your ephys data, 
%   (2) decompress your raw data if it is in .cbin format 
%   (3) run bombcell on your data and save the output and
%   (4) bring up summary plots and a GUI to flip through classified cells.
% The first time, this pipeline will be significantly slower (10-20' more)
% than after because it extracts raw waveforms. Subsequent times these
% pre-extracted waveforms are simply loaded in.
% We recommend running this pipeline on a few datasets and deciding on
% quality metric thresholds depending on the summary plots (histograms 
% of the distributions of quality metrics for each unit) and GUI. 

%% set paths - EDIT THESE 
recording_dir = '/media/jingjie/spike/spk_sorting/sndmap/JLI-R-0042_2024-09-03_09-35-00_g0';
% [recording_info] = extract_xml_recording_config(recording_dir);
% [path_info] = find_oe_recording_dir(recording_dir);
[path_info] = utils.find_npx_recording_dir(recording_dir);
path_info = path_info(1);
ephysKilosortPath = path_info.spk_sorting_path;% path to your kilosort output files 
ephysRawDir = dir(path_info.recording_data); % path to your raw .bin or .dat data
ephysMetaDir = dir(path_info.meta_data_ap); % path to your .meta or .oebin meta file
savePath = ephysKilosortPath; % where you want to save the quality metrics 
decompressDataLocal = fullfile(path_info.recording_data_path,'temp'); % where to save raw decompressed ephys data 
gain_to_uV = NaN; % use this if you are not using spikeGLX or openEphys to record your data. You then must leave the ephysMetaDir 
    % empty(e.g. ephysMetaDir = '')
kilosortVersion = 4;% if using kilosort4, you need to change this value. Otherwise it does not matter. 
%% check MATLAB version 
if exist('isMATLABReleaseOlderThan', 'file') == 0 % function introduced in MATLAB 2020b.
    oldMATLAB = true;
else
    oldMATLAB = isMATLABReleaseOlderThan("R2019a");
end
if oldMATLAB
    error('This MATLAB version is older than 2019a - download a more recent version before continuing')
end

%% load data 
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath);

%% detect whether data is compressed, decompress locally if necessary
rawFile = bc.dcomp.manageDataCompression(ephysRawDir, decompressDataLocal);

%% which quality metric parameters to extract and thresholds 
param = bc.qm.qualityParamValues(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV, kilosortVersion); %for unitmatch, run this:
% param = qualityParamValuesForUnitMatch(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV)

% payparticular attention to param.nChannels
% param.nChannels must correspond to the total number of channels in your raw data, including any sync channels. For Neuropixels probes, this value should typically be either 384 or 385 channels. param.nSyncChannels must correspond to the number of sync channels you recorded. This value is typically 1 or 0.
param.nChannels = 385;
param.nSyncChannels = 1;

%% compute quality metrics 
rerun = 0; %whether to re-run (and save) quality metrics if they are already present
qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

if qMetricsExist == 0 || rerun
    [qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
else
    [param, qMetric] = bc.load.loadSavedMetrics(savePath); 
    unitType = bc.qm.getQualityUnitType(param, qMetric, savePath);
end

%% view units + quality metrics in GUI 
% load data for GUI
loadRawTraces = 0; % default: don't load in raw data (this makes the GUI significantly faster)
bc.load.loadMetricsForGUI;

% GUI guide: 
% left/right arrow: toggle between units 
% g : go to next good unit 
% m : go to next multi-unit 
% n : go to next noise unit
% up/down arrow: toggle between time chunks in the raw data
% u: brings up a input dialog to enter the unit you want to go to

% currently this GUI works best with a screen in portrait mode - we are
% working to get it to handle screens in landscape mode better. 
unitQualityGuiHandle = bc.viz.unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
    param, probeLocation, unitType, loadRawTraces);


%% example: get the quality metrics for one unit
% this is an example to get the quality metric for the unit with the
% original kilosort and phy label of xx (0-indexed), which corresponds to
% the unit with qMetric.clusterID == xx + 1, and to
% qMetric.phy_clusterID == xx . This is *NOT NECESSARILY* the
% (xx + 1)th row of the structure qMetric - some of the  clusters that kilosort
% outputs are empty, because they were dropped in the last stages of the
% algorithm. These empty clusters are not included in the qMetric structure
% there are two ways to do this: 
% 1:
original_id_we_want_to_load = 0;
id_we_want_to_load_1_indexed = original_id_we_want_to_load + 1; 
number_of_spikes_for_this_cluster = qMetric.nSpikes(qMetric.clusterID == id_we_want_to_load_1_indexed);
% or 2:
original_id_we_want_to_load = 0;
number_of_spikes_for_this_cluster = qMetric.nSpikes(qMetric.phy_clusterID == original_id_we_want_to_load);


%% example: get unit labels 
% the output of `unitType = getQualityUnitType(param, qMetric);` gives
% the unitType in a number format. 1 indicates good units, 2 indicates mua units, 3
% indicates non-somatic units and 0 indciates noise units (see below) 
 
goodUnits = unitType == 1;
muaUnits = unitType == 2;
noiseUnits = unitType == 0;
nonSomaticUnits = unitType == 3; 

% example: get all good units number of spikes
all_good_units_number_of_spikes = qMetric.nSpikes(goodUnits);

% (for use with another language: output a .tsv file of labels. You can then simply load this) 
label_table = table(unitType);
writetable(label_table,[savePath filesep 'templates._bc_unit_labels.tsv'],'FileType', 'text','Delimiter','\t');  
      

%% optional: additionally compute ephys properties for each unit and classify cell types 
rerunEP = 0; %whether to re-run (and save) ephys properties if they are already present
region = 'Cortex'; % options include 'Striatum' and 'Cortex'
[ephysProperties, unitClassif] = bc.ep.runAllEphysProperties(ephysKilosortPath, savePath, rerunEP, region);

% example: get good MSN units 
goodMSNs = strcmp(unitClassif, 'MSN') & unitType == 1; 


