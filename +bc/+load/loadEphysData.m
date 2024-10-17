function [spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, ...
    pcFeatures, pcFeatureIdx, channelPositions] = loadEphysData(ephys_path, datasetidx)
% JF, Load ephys data (1-indexed)
% ------
% Inputs
% ------
% ephys_path: character array defining the path to your kilosorted output files 
% datasetidx: 1 x 1 double vector, only use if you have chronically
% recorded stitched datasets. 
% ------
% Outputs
% ------
% spikeTimes_samples: nSpikes × 1 uint64 vector giving each spike time in samples (*not* seconds)
% spikeTemplates: nSpikes × 1 uint32 vector giving the identity of each
%   spike's matched template
% templateWaveforms: nTemplates × nTimePoints × nChannels single matrix of
%   template waveforms for each template and channel
% templateAmplitudes: nSpikes × 1 double vector of the amplitude scaling factor
%   that was applied to the template when extracting that spike
% pcFeatures: nSpikes × nFeaturesPerChannel × nPCFeatures  single
%   matrix giving the PC values for each spike
% pcFeatureIdx: nTemplates × nPCFeatures uint32  matrix specifying which
%   channels contribute to each entry in dim 3 of the pc_features matrix
% channelPositions: nChannels x 2 double matrix, each row gives the x and y 
%   coordinates of each channel
%

% load spike templates (= waveforms)
if exist(fullfile([ephys_path filesep 'spike_clusters.npy']),'file')
   spike_cluster = readNPY([ephys_path filesep 'spike_clusters.npy'])+1; % already manually-curated / Ks4-n KS4, "spike_templates" is called "spike_clusters"
   spikeTemplates = readNPY([ephys_path filesep 'spike_templates.npy'])+1;
   if max(spike_cluster) >= max(spikeTemplates)
       % user did some manual curation cluster split on this recording
       % session, so n_clusters is larger than n_templates
       %
       % to solve this, we need to see how many new clusters were generated
       % first, and then for each new cluster, we need to go back and see
       % which old cluster can be matched. Jingjie Li
       spikeTemplates_uq = unique(spikeTemplates);
       new_cluster_index = max(spikeTemplates):max(spike_cluster);
    
       new_cluster_matched_template_index = 0*new_cluster_index';
       for ii = 1:numel(new_cluster_index)
           item_index = find(spike_cluster==new_cluster_index(ii));
           if isempty(item_index)
               % empty cluster index from manual curation
               new_cluster_matched_template_index(ii) = 1;
           else
               new_cluster_matched_template_index(ii) = spikeTemplates(item_index(1));
           end
       end
       spikeTemplates_uq = [spikeTemplates_uq;new_cluster_matched_template_index];
   end
   spikeTemplates = spike_cluster;
else 
  spike_templates_0idx = readNPY([ephys_path filesep 'spike_templates.npy']); 
  spikeTemplates = spike_templates_0idx + 1;
 end


% load spike times 
if exist(fullfile(ephys_path,'spike_times_corrected.npy'),'file') % When running pyKS stitched you need the 'aligned / corrected' spike times
    spikeTimes_samples = double(readNPY([ephys_path filesep  'spike_times_corrected.npy']));
    spikeTimes_datasets = double(readNPY([ephys_path filesep  'spike_datasets.npy'])) + 1; %  which dataset? (zero-indexed so +1)
else
    spikeTimes_samples = double(readNPY([ephys_path filesep  'spike_times.npy']));
    spikeTimes_datasets = ones(size(spikeTimes_samples));
end

templateAmplitudes = double(readNPY([ephys_path filesep 'amplitudes.npy'])); % ensure double (KS4 saves as single)


% Load and unwhiten templates
templateWaveforms_whitened = readNPY([ephys_path filesep 'templates.npy']);
if exist('spikeTemplates_uq','var')
    templateWaveforms_whitened = templateWaveforms_whitened(spikeTemplates_uq, :, :);
end
winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
templateWaveforms = zeros(size(templateWaveforms_whitened));
for t = 1:size(templateWaveforms,1)
    templateWaveforms(t,:,:) = squeeze(templateWaveforms_whitened(t,:,:))*winv;
end

if exist(fullfile([ephys_path filesep  'pc_features.npy']))
    pcFeatures = readNPY([ephys_path filesep  'pc_features.npy']);
    pcFeatureIdx = readNPY([ephys_path filesep  'pc_feature_ind.npy']) + 1;
    if exist('spikeTemplates_uq','var')
        pcFeatureIdx = pcFeatureIdx(spikeTemplates_uq, :);
    end
else  % not computed in early kilosort3 version - the distance and drift metrics (which are based on the PCs) will not be calculated 
    pcFeatures = NaN;
    pcFeatureIdx = NaN;
end 
channelPositions = readNPY([ephys_path filesep  'channel_positions.npy']) ; 
%goodChannels = readNPY([ephys_path filesep  'channel_map.npy']) + 1;


%% Only use data set of interest - for unit match
if nargin > 2  %- for unit match
   
    spikeTimes_samples = spikeTimes_samples(spikeTimes_datasets == datasetidx);
    spikeTemplates = spikeTemplates(spikeTimes_datasets == datasetidx);
    templateAmplitudes = templateAmplitudes(spikeTimes_datasets == datasetidx);
    pcFeatures=pcFeatures(spikeTimes_datasets == datasetidx,:,:);
end

end