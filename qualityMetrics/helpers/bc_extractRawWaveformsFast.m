function [rawWaveformsFull, rawWaveformsPeakChan] = bc_extractRawWaveformsFast(param, spikeTimes_samples, ...
    spikeTemplates, reExtract, savePath, verbose)
% JF, Get raw waveforms for all templates
% ------
% Inputs
% ------
% param with:
% rawFile: raw .bin or .dat file location
% nChannels: number of recorded channels (including sync), (eg 385)
% nSpikesToExtract: number of spikes to extract per template
% rawFile: string containing the location of the raw .dat or .bin file
%
% spikeTimes_samples: nSpikes × 1 uint64 vector giving each spike time in samples (*not* seconds)
% spikeTemplates: nSpikes × 1 uint32 vector giving the identity of each
%   spike's matched template
% verbose: boolean, display progress bar or not
% savePath: where to save output data
% ------
% Outputs
% ------
% rawWaveformsFull: nUnits × nTimePoints × nChannels single matrix of
%   mean raw waveforms for each unit and channel
% rawWaveformsPeakChan: nUnits x 1 vector of each unit's channel with the maximum
%   amplitude

%% Initialize stuff
% Get spike times and indices
nChannels = param.nChannels; % (385)
nSpikesToExtract = param.nRawSpikesToExtract;
spikeWidth = 82;
halfWid = spikeWidth / 2;
dataTypeNBytes = numel(typecast(cast(0, 'uint16'), 'uint8'));
clustInds = unique(spikeTemplates);
nClust = numel(clustInds);

fprintf('Extracting raw waveforms from %s ... \n', param.rawFile)
% Get binary file name
fid = fopen(param.rawFile, 'r');

%% Interate over spike clusters and find all the data associated with them
rawWaveforms = struct;
rawWaveformsFull = nan(nClust, 384, 82);
rawWaveformsPeakChan = nan(nClust, 1);

for iCluster = 1:nClust
    rawWaveforms(iCluster).clInd = clustInds(iCluster);
    rawWaveforms(iCluster).spkInd = spikeTimes_samples(spikeTemplates == clustInds(iCluster));
    if numel(rawWaveforms(iCluster).spkInd) >= nSpikesToExtract
        spksubi = round(linspace(1, numel(rawWaveforms(iCluster).spkInd), nSpikesToExtract))';
        rawWaveforms(iCluster).spkIndsub = rawWaveforms(iCluster).spkInd(spksubi);
    else
        rawWaveforms(iCluster).spkIndsub = rawWaveforms(iCluster).spkInd;
    end
    nSpkLocal = numel(rawWaveforms(iCluster).spkIndsub);

    rawWaveforms(iCluster).spkMap = nan(384, spikeWidth, nSpkLocal);
    for iSpike = 1:nSpkLocal
        spki = rawWaveforms(iCluster).spkIndsub(iSpike);
        bytei = ((spki - halfWid) * nChannels) * dataTypeNBytes;
        fseek(fid, bytei, 'bof');
        data0 = fread(fid, nChannels*spikeWidth, 'int16=>int16'); % read individual waveform from binary file
        frewind(fid);
        data = reshape(data0, nChannels, []);
        %         if whitenBool
        %             [data, mu, invMat, whMat]=whiten(double(data));
        %         end
        %         if size(data, 2) == spikeWidth
        %             rawWaveforms(iCluster).spkMap(:, :, iSpike) = data;
        %         end
        if size(data, 2) == spikeWidth && nChannels == 385
            rawWaveforms(iCluster).spkMap(:, :, iSpike) = data(1:nChannels-1, :, :); %remove sync channel
        elseif size(data, 2) == spikeWidth
            rawWaveforms(iCluster).spkMap(:, :, iSpike) = data(1:nChannels, :, :);
        end

    end
    rawWaveforms(iCluster).spkMapMean = nanmean(rawWaveforms(iCluster).spkMap, 3);
    rawWaveformsFull(iCluster, :, :) = rawWaveforms(iCluster).spkMapMean - mean(rawWaveforms(iCluster).spkMapMean(:, 1:10), 2);

    spkMapMean_sm = smoothdata(rawWaveforms(iCluster).spkMapMean, 1, 'gaussian', 5);

    [~, rawWaveformsPeakChan(iCluster)] = max(max(spkMapMean_sm, [], 2)-min(spkMapMean_sm, [], 2));
    %      figure();
    %     plot(spkMapMean_sm(rawWaveformsPeakChan(iCluster),:))

    if (mod(iCluster, 100) == 0 || iCluster == nClust) && verbose
        fprintf(['\n   Finished ', num2str(iCluster), ' of ', num2str(nClust), ' units.']);
        figure; imagesc(spkMapMean_sm)
        title(['Unit ID: ', num2str(iCluster)]);
        colorbar;
    end

end

fclose(fid);

rawWaveformFolder = dir(fullfile(savePath, 'templates._bc_rawWaveforms.npy'));
if isempty(rawWaveformFolder) || reExtract
    %save(fullfile(spikeFile.folder, 'rawWaveforms.mat'), 'rawWaveforms', '-v7.3');
    writeNPY(rawWaveformsFull, fullfile(savePath, 'templates._bc_rawWaveforms.npy'))
    writeNPY(rawWaveformsPeakChan, fullfile(savePath, 'templates._bc_rawWaveformPeakChannels.npy'))

    if param.saveMultipleRaw
        writeNPY(spikeMap, fullfile(savePath, 'templates._bc_multi_rawWaveforms.npy'))
    end
end

end