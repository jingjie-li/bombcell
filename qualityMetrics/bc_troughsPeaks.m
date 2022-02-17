function [nPeaks, nTroughs, somatic, peakLocs, troughLocs, cellLikeDuration, spatialDecayPoints, spatialDecaySlope] = bc_troughsPeaks(templateWaveforms, thisUnit, maxChannel, ephys_sample_rate, plotThis, minWvDuration, maxWvDuration)
% JF, Get the number of troughs and peaks for each waveform, and determine
% whether waveform is likely axonal (biggest peak before biggest trough)
% ------
% Inputs
% ------
% thisWaveform: nTemplates × nTimePoints × nChannels single matrix of
%   template waveforms for each template and channel
% ephys_sample_rate: recording sampling rate (eg 30 000)
% plotThis: boolean, whether to plot waveform and detected peaks or not 
% ------
% Outputs
% ------
% nPeaks: number of detected peaks
% nTroughs: number of detected troughs
% somatic: boolean, is largest detected peak after the largest detected
%   trough (indicative of a somatic spike, cf: Deligkaris, K., Bullmann, T. & Frey, U. 
%   Extracellularly recorded somatic and neuritic signal shapes and classification 
%   algorithms for high-density microelectrode array electrophysiology. Front. Neurosci. 10, 421 (2016).)
% 
thisWaveform = templateWaveforms(thisUnit, :, maxChannel);
minProminence = 0.2 * max(abs(squeeze(thisWaveform))); % minimum threshold to detcet peaks/troughs

[PKS, peakLocs] = findpeaks(squeeze(thisWaveform), 'MinPeakProminence', minProminence); % get peaks

[TRS, troughLocs] = findpeaks(squeeze(thisWaveform)*-1, 'MinPeakProminence', minProminence); % get troughs

if isempty(TRS) % if there is no detected trough, just take minimum value as trough
    TRS = min(squeeze(thisWaveform));
    nTroughs = numel(TRS);
    LOCST_all = find(squeeze(thisWaveform) == TRS);
        if numel(TRS) > 1 % if more than one trough, take the first (usually the correct one) %QQ should change to better:
            % by looking for location where the data is most tightly distributed
            TRS = TRS(1);
        end
    
    troughLocs = find(squeeze(thisWaveform) == TRS);
else
    nTroughs = numel(TRS);
end
if isempty(PKS) % if there is no detected peak, just take maximum value as peak
    PKS = max(squeeze(thisWaveform));
    nPeaks = numel(PKS);
    LOCS_all = find(squeeze(thisWaveform) == PKS);
        if numel(PKS) > 1 % if more than one peak, take the first (usually the correct one) %QQ should change to better:
            % by looking for location where the data is most tightly distributed
            PKS = PKS(1);
        end
    peakLocs = find(squeeze(thisWaveform) == PKS);
else
    nPeaks = numel(PKS);
end


peakLoc = peakLocs(PKS == max(PKS)); %QQ should change to better:
            % by looking for location where the data is most tightly distributed
if numel(peakLoc) > 1
    peakLoc = peakLoc(end);

end
troughLoc = troughLocs(TRS == max(TRS)); %QQ should change to better:
            % by looking for location where the data is most tightly distributed
if numel(troughLoc) > 1
    troughLoc = troughLoc(1);
end

if peakLoc > troughLoc
    somatic = 1;
else
    somatic = 0;
end

if any(abs(thisWaveform(1:10))> 0.05) || abs(troughLoc-peakLoc)*0.0333 < minWvDuration*0.001 || abs(troughLoc-peakLoc)*0.0333 > maxWvDuration*0.001 
    cellLikeDuration = 0;
else
    cellLikeDuration = 1;
end
if maxChannel > 10
    spatialDecayPoints = max(abs(squeeze(templateWaveforms(thisUnit, :, maxChannel:-2:maxChannel-10))));
else
    spatialDecayPoints = max(abs(squeeze(templateWaveforms(thisUnit, :, maxChannel:2:10))));
end
spatialDecaySlope = polyfit(spatialDecayPoints, 1:5,1); 
spatialDecaySlope = spatialDecaySlope(1);
if plotThis
    figure();
    clf;
    % hacky way of plotting - but hey it works
    pbad = plot(1e3*((0:size(thisWaveform, 2) - 1) / ephys_sample_rate), ...
        thisWaveform);
    hold on;
    mm = max(thisWaveform);
    pbb = pbad;
    pbb.XData(1:end) = NaN;
    pbad = plot(1e3*((0:size(thisWaveform, 2) - 1) / ephys_sample_rate), ...
        thisWaveform);

    pbb.XData([peakLocs, troughLocs]) = pbad.XData([peakLocs, troughLocs]);
    %pbb.YData(~[PKS,-TRS])=NaN;

    set(pbb, 'Marker', 'v');
    xlabel('time (ms)')
    ylabel('amplitude (a.u.)')
    legend('detected peaks/troughs')
    makepretty;
end
end

