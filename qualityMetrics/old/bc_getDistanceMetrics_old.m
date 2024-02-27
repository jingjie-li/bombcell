function [isoD, Lratio, silhouetteScore, d2_mahal, Xplot, Yplot] = bc_getDistanceMetrics_old(pc_features, ...
    pc_feature_ind, thisUnit, numberSpikes, spikesIdx, allSpikesIdx, nChansToUse, plotThis)
% JF, Get distance metrics
% ------
% Inputs
% ------
% pc_features: nSpikes × nFeaturesPerChannel × nPCFeatures  single
%   matrix giving the PC values for each spike.
% pc_feature_ind: nTemplates × nPCFeatures uint32  matrix specifying which
%   channels contribute to each entry in dim 3 of the pc_features matrix
% thisUnit: unit number
% numberSpikes: number of spikes for that unit
% spikesIdx: boolean vector indicating which spikes belong to thisUnit
% allSpikesIdx: nSpikes × 1 uint32 vector giving the identity of each
%   spike's matched template, for all templates
% nChansToUse: number of channels to use to compute distance metrics (eg 4)
% plotThis: boolean, whether to plot the mahalobnis distance between spikes
%   of thisUnit and otherUnits on the nChansToUse closest channels
% ------
% Outputs
% ------
% isoD: isolation distance
% Lratio: l-ratio
% silhouetteScore: silhouette score
% d2_mahal QQ describe
% Xplot QQ describe
% Yplot QQ describe
%
% based on functions in https://github.com/cortex-lab/sortingQuality

nPCs = size(pc_features, 2); %should be 3 PCs


% features for this cluster
theseFeatures = reshape(pc_features(spikesIdx, :, 1:nChansToUse), numberSpikes, []);
theseChannels = pc_feature_ind(thisUnit, 1:nChansToUse);

nCount = 1;

% features for all other clusters on 4 first channels of this cluster
uC = unique(pc_feature_ind(:, 1));
theseMahalD = nan(numel(uC), 1);
theseOtherFeaturesInd = zeros(0, size(pc_features, 2), nChansToUse);
theseOtherFeatures = zeros(0, size(pc_features, 2), nChansToUse);
for iOtherUnit = 1:numel(uC)
    thisOtherUnit = uC(iOtherUnit);
    if thisOtherUnit ~= thisUnit
        theseOtherChans = pc_feature_ind(iOtherUnit, :);
        for iChannel = 1:nChansToUse
            if ismember(theseChannels(iChannel), theseOtherChans)
                theseOtherSpikes = allSpikesIdx == thisOtherUnit;
                thisCommonChan = find(theseOtherChans == theseChannels(iChannel), 1);
                theseOtherFeatures(nCount:nCount+sum(theseOtherSpikes)-1, :, iChannel) = ...
                    pc_features(theseOtherSpikes, :, thisCommonChan);
                theseOtherFeaturesInd(nCount:nCount+sum(theseOtherSpikes)-1, :, iChannel) = ...
                    ones(size(pc_features(theseOtherSpikes, :, thisCommonChan), 1), ...
                    size(pc_features(theseOtherSpikes, :, thisCommonChan), 2), 1) * double(thisOtherUnit);
                %thisOtherFeaturesInd = ones(numel(nCount:nCount+sum(theseOtherSpikes)-1),1)*double(thisOtherUnit);
                %theseOtherFeaturesInd = [theseOtherFeaturesInd, thisOtherFeaturesInd];

            end
        end
        if any(ismember(theseChannels(:), theseOtherChans))
            %keyboard;
            nCount = nCount + sum(theseOtherSpikes);
            %if ~isempty(theseOtherFeaturesSingle(iOtherUnit).THIS)
            [r, ~, ~] = ind2sub(size(theseOtherFeaturesInd), find(theseOtherFeaturesInd == double(thisOtherUnit)));
            %find(theseOtherFeaturesInd==double(thisOtherUnit))
            if size(theseFeatures, 1) > size(theseFeatures, 2) && size(r, 1) > size(theseFeatures, 2)
                theseMahalD(iOtherUnit) = nanmean(mahal(reshape(theseOtherFeatures(r, :, :), size(r, 1), nPCs*nChansToUse), theseFeatures));
                thesU(iOtherUnit) = double(thisOtherUnit);
            else
                theseMahalD(iOtherUnit) = NaN;
                thesU(iOtherUnit) = double(thisOtherUnit);
            end
            %end
        end
    end
end

if nCount > numberSpikes && numberSpikes > nChansToUse * nPCs 
    %isolation distance
    halfWayPoint = size(theseFeatures, 1);

    %l-ratio
    theseOtherFeatures = reshape(theseOtherFeatures, size(theseOtherFeatures, 1), []);
    mahalD = mahal(theseOtherFeatures, theseFeatures); %squared mahal. distance
    mahalD = sort(mahalD);
    isoD = mahalD(halfWayPoint);
    L = sum(1-chi2cdf(mahalD, nPCs*nChansToUse)); % assumes a chi square distribution - QQ add test for multivariate data
    Lratio = L / numberSpikes;

    %silhouette score
    closestCluster = thesU(find(theseMahalD == min(theseMahalD)));
    mahalDself = mahal(theseFeatures, theseFeatures); %squared mahal. distance
    [r, ~, ~] = ind2sub(size(theseOtherFeaturesInd), find(theseOtherFeaturesInd == double(closestCluster)));
    % find(theseOtherFeaturesInd==double(thisOtherUnit))

    mahalDclosest = mahal(reshape(theseOtherFeatures(r, :, :), size(r, 1), nPCs*nChansToUse), theseFeatures); %squared mahal. distance
    silhouetteScore = (nanmean(mahalDclosest) - nanmean(mahalDself)) / max([mahalDclosest; mahalDself]);


elseif ~isempty(theseOtherFeatures) && numberSpikes > nChansToUse * nPCs %isolation distance not defined
    %isolation distance
    halfWayPoint = NaN;
    isoD = NaN;

    %l-ratio
    theseOtherFeatures = reshape(theseOtherFeatures, size(theseOtherFeatures, 1), []);
    mahalD = mahal(theseOtherFeatures, theseFeatures); %squared mahal. distance
    mahalD = sort(mahalD);
    L = sum(1-chi2cdf(mahalD, nPCs*nChansToUse)); % assumes a chi square distribution - QQ add test for multivariate data
    Lratio = L / numberSpikes;

    %silhouette score
    closestCluster = thesU(find(theseMahalD == min(theseMahalD)));
    mahalDself = mahal(theseFeatures, theseFeatures); %squared mahal. distance
    [r, ~, ~] = ind2sub(size(theseOtherFeaturesInd), find(theseOtherFeaturesInd == double(closestCluster)));
    % find(theseOtherFeaturesInd==double(thisOtherUnit))

    mahalDclosest = mahal(reshape(theseOtherFeatures(r, :, :), size(r, 1), nPCs*nChansToUse), theseFeatures);
    silhouetteScore = (nanmean(mahalDself) - nanmean(mahalDclosest)) / max(mahalDclosest);


else
    %isolation distance
    halfWayPoint = NaN;
    mahalD = NaN;
    isoD = NaN;

    %l-ratio
    L = NaN;
    Lratio = NaN;

    %silhouette score
    silhouetteScore = NaN;

end

if numberSpikes > nChansToUse * nPCs && exist('r', 'var')
    Yplot = reshape(theseOtherFeatures(r, :, :), size(r, 1), nPCs*nChansToUse);
    Xplot = theseFeatures;
    d2_mahal = mahal(Yplot, Xplot);

    mahal_noise = mahal(Xplot,Xplot);
    mahal_unit = mahal(Yplot,Yplot);

    if plotThis

        figure();
        % Calculate the number of subplots needed
        nDims = 12;
        nSubplots = nDims * (nDims - 1) / 2; % n choose 2 for combinations

        % Counter for subplot indexing
        subplotIdx = 1;

        for iDimX = 1:nDims
            for iDimY = iDimX + 1:nDims % Start from iDimX+1 to avoid diagonal and ensure unique pairs
                subplot(6, 11, subplotIdx); % Arrange in 6x11 grid to fit all combinations
                scatter(Xplot(:, iDimX), Xplot(:, iDimY), 10, [0.7, 0.7, 0.7], 'x'); % Xplot scatter plot
                hold on;
                scatter(Yplot(:, iDimX), Yplot(:, iDimY), 10, d2_mahal, 'o', 'filled'); % Yplot scatter plot with Mahalanobis distance coloring
                % Configure colorbar
                if iDimX==1 && iDimY ==2
                hb = colorbar;
                hb.Color = [0.7, 0.7, 0.7];
                ylabel(hb, 'Mahalanobis Distance');
                end

                % Labels and titles for clarity
                xlabel(['PC',  num2str(iDimX -(floor(iDimX/4)*4)+1), ', channel ' num2str(ceil(iDimX/4))]);
                ylabel(['PC',  num2str(iDimY -(floor(iDimY/4)*4)+1), ', channel ' num2str(ceil(iDimY/4))]);

                % Increment subplot index
                subplotIdx = subplotIdx + 1;

                hold off; % Ready for next plot
            end
        end

        figure(); 
        subplot(3,1,1)
        hold on;
        
        % Adjusting histograms to be outlines with a larger bin size for smoothing
        histogram(mahal_noise, 'BinWidth', 2, 'FaceColor', 'none', 'EdgeColor', 'blue', 'DisplayStyle', 'stairs', 'LineWidth', 2);
        histogram(mahal_unit, 'BinWidth', 2, 'FaceColor', 'none', 'EdgeColor', 'red', 'DisplayStyle', 'stairs', 'LineWidth', 2);
        
        legend('Noise (other units)', 'This unit');
        xlabel('Squared Mahalanobis distance');
        ylabel('Count');
        
        hold off;

        subplot(3,1,2)
        [Fx, Xx] = ecdf(mahal_noise);
        [Fy, Xy] = ecdf(mahal_unit);
        
        % Plot cumulative distributions
        plot(Xx, Fx, 'b', 'LineWidth', 2); % Xplot in blue
        plot(Xy, Fy, 'r', 'LineWidth', 2); % Yplot in red
        
        % Find intersection
        [~, intersectIdx] = min(abs(Fx - Fy)); % Simple approach, might need refinement
        xIntersection = Xx(intersectIdx);
        yIntersection = Fx(intersectIdx);
        
        % Plot intersection
        plot(xIntersection, yIntersection, 'ko', 'MarkerFaceColor', 'k'); % Intersection point
        plot([xIntersection, xIntersection], [0, yIntersection], 'k--'); % Dashed line

legend('Xplot CDF', 'Yplot CDF', 'Intersection');
xlabel('Value');
ylabel('Cumulative Probability');
title('Cumulative Distributions of Xplot and Yplot');


        % Add legends and other global annotations outside the loop, if they apply globally
        legend({'this cluster''s spikes', 'nearby clusters'' spikes', ['isolation distance = ', num2str(isoD)], ...
            ['silhouette score = ', num2str(silhouetteScore)], ['l-ratio = ', num2str(Lratio)]}, 'Location', 'bestoutside', 'TextColor', [0.7, 0.7, 0.7]);

        % Apply plot beautification if available
        if exist('prettify_plot', 'file')
            prettify_plot('FigureColor', 'w');
        else
            warning('https://github.com/Julie-Fabre/prettify-matlab repo missing - download it and add it to your matlab path to make plots pretty');
            % Optionally, define makepretty function or similar fallback here
        end

    end
else
    d2_mahal = NaN;
    Yplot = NaN;
    Xplot = NaN;
end
end