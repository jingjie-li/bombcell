function prettify_legend(ax, LegendReplace, LegendLocation, LegendBox)


if ~LegendReplace 
    set(ax.Legend, 'Location', LegendLocation)
    set(ax.Legend, 'Box', LegendBox)

    % based on location, order the legend in the same order the lines/ dots
    % are in 
    xvalues = [ax.Legend.Position(1) .* ax.XLim(2) + ax.XLim(1),...
        (ax.Legend.Position(1) +  ax.Legend.Position(3)) .* ax.XLim(2) + ax.XLim(1)];

    handles = findall(ax.Children);

    for iHandle = 1:length(handles)

        h = handles(iHandle);
        % Extract color and display name
        if strcmp(get(h, 'Type'), 'line')
            color = h.Color;
        elseif strcmp(get(h, 'Type'), 'scatter')
            color = h.CData;
        end
        name = h.DisplayName;

        % Calculate label position based on average y-value
        ydata = h.YData;
        if contains(ax.Legend.Location, 'eastoutside')
            yavg(iHandle) = mean(ydata(h.XData(end)));
        elseif contains(ax.Legend.Location, 'westoutside')
            yavg(iHandle) = mean(ydata(h.XData(1)));
        elseif contains(ax.Legend.Location, 'bestoutside') && contains(ax.Legend.Orientation, 'vertical')
            yavg(iHandle) = mean(ydata(h.XData >= xvalues(end)));
        else
            yavg(iHandle) = mean(ydata(h.XData <= xvalues(1)));
        end
        
    end

    [~, order] = sort(yavg, 'descend');
    legend(ax, handles(order)); 



else
    yRange = ax.YLim(2) - ax.YLim(1);
    offset = yRange * 0.05;  % Adjust offset as needed for better spacing

    handles = findall(ax.Children);
    
    lines = handles(ishandle(handles) & strcmp(get(handles, 'Type'), 'line'));
    points = handles(ishandle(handles) & strcmp(get(handles, 'Type'), 'scatter'));

    % Remove any single points from legend
    lines = lines(arrayfun(@(x) numel(x.XData) > 1, lines));

    for h = lines'
        % Extract color and display name
        color = h.Color;
        name = h.DisplayName;

        % Calculate label position based on average y-value
        ydata = h.YData;
        yavg = mean(ydata);

        % Place the text
        text(ax, max(h.XData), yavg + offset, name, 'Color', color, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    end

    for h = points'
        % Extract color and display name
        color = h.CData;
        name = h.DisplayName;

        % Calculate label position based on average y-value
        ydata = h.YData;
        yavg = mean(ydata);

        % Place the text
        text(h.XData, yavg + offset, name, 'Color', color, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end

    % Remove legend
    ax.Legend.Visible = 'off';
end
end
