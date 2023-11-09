function prettify_legend(ax)
    yRange = ax.YLim(2) - ax.YLim(1);
    offset = yRange * 0.05;  % This is 5% of the y-axis range; adjust as needed for a larger/smaller offset
    
    handles = findall(ax);

    lines = handles(ishandle(handles) & strcmp(get(handles, 'Type'), 'line'));
    points = handles(ishandle(handles) & strcmp(get(handles, 'Type'), 'scatter'));

    %remove any points
    for h = lines'
        keepLine = length(h.XData) > 1;
    end
    lines = lines(keepLine);

    for h = lines'
        % Extract color and display name
        color = h.Color;
        name = h.DisplayName;

        % Identify the position with maximum separation from other lines
        xdata = h.XData;
        ydata = h.YData;
        separation = zeros(size(ydata));
        
        for i = length(ydata):-1:1
            for hl = lines'
                if hl ~= h
                    if size(xdata) > 1
                        separation(i) = separation(i) + abs(ydata(i) - interp1(hl.XData, hl.YData, xdata(i), 'linear', 'extrap'));
                    else
                        separation(i) =0;
                    end
                end
            end
        end

        if separation == 0
            [~, idx] = max(ydata);
        else
            [~, idx] = max(separation(end-5:end)); % considering last 6 points for better labeling
            idx = idx + length(ydata) - 6;
        end

       

        % Adjust label position based on the position of other lines
        
            yOffset = sign(median([h.YData(idx) - ax.YLim(1), ax.YLim(2) - h.YData(idx)])) * offset;
       
        
        % Place the text
        text(xdata(idx), ydata(idx) + yOffset, name, 'Color', color, 'FontWeight', 'bold');
    end
    
    for h = points'
        % Extract color and display name
        color = h.CData;
        name = h.DisplayName;
        xdata = h.XData;
        ydata = h.YData;

        % Place the text
        text(xdata, ydata + offset, name, 'Color', color, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end

    % remove legend
end
