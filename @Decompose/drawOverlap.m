function drawOverlap(obj, figNum, ifrag, piOnly, xs, threshold)
    % draw overlap lines for a fragment
    % figNum the number of the figure to draw on
    % ifrag the number of the fragment to use
    % piOnly a boolean value to toggle showing the pi orbitals
    %        if this is true then only the pi orbitals will be drawn
    % xs is the left and right x coords for drawing the dotted lines
    % threshold is an optional argument to set the limit on what is counted
    %           as an overlap
    figure(figNum);
    hold on;
    if nargin < 5
        threshold = 0.25;
    end
    
    ol = obj.overlap{ifrag}.^2;
    Efrag = obj.frags{ifrag}.Eorb;
    Efull = obj.full.Eorb;
    for j = 1:length(Efrag)
        for k = 1:length(Efull)
            if (~piOnly || (obj.full.piCharacter(k)>0.1))
                if ol(j,k)>threshold
                    e1 = Efrag(j);
                    e2 = Efull(k);
                    if (k <= obj.full.Nelectrons/2)
                        format = 'b';
                    else
                        format = 'g';
                    end
                    plot(xs,[e1 e2],[format,':']);
                end
            end
        end
    end
end

