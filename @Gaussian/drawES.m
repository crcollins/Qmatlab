function drawES(obj,figNum,width,xoffset)
    % width the width of the line ot draw for the energy level
    % xoffset the offset of where the line is located
    if (nargin < 2)
       figNum = 1;
    end
    figure(figNum); %should match figure number from draw.m
    nlevels = length(obj.Ees);
    spacing = width/nlevels;
    for i = 1:nlevels
        x = (i-.5)*spacing+xoffset;
        ncomp = length(obj.Ecomp{i});
        t = 0;
        spacing2 = width / max([35 ncomp]); % keep clusters together
        for ic = 1:ncomp
            t1 = obj.Ecomp{i}{ic};
            if (t1.amp.^2 > 0.2)     
                hold on;
                plot([x+t*spacing2 x+t*spacing2],[obj.Eorb(t1.filled) obj.Eorb(t1.empty)],'r-','lineWidth',10*t1.amp.^2);
                t = t + 1;
            end
        end
    end
end



