function drawPercents(obj, figNum, orb, center, xoffset)
    figure(figNum);
    hold on;
    
    S = obj.full.overlap;
    for i = 1:length(obj.fragList)
        ranges{i} = obj.fragList{i};
    end
    popFrags = zeros(length(ranges));
    for j=1:length(popFrags)
        for k = 1:length(popFrags)
            popFrags(j,k) = obj.full.orb(ranges{j},orb)' * S(ranges{j},ranges{k}) * obj.full.orb(ranges{k},orb);
        end
    end
    left = popFrags(1,1) + 0.5 * popFrags(1,2) + 0.5 * popFrags(2,1);
    right = popFrags(2,2) + 0.5 * popFrags(1,2) + 0.5 * popFrags(1,2);
    text(center(1)-xoffset, center(2), sprintf('%.2f',left), 'horizontalalignment', 'right');
    text(center(1)+xoffset, center(2), sprintf('%.2f',right), 'horizontalalignment', 'left');
end