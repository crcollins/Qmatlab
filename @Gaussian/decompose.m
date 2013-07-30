function values = decompose(obj, atomLists, orbital)
    % Calculate the percent contribution of a list of atoms to a given
    % atomLists a cell array of arrays of atoms on each fragment 
    % orbital the orbital to do the decomposition on
    ranges = cell(length(atomLists),1);
    for i = 1:length(atomLists)
        ranges{i} = ismember(obj.atom,atomLists{i});
    end
    
    S = obj.overlap;
    pop = zeros(length(ranges));
    for j=1:length(ranges)
        for k = 1:length(ranges)
            pop(j,k) = obj.orb(ranges{j},orbital)' * S(ranges{j},ranges{k}) * obj.orb(ranges{k},orbital);
        end
    end
    values = diag(pop) + (sum(sum(pop)) - sum(diag(pop))) / length(ranges);
end