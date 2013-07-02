function initialize(obj)

% Calculation on the full molecule
fullList = [];
for i = 1:length(obj.fragList)
    fullList = [fullList; obj.fragList{i}(:)];
end

keywords = ['td ', obj.keywords];
name = obj.fullIn.filename;
tempdir = obj.fullIn.writeTPL(name,fullList,keywords);
obj.full = Gaussian(tempdir,obj.fullIn.filename,struct);
obj.full.run();

for i = 1:(length(obj.fragList)-1)
    % fragment calcs
    %  vector pointing from link1 to link2
    direction = obj.fullIn.rcart(:,obj.links{i}(2)) - ...
       obj.fullIn.rcart(:,obj.links{i}(1));
    % hydrogen on frag1 is directed away from link 1 along direction
    rLink{i}{1} = obj.fullIn.rcart(:,obj.links{i}(1)) + ...
       obj.rlinks{i}(1) * direction/norm(direction);
    % hydrogen on frag2 is directed away from link 2 along -direction
    rLink{i}{2} = obj.fullIn.rcart(:,obj.links{i}(2)) - ...
       obj.rlinks{i}(2) * direction/norm(direction);
end

for ifrag = 1:length(obj.fragList)
   name = [obj.fullIn.filename,'-',int2str(ifrag)];
   if ifrag == 1
       use = {rLink{1}{1}};
   else
       if ifrag == length(obj.fragList)
           use = {rLink{ifrag-1}{2}};
       else
           use = {rLink{ifrag-1}{2}, rLink{ifrag}{1}};
       end
   end
   tempdir = obj.full.writeTPL(name,obj.fragList{ifrag},obj.keywords,use);
   obj.frags{ifrag} =  Gaussian(tempdir,name,struct);
   obj.frags{ifrag}.run();
end

% map fragment AO basis to full AO basis
icount = 1;
for ifrag = 1:length(obj.fragList)
   ft = obj.frags{ifrag};
   if ifrag ~= length(obj.fragList) && ifrag ~= 1
       obj.nonLink{ifrag} = find(ft.atom ~= ft.atom(end) & ft.atom ~= ft.atom(end-1));
   else
       obj.nonLink{ifrag} = find(ft.atom ~= ft.atom(end));
   end
   lengthNonLink = length(obj.nonLink{ifrag});
   % assume same ordering of orbs in full as frag, without link
   obj.maps{ifrag} = icount:(icount + lengthNonLink - 1);
   icount = icount + lengthNonLink;
end

% calculate overlaps
for ifrag = 1:length(obj.fragList)
   % frag(ao,mo)' * S(ao,ao) * full(ao,mo)
   noLink = obj.nonLink{ifrag};
   ofrag = obj.frags{ifrag}.orb(noLink,:);
   Stemp = obj.full.overlap(obj.maps{ifrag},:);
   ofull = obj.full.orb(:,:);
   obj.overlap{ifrag} = ofrag' * Stemp * ofull;
end
end