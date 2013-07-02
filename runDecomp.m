qmatlab = pwd;

mols = {    ...
%     '1A', {1:13, 14:24}, {[13 14]}; ...
%     '1B', {1:13, 14:21}, {[13 14]}; ...
%     '1C', {1:13, 14:21}, {[13 14]}; ...
%     '2A', {1:15, 16:26}, {[15 16]}; ...
%     '2B', {1:15, 16:23}, {[15 16]}; ...
%     '2C', {1:15, 16:23}, {[15 16]}; ...
%     '2AA', {1:22, 23:33, 34:44}, {[11 30], [8 34]}; ...
    'P3', {1:11, 12:21, 22:32}, {[3 12], [19 22]};
};
objs = cell(size(mols,1), 1);
e = zeros(size(objs,1)*2,1);

dataPath = 'C:\Users\ccollins\Desktop\';
% dataPath = 'C:\Users\ccollins\Desktop\start\ordered\cart\';
for i=1:size(mols,1)
    disp(mols{i, 1});
    gstart = Gaussian(dataPath,mols{i,1},struct);
    gstart.run();
    
    fragList = mols{i,2};
    links = mols{i,3};
    rlinks = {[1.1 1.1], [1.1 1.1]};
    keywords = 'b3lyp/sto-3g';
    obj = Decompose(gstart,fragList,links,rlinks,keywords);
    obj.initialize();
    objs{i} = obj;

    homo = objs{i}.full.Nelectrons/2;
    e(i*2-1) = objs{i}.full.Eorb(homo+6);
    e(i*2) = objs{i}.full.Eorb(homo-5);
end

ymax = max(e);
ymin = min(e);

for i=1:size(objs,1)
    objs{i}.draw(1,i);
    figure(i * 2 - 1);
    axis([0 5 ymin ymax]);
end