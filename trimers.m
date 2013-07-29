qmatlab = pwd;
mols = {
    '2OAA', {34:44, 1:22, 23:33}, {[35, 21], [22, 27]} ; ...
    '2OAAAO', {34:44, 1:22, 23:33}, {[35, 21], [22, 27]} ; ...
    '2OAACO', {34:44, 1:22, 23:33}, {[35, 21], [22, 27]} ; ...
    '2OAOAO', {38:52, 1:22, 23:37}, {[38, 22], [21 23]} ; ...
    '2OAOAOAO', {38:52, 1:22, 23:37}, {[38, 22], [21, 23]} ; ...
    '2OAOAOCO', {38:52, 1:22, 23:37}, {[38, 22], [21, 23]} ; ...
    '2OBB', {31:38, 1:22, 23:30}, {[31, 21], [22, 23]} ; ...
    '2OBBAO', {31:38, 1:22, 23:30}, {[31, 21], [22, 23]} ; ...
    '2OBBCO', {31:38, 1:22, 23:30}, {[31, 21], [22, 23]} ; ...
    '2OCC',  {31:38, 1:22, 23:30}, {[31, 22], [21, 23]} ; ...
    '2OCCAO', {31:38, 1:22, 23:30}, {[31, 22], [21, 23]} ; ...
    '2OCCCO', {31:38, 1:22, 23:30}, {[31, 22], [21, 23]} ; ...
};

objs = cell(size(mols,1), 1);
vals = zeros(3,size(objs,1));
dataPath = 'C:\Users\ccollins\Downloads\documents-export-2013-07-09\test\';
for i=1:size(mols,1);
    disp(mols{i, 1});
    gstart = Gaussian(dataPath,mols{i,1},struct);
    gstart.run();
    objs{i} = gstart;

    figure(i);
    scale = [1 1];
    sy = 1.1;
    yoffset = [-1.5 -.5 .5 1.5] * sy;

    objs{i}.reorient();
    bb = objs{i}.boundingBox();
    homo = ceil(objs{i}.Nelectrons/2);
    for j = -1:2
        center(1) = -bb.minx - (bb.width/2);
        center(2) = -bb.miny - (bb.height/2) + bb.height * (yoffset(j+2));
        objs{i}.drawStructureOrb(i, homo+j, center, scale);
        values = objs{i}.decompose(mols{i,2}, homo+j);

        text(center(1)-bb.width/3, center(2)+.25*bb.height, sprintf('%.2f',values(1)), 'horizontalalignment', 'right');
        text(center(1), center(2)+.25*bb.height, sprintf('%.2f',values(2)), 'horizontalalignment', 'center');
        text(center(1)+bb.width/3, center(2)+.25*bb.height, sprintf('%.2f',values(3)), 'horizontalalignment', 'left');
    end
end
