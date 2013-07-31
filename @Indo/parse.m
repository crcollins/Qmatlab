function parse(obj)

% Will have utils.findText() issue errors for us as appropriate
issueErrors = true;
ido_file = [obj.dataPath, obj.filename, '.ido'];
fid1 = fopen(ido_file,'r');

if fid1 == -1
    error(['in Indo.readIndo, could not find file: ',ido_file]);
    fclose('all');
end

t1 = textscan(fid1,'%s');
fclose(fid1);

obj.norb = fread(fid1,1,'integer*4');
obj.aorbAtom = fread(fid1,[1,obj.norb],'integer*4');
obj.aorbAtom = obj.aorbAtom +1; % C++ starts count at 0 instead of 1
obj.aorbType = fread(fid1,[1,obj.norb],'integer*4');

ntest = fread(fid1,1,'integer*4');
if (ntest ~= obj.norb)
    error('atomic and fock basis sizes differ');
end

obj.nfilled = fread(fid1,1,'integer*4');
obj.hfE  = fread(fid1,1,'real*8');
obj.orbE = fread(fid1,[1,obj.norb],'real*8');
obj.orb = fread(fid1,[obj.norb,obj.norb],'real*8');

obj.nsci = fread(fid1,1,'integer*4');
obj.nscibasis = fread(fid1,1,'integer*4');
obj.esci = fread(fid1,[1,obj.nsci],'real*8');
ntest = fread(fid1,[1,2],'integer*4');
obj.r = zeros(obj.nsci,obj.nsci,3);
obj.r(:,:,1) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
ntest = fread(fid1,[1,2],'integer*4');
obj.r(:,:,2) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
ntest = fread(fid1,[1,2],'integer*4');
obj.r(:,:,3) = fread(fid1,[obj.nsci,obj.nsci],'real*8');
temp = fread(fid1,[2,obj.nscibasis],'integer*4');
obj.ehsci = temp' +1; % +1 fixes the counting from 0 issue
obj.wfsci = fread(fid1,[obj.nscibasis,obj.nsci],'real*8');

end