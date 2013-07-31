function parameters(obj, fileprefix)

fields = obj.config.field;
charge = obj.config.charge;
nstates= obj.config.nstates;
norbs = obj.config.norbs;

fid0 = fopen([obj.dataPath,'parameters.txt'],'wt');

fprintf(fid0, 'jobname = %s\n', fileprefix);
fprintf(fid0, 'charge = %i\n', charge);
fprintf(fid0, 'norbs = %i\n', norbs);
fprintf(fid0, 'nstates = %i\n', nstates);
fprintf(fid0, 'efieldx = %8.3f\n', fields(1));
fprintf(fid0, 'efieldy = %8.3f\n', fields(2));
fprintf(fid0, 'efieldz = %8.3f\n', fields(3));

if (~isempty(obj.config.potfile))
    potfile = obj.config.potfile
    fprintf(fid0,'pot_file = %s\n', potfile);
end

fclose(fid0);

end