function tempDir = writeTPL(obj,jobname,atoms,keywords,rLink)
    if isempty(atoms)
        atoms = 1:length(obj.Z);
    end
    newline = char(10);
    % A collection of possible Atomic symbol - Atomic number pairs
    % TODO replace with a more formalized multi use list
    syms{1} = 'H'; syms{6} = 'C'; syms{7} = 'N'; syms{8} = 'O';
    syms{15} = 'P'; syms{16} = 'S';
    
    tpl_file = [jobname,'.tpl'];
    % temp directory to make all calculations seperate
    tempDir = [tempname('C:\G09W\Scratch'), '\'];
    mkdir(tempDir);
    fid1 = fopen([tempDir,tpl_file],'w');

    fwrite(fid1,['%chk=temp.chk',newline]);
    fwrite(fid1,['# ',keywords,' NoSymmetry iop(3/33=4) pop=regular',newline]);
    fwrite(fid1,newline);
    fwrite(fid1,jobname);
    fwrite(fid1,newline);
    fwrite(fid1,newline);
    fwrite(fid1,sprintf('%d %d', obj.charge, obj.multiplicity));
    fwrite(fid1,newline);
    % write out all of the coords of the atoms
    for iatom = atoms(:)'
       fwrite(fid1,[' ',syms{obj.Z(iatom)},' ']);
       for ic = 1:3
          fwrite(fid1,sprintf('%.4f ',obj.rcart(ic,iatom)));
       end
       fwrite(fid1,newline);
    end
    % Add a final Hydrogen if there is a link defined
    if (nargin > 4)
       fwrite(fid1,[' H ',num2str(rLink(:)'),newline]);
    end
    % A newline is required at the end of the file for gaussian to take it
    fwrite(fid1,newline);

    fclose(fid1);
end

