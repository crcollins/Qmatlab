function res = piCharacter(obj, iorb)
    % select only the pi orbitals on the Z axis
    a1 = obj.orb((obj.type == 1) & (obj.subtype == 3) , iorb);
    res = sum(a1.^2);
end