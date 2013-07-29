function drawEnergyLevels(obj, figNum, xs, piOnly)
    % xs the values of the left and right x coords of the line
    % piOnly a boolean value to show pi orbitals only
    %           if piOnly is set to true all the non pi orbitals will be
    %           yellow. Otherwise all the orbitals are colored the same
    %           based on whether or not they are greater than or equal to
    %           the HOMO
    figure(figNum);
    hold on;
    for iorb = 1:size(obj.orb,1)
      e = obj.Eorb(iorb);
      if (iorb <= obj.Nelectrons/2)
         format = 'b'; % filled orbs are blue
      else
         format = 'g'; % empty orbs are green
      end
      if (piOnly && (obj.piCharacter(iorb) < 0.1))
         format = 'y'; % non-pi will be yellow
      end
      plot(xs, [e e], format);
    end
end

