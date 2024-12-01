function [thigh,knee] = IK_for_serialLeg(linksLength,eqLength,eqAngle)

thigh = eqAngle + acos((linksLength(1).^2+eqLength.^2-linksLength(2).^2)/(2*linksLength(1)*eqLength));
delx = eqLength*cos(eqAngle)-linksLength(1)*cos(thigh);
dely = eqLength*sin(eqAngle)-linksLength(1)*sin(thigh);
knee = atan2(dely,delx);

end