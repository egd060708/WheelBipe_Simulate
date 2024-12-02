function [anglef,angleb] = IK_for_serialLeg1(linksLength,eqLength,eqAngle)

angleb = eqAngle + acos((linksLength(1).^2+eqLength.^2-linksLength(2).^2)/(2*linksLength(1)*eqLength));
delx = eqLength*cos(eqAngle)-linksLength(1)*cos(angleb);
dely = eqLength*sin(eqAngle)-linksLength(1)*sin(angleb);
anglef = atan2(dely,delx);

end