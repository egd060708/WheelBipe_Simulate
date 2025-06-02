function [anglef,angleb] = IK_for_serialLeg1_fact(linksLength,eqLength,eqAngle)
x = eqLength * sin(eqAngle);
z = eqLength * cos(eqAngle);


[L1,L2,L3,L4,L5]=deal(linksLength(1),linksLength(2),linksLength(2),linksLength(1),0);
if (x+L5/2)<0
    angleb=pi-atan(-z/(-x-L5/2))+acos((L4^2-L3^2+(-0.5*L5-x)^2+z^2)/(2*L1*sqrt((-0.5*L5-x)^2+z^2)));
else
    angleb=atan(-z/(x+L5/2))+acos((L4^2-L3^2+(-0.5*L5-x)^2+z^2)/(2*L1*sqrt((-0.5*L5-x)^2+z^2)));
end

if (L5/2-x)<0
    anglef=atan(-z/(x-L5/2))-acos((L1^2-L2^2+(0.5*L5-x)^2+z^2)/(2*L1*sqrt((0.5*L5-x)^2+z^2)));
else
    anglef=pi-atan(-z/(L5/2-x))-acos((L1^2-L2^2+(0.5*L5-x)^2+z^2)/(2*L1*sqrt((0.5*L5-x)^2+z^2)));
end

end