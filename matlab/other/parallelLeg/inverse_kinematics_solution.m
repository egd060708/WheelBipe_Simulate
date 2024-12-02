function [Alpha1,Alpha2] = inverse_kinematics_solution(x,z)
[L1,L2,L3,L4,L5]=deal(0.15,0.288,0.288,0.15,0.15);
if (x+L5/2)<0
    Alpha2=pi-atan(-z/(-x-L5/2))+acos((L4^2-L3^2+(-0.5*L5-x)^2+z^2)/(2*L1*sqrt((-0.5*L5-x)^2+z^2)));
else
    Alpha2=atan(-z/(x+L5/2))+acos((L4^2-L3^2+(-0.5*L5-x)^2+z^2)/(2*L1*sqrt((-0.5*L5-x)^2+z^2)));
end

if (L5/2-x)<0
    Alpha1=atan(-z/(x-L5/2))-acos((L1^2-L2^2+(0.5*L5-x)^2+z^2)/(2*L1*sqrt((0.5*L5-x)^2+z^2)));
else
    Alpha1=pi-atan(-z/(L5/2-x))-acos((L1^2-L2^2+(0.5*L5-x)^2+z^2)/(2*L1*sqrt((0.5*L5-x)^2+z^2)));
end

end