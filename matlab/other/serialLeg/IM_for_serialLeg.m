function [Pos,IM] = IM_for_serialLeg(Llinks,Mlinks,Ilinks,eqLength,eqAngle)

Pos = zeros(2);
[knee,thigh] = IK_for_serialLeg1(Llinks,eqLength,eqAngle);

PosLinks = zeros(4,2);
PosLinks(1,1) = -Llinks(1)*0.5*sin(thigh);
PosLinks(1,2) = -Llinks(1)*0.5*cos(thigh);
PosLinks(2,1) = -Llinks(1)*sin(thigh) - 0.5*Llinks(2)*sin(knee);
PosLinks(2,2) = -Llinks(1)*cos(thigh) - 0.5*Llinks(2)*cos(knee);
PosLinks(3,1) = -Llinks(1)*0.5*sin(thigh) - Llinks(4)*sin(knee);
PosLinks(3,2) = -Llinks(1)*0.5*cos(thigh) - Llinks(4)*cos(knee);
PosLinks(4,1) = -0.5*Llinks(4)*sin(knee);
PosLinks(4,2) = -0.5*Llinks(4)*cos(knee);

Mall = Mlinks(1) + Mlinks(2) + Mlinks(3) + Mlinks(4);
Pos(1) = (PosLinks(1,1)*Mlinks(1) + PosLinks(2,1)*Mlinks(2) + PosLinks(3,1)*Mlinks(3) + PosLinks(4,1)*Mlinks(4))/Mall;
Pos(2) = (PosLinks(1,2)*Mlinks(1) + PosLinks(2,2)*Mlinks(2) + PosLinks(3,2)*Mlinks(3) + PosLinks(4,2)*Mlinks(4))/Mall;

IM = 0;
eqX = eqLength*sin(eqAngle);
eqZ = eqLength*cos(eqAngle);
for i=1:1:4
    IM = IM + Ilinks(i) + Mlinks(i)*((eqX-PosLinks(i,1)).^2+(eqZ-PosLinks(i,2)).^2);
end

end