function [Pos,IM] = IM_for_serialLeg_fact(Llinks,Mlinks,Ilinks,eqLength,eqAngle)

Pos = zeros(2);
[knee,thigh] = IK_for_serialLeg1_fact(Llinks,eqLength,eqAngle);

PosLinks = zeros(8,2);
PosLinks(1,1) = -Llinks(1)*0.5*sin(thigh);
PosLinks(1,2) = -Llinks(1)*0.5*cos(thigh);
PosLinks(2,1) = -Llinks(1)*sin(thigh) - 0.5*Llinks(2)*sin(knee);
PosLinks(2,2) = -Llinks(1)*cos(thigh) - 0.5*Llinks(2)*cos(knee);

PosLinks(3,1) = Llinks(3)*0.5*sin(thigh);
PosLinks(3,2) = -Llinks(3)*0.5*cos(thigh);
PosLinks(4,1) = Llinks(3)*sin(thigh) + 0.5*Llinks(4)*sin(knee);
PosLinks(4,2) = -Llinks(3)*cos(thigh) - 0.5*Llinks(3)*cos(knee);

PosLinks(5,1) = -Llinks(3)*sin(thigh);
PosLinks(5,2) = -Llinks(3)*cos(thigh);
PosLinks(6,1) = -(Llinks(3) + Llinks(6)/2)*sin(thigh) + (Llinks(5) - Llinks(4))*sin(knee);
PosLinks(6,2) = -(Llinks(3) + Llinks(6)/2)*cos(thigh) + (Llinks(5) - Llinks(4))*cos(knee) ;

PosLinks(7,1) = -Llinks(1)*sin(thigh);
PosLinks(7,2) = -Llinks(1)*cos(thigh);
PosLinks(8,1) = -Llinks(1)*sin(thigh) - Llinks(2)*sin(knee);
PosLinks(8,2) = -Llinks(1)*cos(thigh) - Llinks(2)*cos(knee);


Mall = Mlinks(1) + Mlinks(2) + Mlinks(3) + Mlinks(4) + Mlinks(5) + Mlinks(6) + Mlinks(7) + Mlinks(8);
Pos(1) = (PosLinks(1,1)*Mlinks(1) + PosLinks(2,1)*Mlinks(2) + PosLinks(3,1)*Mlinks(3) + PosLinks(4,1)*Mlinks(4) + PosLinks(5,1)*Mlinks(5) + PosLinks(6,1)*Mlinks(6) + PosLinks(7,1)*Mlinks(7) + PosLinks(8,1)*Mlinks(8))/Mall;
Pos(2) = (PosLinks(1,2)*Mlinks(1) + PosLinks(2,2)*Mlinks(2) + PosLinks(3,2)*Mlinks(3) + PosLinks(4,2)*Mlinks(4) + PosLinks(5,2)*Mlinks(5) + PosLinks(6,2)*Mlinks(6) + PosLinks(7,2)*Mlinks(7) + PosLinks(8,2)*Mlinks(8))/Mall;

IM = 0;
% eqX = eqLength*sin(eqAngle);
% eqZ = eqLength*cos(eqAngle);
for i=1:1:8
    IM = IM + Ilinks(i) + Mlinks(i)*((Pos(1)-PosLinks(i,1)).^2+(Pos(2)-PosLinks(i,2)).^2);
end

end
