function L= plucker_calculation(point1, point2)

p01 = point1(1)*point2(2) - point2(1)*point1(2); %w1x2 - w2x1
p02 = point1(1)*point2(3) - point2(1)*point1(3);
p03 = point1(1)*point2(4) - point2(1)*point1(4);
p23 = point1(3)*point2(4) - point2(4)*point1(3);
p31 = point1(4)*point2(2) - point2(4)*point1(2);
p12 = point1(2)*point2(3) - point2(3)*point1(2);
L = [p01, p02, p03, p23, p31, p12];
end