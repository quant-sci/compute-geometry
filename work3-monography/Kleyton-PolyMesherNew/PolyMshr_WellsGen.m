function [W,NumAdded] = PolyMshr_WellsGen(Ctr,W,Wells,weight)
%Wells is a nx2 matrix of coordinates of the wells on the region in the form of [x1 y1;x2 y2;x3 y3;...], where n is the number of wells;
%weight is a number defined by the user in order to concentrate more or less the polygons around the wells.
NWells = size(Wells,1);
wx = Wells(:,1);
wy = Wells(:,2);
distribute = 0;
for i = 1:NWells
  distribute = distribute + exp(-weight*((((W(:,1)-wx(i))).^2)+((W(:,2)-wy(i)).^2)));
end
r0 = distribute.^2;
logic = r0/max(r0)>0.8;
W = logic.*W;
aux = 1;
for i = 1:size(W,1)
  if W(i,1) > 0
    AP(aux,1) = W(i,1);
    AP(aux,2) = W(i,2);
    aux = aux + 1;
  end
end
W = AP;
Dist = zeros(NWells,size(W,1));
for i = 1:size(W,1)
  for j = 1:NWells
    Dist(j,i) = sqrt(((W(i,1)-Wells(j,1))^2) + ((W(i,2)-Wells(j,2))^2));
  end
end
if NWells == 1
  logic = transpose(Dist<0.001);
else
  logic = transpose(min(Dist)<0.001);
end
W = logic.*W;
aux = 1;
for i = 1:size(W,1)
  if W(i,1) > 0
    AP(aux,1) = W(i,1);
    AP(aux,2) = W(i,2);
    aux = aux + 1;
  end
end
W = AP;
W = table2array(unique(table(W)));
NumAdded = size(W,1) - Ctr;
end