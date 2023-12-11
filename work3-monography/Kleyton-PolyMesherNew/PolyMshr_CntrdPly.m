function [Pc,A] = PolyMshr_CntrdPly(Domain,Element,Node,NElem,P)
Pc=zeros(NElem,2); A=zeros(NElem,1); BdBox = Domain('BdBox');
for el = 1:NElem
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(Element{el}); 
  vxS=vx([2:nv 1]); vyS=vy([2:nv 1]); %Shifted vertices
  temp = vx.*vyS - vy.*vxS;
  A(el) = 0.5*sum(temp);
  Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end
for i = 1:length(Pc)
  if Pc(i,1) >= BdBox(2) || Pc(i,1) <= BdBox(1)
    Pc(i,1) = P(i,1);
  elseif Pc(i,2) >= BdBox(4) || Pc(i,2) <= BdBox(3)
    Pc(i,2) = P(i,2);
  end
end
for i = 1:length(Pc)
  if isnan(Pc(i,1))
    Pc(i,1) = P(i,1);
    Pc(i,2) = P(i,2);
  end
end
for j = 1:length(A)
  if isnan(A(j))
    A(j) = 0;
  end
end
end