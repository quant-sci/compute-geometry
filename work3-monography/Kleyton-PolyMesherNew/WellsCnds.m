function x = WellsCnds(Node,Wells)
  WellNodes = zeros(size(Wells,1),1);
  for i = 1:size(Wells,1)
    v = sqrt(((Node(:,1)-Wells(i,1)).^2) + ((Node(:,2)-Wells(i,2)).^2));
    WellNodes(i,1) = find(v==min(v));
  end
  FixedNodes = [WellNodes];
  SuppWell = zeros(length(FixedNodes),3); 
  SuppWell(:,1)=FixedNodes; SuppWell(1:end-1,2)=1; SuppWell(end,3)=1;
  x = SuppWell;
end