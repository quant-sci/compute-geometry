function x = FaultsCnds(Node,FL)
  a = zeros(size(Node,1),1); v = zeros(size(FL,1),1);
  for i = 1:size(Node,1)
    for j = 1:size(FL,1)
      v(j) = sqrt(((Node(i,1)-FL(j,1)).^2) + ((Node(i,2)-FL(j,2)).^2));
    end
    a(i) = min(v);
  end
  FaultNodes = find(a<0.01); FixedNodes = [FaultNodes];
  SuppFault = zeros(length(FixedNodes),3); 
  SuppFault(:,1)=FixedNodes; SuppFault(1:end-1,2)=1; SuppFault(end,3)=1;
  x = SuppFault;
end