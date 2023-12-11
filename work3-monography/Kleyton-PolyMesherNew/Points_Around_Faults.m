FA=zeros(NElem/4,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem/4
  Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem/4,1)+BdBox(1);
  Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem/4,1)+BdBox(3);
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);                   %Index of seeds inside the domain
  NumAdded = min(NElem/4-Ctr,length(I));  %Number of seeds that can be added
  FA(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  D = zeros(NElem/4,1);
  for i = (Ctr+1):length(FA)
    t = FA(i,1);
    aux = 0;
    if t <= x1
      if t >= x0
        Err = 1;
        while Err > 1e-4
          t = FA(i,1) - df(t)*(f(t) - FA(i,2));
          Err = t - aux;
          aux = t;
        end
      else
        t = x0;
      end
    else
      t = x1;
    end
    D(i,1) = sqrt(((FA(i,1)-t)^2) + ((FA(i,2)-f(t))^2));
  end
  logic = D<0.2;
  FA = logic.*FA;
  aux = 1;
  for i = 1:size(FA,1)
    if FA(i,1) > 0
      AP(aux,1) = FA(i,1);
      AP(aux,2) = FA(i,2);
      aux = aux + 1;
    end
  end
  FA = AP;
  FA = table2array(unique(table(FA)));
  NumAdded = size(FA,1) - Ctr;
  Ctr = Ctr + NumAdded;
end
for i = 1:length(FA)
  P(i+NElem/2,1) = FA(i,1);
  P(i+NElem/2,2) = FA(i,2);
end