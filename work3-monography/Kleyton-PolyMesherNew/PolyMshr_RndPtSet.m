function [P,F,ndf] = PolyMshr_RndPtSet(NElem,Domain,Wells,weight,Faults,Faults_Description,Parameters)
%nw is the number of elements desired around the well or wells. It is chosen by the user. It must be less than NElem and, if there is any fault in the region, ndf and naf must be considered (see PolyMshr_FaultGen).
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem
  Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
  Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);                 %Index of seeds inside the domain
  NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
  W(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  if isempty(Wells) == 0
    nw = 100;
    if Ctr < nw
      [W,NumAdded] = PolyMshr_WellsGen(Ctr,W,Wells,weight);
    end
  end
  Ctr = Ctr+NumAdded;
end
for i = 1:length(W)
  P(i,1) = W(i,1);
  P(i,2) = W(i,2);
end
if isempty(Faults) == 0
  [F,ndf] = PolyMshr_FaultGen(Domain,Faults,Faults_Description,Parameters);
  for i = 1:length(F)
    P(length(P)-length(F)+i,1) = F(i,1);
    P(length(P)-length(F)+i,2) = F(i,2);
  end
else
  F = [];
end
end