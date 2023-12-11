%------------------ PolyMesher  version: 1.1 (Aug13) ---------------------%
% Ref1: C Talischi, GH Paulino, A Pereira, IFM Menezes,                   %
%      "PolyMesher: A general-purpose mesh generator for polygonal        %
%      elements written in Matlab", Struct Multidisc Optim, 2012,         %
%      DOI 10.1007/s00158-011-0706-z                                      %
%                                                                         %
% Ref2: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho,      %
%      "Implementation of fluid flow topology optimization in PolyTop",   %
%      Struct Multidisc Optim, 2013, DOI XX.XXXX/XXXXXX-XXX-XXX-X         %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load,P] = PolyMesher(Domain,NElem,MaxIter,Wells,weight,Faults,Faults_Description,Parameters,P)
if ~exist('P','var'), [P,F,FL,ndf]=PolyMshr_RndPtSet(NElem,Domain,Wells,weight,Faults,Faults_Description,Parameters); end
NElem = size(P,1);
Tol=5e-6; It=0; Err=1; c=1.5;
BdBox = Domain('BdBox'); PFix = Domain('PFix');
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
Pc = P; figure;
while(It<=MaxIter && Err>Tol)
  Alpha = c*sqrt(Area/NElem);
  for i = 1:length(P)
    if i > (length(P)-ndf)
      if isempty(F) == 0
        P(i,1) = F(i-length(P)+length(F),1);
        P(i,2) = F(i-length(P)+length(F),2);
      else
        P(i,1) = Pc(i,1); %Lloyd's update
        P(i,2) = Pc(i,2);
      end
    else
      P(i,1) = Pc(i,1); %Lloyd's update
      P(i,2) = Pc(i,2);
    end
  end
  R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha);   %Generate the reflections
  [P,R_P] = PolyMshr_FixedPoints(P,R_P,PFix);   %Fixed Points
  [Node,Element] = voronoin([P;R_P]);           %Construct Voronoi diagram
  for i = 1:length(Node)
    if Node(i,1) < BdBox(1)
      Node(i,1) = BdBox(1);
    elseif Node(i,1) > BdBox(2)
      Node(i,1) = BdBox(2);
    elseif Node(i,2) < BdBox(3)
      Node(i,2) = BdBox(3);
    elseif Node(i,2) > BdBox(4)
      Node(i,2) = BdBox(4);
    end
  end
  [Pc,A] = PolyMshr_CntrdPly(Domain,Element,Node,NElem,P);
  Area = sum(abs(A));
  Err = sqrt(sum((A.^2).*sum((Pc-P).*(Pc-P),2)))*NElem/Area^1.5;
  fprintf('It: %3d   Error: %1.3e\n',It,Err); It=It+1;
  if NElem<=2000, PolyMshr_PlotMsh(Node,Element,NElem); end
end
[Node,Element] = PolyMshr_ExtrNds(NElem,Node,Element);  %Extract node list
[Node,Element] = PolyMshr_CllpsEdgs(Node,Element,0.1);  %Remove small edges
[Node,Element] = PolyMshr_RsqsNds(Node,Element);        %Reorder Nodes
BC=Domain('BC',{Node,Element}); Supp=BC{1}; Load=BC{2}; %Recover BC arrays
if isempty(Wells) == 0
  WC = WellsCnds(Node,Wells);
else
  WC = [];
end
if isempty(Faults) == 0
  FC = FaultCnds(Node,FL);
else
  FC = [];
end
PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load,WC,FC);    %Plot mesh and BCs
end
%------------------------------------------------- GENERATE RANDOM POINTSET
function [P,F,FL,ndf] = PolyMshr_RndPtSet(NElem,Domain,Wells,weight,Faults,Faults_Description,Parameters)
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
    nw = 200;
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
  [F,FL,ndf] = PolyMshr_FaultGen(Domain,Faults,Faults_Description,Parameters);
  for i = 1:length(F)
    P(length(P)-length(F)+i,1) = F(i,1);
    P(length(P)-length(F)+i,2) = F(i,2);
  end
else
  F = [];
end
end
%------------------------------------------------------------- FIXED POINTS
function [P,R_P] = PolyMshr_FixedPoints(P,R_P,PFix)
PP = [P;R_P];
for i = 1:size(PFix,1)
  [B,I] = sort(sqrt((PP(:,1)-PFix(i,1)).^2+(PP(:,2)-PFix(i,2)).^2));
  for j = 2:4
    n = PP(I(j),:) - PFix(i,:); n = n/norm(n);
    PP(I(j),:) = PP(I(j),:)-n*(B(j)-B(1));
  end
end
P = PP(1:size(P,1),:); R_P = PP(1+size(P,1):end,:);
end
%--------------------------------------------------------- REFLECT POINTSET
function R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha)
eps=1e-8; eta=0.9;
d = Domain('Dist',P);  
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain('Dist',P+repmat([eps,0],NElem,1))-d)/eps;
n2 = (Domain('Dist',P+repmat([0,eps],NElem,1))-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,2)
R_P(:,1) = P1(I)-2*n1(I).*d(I);  
R_P(:,2) = P2(I)-2*n2(I).*d(I);
d_R_P = Domain('Dist',R_P);
J = abs(d_R_P(:,end))>=eta*abs(d(I)) & d_R_P(:,end)>0;
R_P=R_P(J,:); R_P=unique(R_P,'rows');
end
%---------------------------------------------- COMPUTE CENTROID OF POLYGON
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
%------------------------------------------------------- EXTRACT MESH NODES
function [Node,Element] = PolyMshr_ExtrNds(NElem,Node0,Element0)
map = unique([Element0{1:NElem}]);
cNode = 1:size(Node0,1);
cNode(setdiff(cNode,map)) = max(map);
[Node,Element] = PolyMshr_RbldLists(Node0,Element0(1:NElem),cNode);
end
%----------------------------------------------------- COLLAPSE SMALL EDGES
function [Node0,Element0] = PolyMshr_CllpsEdgs(Node0,Element0,Tol)
while(true)
  cEdge = [];
  for el=1:size(Element0,1)
    if size(Element0{el},2)<4, continue; end  %Cannot collapse triangles
    vx=Node0(Element0{el},1); vy=Node0(Element0{el},2); nv=length(vx);
    beta = atan2(vy-sum(vy)/nv, vx-sum(vx)/nv);
    beta = mod(beta([2:end 1]) -beta,2*pi);
    betaIdeal = 2*pi/size(Element0{el},2);
    Edge = [Element0{el}',Element0{el}([2:end 1])'];
    cEdge = [cEdge; Edge(beta<Tol*betaIdeal,:)];
  end
  if (size(cEdge,1)==0), break; end
  cEdge = unique(sort(cEdge,2),'rows');
  cNode = 1:size(Node0,1);
  for i=1:size(cEdge,1)
    cNode(cEdge(i,2)) = cNode(cEdge(i,1));
  end
  [Node0,Element0] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
end
%--------------------------------------------------------- RESEQUENSE NODES
function [Node,Element] = PolyMshr_RsqsNds(Node0,Element0)
NNode0=size(Node0,1); NElem0=size(Element0,1);
ElemLnght=cellfun(@length,Element0); nn=sum(ElemLnght.^2); 
i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
for el = 1:NElem0
  eNode=Element0{el}; ElemSet=index+1:index+ElemLnght(el)^2;
  i(ElemSet) = kron(eNode,ones(ElemLnght(el),1))';
  j(ElemSet) = kron(eNode,ones(1,ElemLnght(el)))';
  s(ElemSet) = 1;
  index = index+ElemLnght(el)^2;
end
K = sparse(i,j,s,NNode0, NNode0);
p = symrcm(K);
cNode(p(1:NNode0))=1:NNode0;
[Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
%------------------------------------------------------------ REBUILD LISTS
function [Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode)
Element = cell(size(Element0,1),1);
[foo,ix,jx] = unique(cNode);
if ~isequal(size(jx),size(cNode)), jx=jx'; end % +R2013a compatibility fix
if size(Node0,1)>length(ix), ix(end)=max(cNode); end
Node = Node0(ix,:); 
for el=1:size(Element0,1)
  Element{el} = unique(jx(Element0{el}));
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(vx);
  [foo,iix] = sort(atan2(vy-sum(vy)/nv,vx-sum(vx)/nv));
  Element{el} = Element{el}(iix);
end
end
%---------------------------------------------------------------- PLOT MESH
function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load,WC,FC)
clf; axis equal; axis off; hold on;
Element = Element(1:NElem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
  plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
end
if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
  plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
end
if exist('WC','var')&&~isempty(WC) %Plot Wells BC if specified
  plot(Node(WC(:,1),1),Node(WC(:,1),2),'b*','MarkerSize',8);
end
if exist('FC','var')&&~isempty(FC) %Plot Faults BC if specified
  plot(Node(FC(:,1),1),Node(FC(:,1),2),'r','MarkerSize',8);
end
end
%--------------------------------------------- GENERATE POINTS AROUND WELLS
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
logic = r0/max(r0)>0.9;
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
%-------------------------------------------- GENERATE POINTS AROUND FAULTS
function [F,FL,ndf] = PolyMshr_FaultGen(Domain,Faults,Faults_Description,Parameters)
%Faults is a matrix of coordinates which represents the faults inside the domain;
%If Faults_Description = 1, then Fault = [x0 x1] is a senoid line iniciating in x0 and ending in x1 in the form y = a*sin(b*x+c)+d;
%If Faults_Description = 2, then Fault = [x0 x1] is a cosenoid line iniciating in x0 and ending in x1 in the form y = a*cos(b*x+c)+d;
%If Faults_Description = 3, then Fault = [x0 x1] is a logarithm line iniciating in x0 and ending in x1 in the form y = (a/c)*log(b,x)+d, where we have logarithm of x on base b;
%If Faults_Description = 4, then Fault = [x0 x1] is a polinomial line iniciating in x0 and ending in x1 in the form y = a*(x^3)+b*(x^2)+c*x+d;
%Parameters = [a b c d] must be filled for Faults_Description >= 1;
%ndf is the number of elements to delimitate the fault and it is chosen by the user. The number should be even in order to distribute the same amount of elements on both sides of the fault and must greater than zero; and
%naf is the number of elements around the fault in order to concentrate a few elements for better understanding of the behaviour around the fault and it is chosen by the user. It must be greater than zero.
ndf = 60;
naf = 300;
Fd = zeros(ndf,2);
Fa = zeros(naf,2);
F = zeros(ndf+naf,2);
syms f(x) 
syms df(x)
if Faults_Description == 1
  a = Parameters(1,1);
  b = Parameters(1,2);
  c = Parameters(1,3);
  d = Parameters(1,4);
  f(x) = a*sin(b*x+c) + d;
  df(x) = a*b*cos(b*x+c);
end
if Faults_Description == 2
  a = Parameters(1,1);
  b = Parameters(1,2);
  c = Parameters(1,3);
  d = Parameters(1,4);
  f(x) = a*cos(b*x+c) + d;
  df(x) = -a*b*sin(b*x+c);
end
if Faults_Description == 3
  a = Parameters(1,1);
  b = Parameters(1,2);
  c = Parameters(1,3);
  d = Parameters(1,4);
  f(x) = (a/c)*(log(x)/log(b)) + d;
  df(x) = (a/(log(b)*c))*(1/x);
end
if Faults_Description == 4
  a = Parameters(1,1);
  b = Parameters(1,2);
  c = Parameters(1,3);
  d = Parameters(1,4);
  f(x) = a*(x^3) + b*(x^2) + c*x + d;
  df(x) = 3*a*(x^2) + 2*b*x + c;
end
x0 = Faults(1,1);
x1 = Faults(1,2);
dp = 0.001;
dx = (x1 - x0)/(ndf/2);
x = x0 + dx;
y0 = f(x0);
y = f(x);
for i = 1:ndf/2
  if isnan(df(x))
    if f(x) > f(x-dx) 
      beta = pi/2;   %In this case, the function tendency is to +inf
    else
      beta = 3*pi/2; %In this case, the function tendency is to -inf
    end
  else
    beta = atan2(df(x),1);
  end
  alfa = atan2((y - y0),(x - x0));
  gamma = alfa - beta;
  dw = sqrt(((x - x0)^2) + ((y - y0)^2));
  dz = sqrt((dp^2) + (dw^2) + (2*dp*dw*sin(gamma)));
  theta = asin((dp/dz)*cos(gamma));
  deltax = dz*cos(alfa + theta);
  deltay = dz*sin(alfa + theta);
  Fd(2*i-1,1) = x0 + deltax;
  Fd(2*i-1,2) = y0 + deltay;
  Fd(2*i,1) = x0 + 2*(x-x0) - deltax;
  Fd(2*i,2) = y0 + 2*(y-y0) - deltay;
  x = x + dx;
  y = f(x);
end
FL = zeros((round((x1-x0)/0.01)) + 1,2); x = x0;
for i = 1:(round(((x1-x0)/0.01)+1))
  FL(i,1) = x;
  FL(i,2) = f(x);
  x = x + 0.01;
end
Ctr=0;
while Ctr<naf
  Y(:,1) = (x1 - x0)*rand(naf,1) + x0;
  Y(:,2) = (f(x) - f(x0))*rand(naf,1) + f(x0);
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);               %Index of seeds inside the domain
  NumAdded = min(naf-Ctr,length(I));  %Number of seeds that can be added
  Fa(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  D = zeros(naf,1);
  for i = (Ctr+1):length(Fa)
    t = Fa(i,1);
    aux = 0;
    if t <= x1
      if t >= x0
        Err = 1;
        while Err > 1e-3
          t = Fa(i,1) - df(t)*(f(t) - Fa(i,2));
          Err = t - aux;
          aux = t;
        end
      else
        t = x0;
      end
    else
      t = x1;
    end
    D(i,1) = eval(sqrt(((Fa(i,1)-t)^2) + ((Fa(i,2)-f(t))^2)));
  end
  logic = D<0.1;
  Fa = logic.*Fa;
  aux = 1;
  AP = [];
  for i = 1:size(Fa,1)
    if Fa(i,1) > 0
      AP(aux,1) = Fa(i,1);
      AP(aux,2) = Fa(i,2);
      aux = aux + 1;
    end
  end
  Fa = AP;
  Fa = table2array(unique(table(Fa)));
  NumAdded = size(Fa,1) - Ctr;
  Ctr = Ctr + NumAdded;
end
for i = 1:(naf+ndf)
  if i <= naf
    F(i,1) = Fa(i,1);
    F(i,2) = Fa(i,2);
  else
    F(i,1) = Fd((i-naf),1);
    F(i,2) = Fd((i-naf),2);
  end
end  
end
%---------------------------------------------- SPECIFY WELLS CONDITIONS
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
%---------------------------------------------- SPECIFY FAULTS CONDITIONS
function x = FaultCnds(Node,FL)
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