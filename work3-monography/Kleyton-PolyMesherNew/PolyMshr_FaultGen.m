function [F,ndf] = PolyMshr_FaultGen(Domain,Faults,Faults_Description,Parameters)
%Faults is a matrix of coordinates which represents the faults inside the domain;
%If Faults_Description = 1, then Fault = [x0 y0;x1 y1] is a straight line iniciating in (x0,y0) and ending in (x1,y1);
%If Faults_Description = 1, then Parameters = [];
%If Faults_Description = 2, then Fault = [x0 x1] is a senoid line iniciating in x0 and ending in x1 in the form y = a*sin(b*x+c)+d;
%If Faults_Description = 3, then Fault = [x0 x1] is a cosenoid line iniciating in x0 and ending in x1 in the form y = a*cos(b*x+c)+d;
%If Faults_Description = 4, then Fault = [x0 x1] is a logarithm line iniciating in x0 and ending in x1 in the form y = (a/c)*log(b,x)+d, where we have logarithm of x on base b;
%If Faults_Description = 5, then Fault = [x0 x1] is a polinomial line iniciating in x0 and ending in x1 in the form y = a*(x^3)+b*(x^2)+c*x+d;
%Parameters = [a b c d] must be filled for Faults_Description > 1;
%ndf is the number of elements to delimitate the fault and it is chosen by the user. The number should be even in order to distribute the same amount of elements on both sides of the fault and must greater than zero; and
%naf is the number of elements around the fault in order to concentrate a few elements for better understanding of the behaviour around the fault and it is chosen by the user. It must be greater than zero.
ndf = 60;
naf = 300;
Fd = zeros(ndf,2);
Fa = zeros(naf,2);
F = zeros(ndf+naf,2);
if Faults_Description == 1
  x0 = Faults(1,1);
  y0 = Faults(1,2);
  x1 = Faults(2,1);
  y1 = Faults(2,2);
  D = sqrt(((x1-x0)^2)+((y1-y0)^2));
  alfa = atan2((y1-y0),(x1-x0));
  dx = 2*D*cos(alfa)/ndf;
  dy = 2*D*sin(alfa)/ndf;
  dp = 0.001;
  for i = 1:ndf/2
    dw = (i-1)*(sqrt((dx^2) + (dy^2)));
    theta = atan2(dp,((i-1)*dw));
    dz = sqrt((dp^2) + (dw^2));
    deltax = dz*cos(alfa + theta);
    deltay = dz*sin(alfa + theta);
    Fd(2*i-1,1) = x0 + deltax;
    Fd(2*i-1,2) = y0 + deltay;
    deltax = dz*cos(alfa - theta);
    deltay = dz*sin(alfa - theta);  
    Fd(2*i,1) = x0 + deltax;
    Fd(2*i,2) = y0 + deltay;
  end
else
  syms f(x) 
  syms df(x)
  if Faults_Description == 2
    a = Parameters(1,1);
    b = Parameters(1,2);
    c = Parameters(1,3);
    d = Parameters(1,4);
    f(x) = a*sin(b*x+c) + d;
    df(x) = a*b*cos(b*x+c);
  end
  if Faults_Description == 3
    a = Parameters(1,1);
    b = Parameters(1,2);
    c = Parameters(1,3);
    d = Parameters(1,4);
    f(x) = a*cos(b*x+c) + d;
    df(x) = -a*b*sin(b*x+c);
  end
  if Faults_Description == 4
    a = Parameters(1,1);
    b = Parameters(1,2);
    c = Parameters(1,3);
    d = Parameters(1,4);
    f(x) = (a/c)*(log(x)/log(b)) + d;
    df(x) = (a/(log(b)*c))*(1/x);
  end
  if Faults_Description == 5
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
end
Ctr=0;
while Ctr<naf
  Y(:,1) = (x1 - x0)*rand(naf,1) + x0;
  if Faults_Description == 1
    Y(:,2) = (y1 - y0)*rand(naf,1) + y0;
  else
    Y(:,2) = (f(x) - f(x0))*rand(naf,1) + f(x0);
  end
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);               %Index of seeds inside the domain
  NumAdded = min(naf-Ctr,length(I));  %Number of seeds that can be added
  Fa(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  D = zeros(naf,1);
  for i = (Ctr+1):length(Fa)
    if Faults_Description == 1
      a = y0 - y1;
      b = x1 - x0;
      c = x0*y1 - x1*y0;
      D(i,1) = (abs(a*Fa(i,1)+b*Fa(i,2)+c))/(sqrt((a^2)+(b^2)));
    else
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