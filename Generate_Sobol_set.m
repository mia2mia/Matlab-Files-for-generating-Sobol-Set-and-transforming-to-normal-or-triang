%%% Sensitivity size
%{  
Generate Sobol Set (1,000 by n) 10 times
p = sobolset(1000,'Skip',1e3,'Leap',1e2);
p = scramble(p,'MatousekAffineOwen');
X0 = net(p,29122);
%}
%{
p1 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p1 = scramble(p1,'MatousekAffineOwen');
X01 = net(p1,29122);

p2 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p2 = scramble(p2,'MatousekAffineOwen');
X02 = net(p2,29122);

p3 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p3 = scramble(p3,'MatousekAffineOwen');
X03 = net(p3,29122);

p4 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p4 = scramble(p4,'MatousekAffineOwen');
X04 = net(p4,29122);

p5 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p5 = scramble(p5,'MatousekAffineOwen');
X05 = net(p5,29122);

p6 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p6 = scramble(p6,'MatousekAffineOwen');
X06 = net(p6,29122);

p7 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p7 = scramble(p7,'MatousekAffineOwen');
X07 = net(p7,29122);

p8 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p8 = scramble(p8,'MatousekAffineOwen');
X08 = net(p8,29122);

p9 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p9 = scramble(p9,'MatousekAffineOwen');
X09 = net(p9,29122);

p10 = sobolset(1000,'Skip',1e3,'Leap',1e2);
p10 = scramble(p10,'MatousekAffineOwen');
X10 = net(p10,29122);


%{
Concatenate the sobol set to one array
%}


X0 = [X01, X02, X03, X04, X05, X06, X07, X08, X09, X10];

save GcamSobolGen.mat X0 -v7.3;
clear
%}

load('GcamSobolGen.mat', 'X0');
SensOrder = 1/3;
PreAllocatedMean = ones(29122,10000);
load('GcamMeanValues.mat','unnamed'); GcamMeanValues=unnamed;
%   UniformGcamMeanValues=bsxfun(@times,GcamMeanValues,X0);
UniformGcamMeanValues=bsxfun(@times,GcamMeanValues,PreAllocatedMean);
%% Correct for negative numbers
% UniformGcamMeanValNeg=(-1).*(UniformGcamMeanValues<0);
% SignGCAMValue = ((-1).*(UniformGcamMeanValues<0))+PreAllocatedMean;
% SignGCAMValue = sign(UniformGcamMeanValues);
%{
absUniformGcamMeanValues = abs(UniformGcamMeanValues);

%% MISTAKE Edit for negative numbers low UniformGcamMeanValues - SensOrder.*|UniformGcamMeanValues|
ALower = (1-SensOrder).*absUniformGcamMeanValues;
BUpper = (1+SensOrder).*absUniformGcamMeanValues;
CMode  = absUniformGcamMeanValues;
%}
ALower = (1-SensOrder).*(abs(UniformGcamMeanValues));
BUpper = (1+SensOrder).*(abs(UniformGcamMeanValues));
CMode  = (abs(UniformGcamMeanValues));
%{
CModeCDF = (CMode - ALower)./(BUpper - ALower);
LogicHiThanMode    = (X0>CModeCDF);
%}
%CModeCDF = (CMode - ALower)./(BUpper - ALower);
LogicHiThanMode    = (X0>((CMode - ALower)./(BUpper - ALower)));

%{
LogicLowOrThanMode = ~LogicHiThanMode;
TriangDistLowOrThanMode = LogicLowOrThanMode.*(ALower + ((X0.*((BUpper - ALower).*(CMode - ALower))).^(1/2)));
%}

TriangDistLowOrThanMode = (~LogicHiThanMode).*(ALower + ((X0.*((BUpper - ALower).*(CMode - ALower))).^(1/2)));
TriangDistHiThanMode    = LogicHiThanMode   .*(BUpper - (((1-X0).*((BUpper - ALower).*(BUpper - CMode))).^(1/2)));
absTriangDist = TriangDistLowOrThanMode + TriangDistHiThanMode;
save AbsTriangData.mat absTriangDist -v7.3;
% TriangDist = absTriangDist.*SignGCAMValue;
TriangDist = absTriangDist.*(sign(UniformGcamMeanValues));
save GcamTrianGen.mat TriangDist -v7.3;
%{
double triangular(double a,double b,double c) {
   double U = rand() / (double) RAND_MAX;
   double F = (c - a) / (b - a);
   if (U <= F)
      return a + sqrt(U * (b - a) * (c - a));
   else
      return b - sqrt((1 - U) * (b - a) * (b - c));
}

%}
%{
% p = 
%     Sobol point set in 3 dimensions (8.918019e+013 points)
%     Properties:
%               Skip : 1000
%               Leap : 100
%     ScrambleMethod : none
%         PointOrder : standard
% Use scramble to apply a random linear scramble combined with a random digital shift:

%% p = scramble(p,'MatousekAffineOwen');
% p = 
%     Sobol point set in 3 dimensions (8.918019e+013 points)
%     Properties:
%               Skip : 1000
%               Leap : 100
%     ScrambleMethod : MatousekAffineOwen
%         PointOrder : standard
% Use net to generate the first four points:

%% X0 = net(p,29122)
% X0 =
%     0.7601    0.5919    0.9529
%     0.1795    0.0856    0.0491
%     0.5488    0.0785    0.8483
%     0.3882    0.8771    0.8755
% Use parenthesis indexing to generate every third point, up to the 11th point:

%% X = p(1:3:11,:)
% X =
%     0.7601    0.5919    0.9529
%     0.3882    0.8771    0.8755
%     0.6905    0.4951    0.8464
%     0.1955    0.5679    0.3192

%}