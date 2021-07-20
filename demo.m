clc;clear;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ConvexVST is a variance -stabilization algorithm for discrete random
% variables based on convex optimization. This code is a demo to show
% how ConvexVST works.
% Please gp to http://proceedings.mlr.press/v139/wang21p.html for more info

% ----- Input  -----
% pdf:         The distribution of the family of random variables. (M x N)
%              M is the number of choices in the sample space and N is the
%              number of bins. 

% ----- Output -----
% stabFunction:Transform function. (1 x M)
% variance:    Stabilized variance after transform. (1 x M-2)(without ends)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
% Everyone is permitted to copy and distribute verbatim copies
% of this license document, but changing it is not allowed.
%
% AUTHORS:
%     Mengfan Wang
%     email: mengfanw@vt.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ./pdf.mat
plotFlag = 1;                 % 1 for plotting transform function and variance curve
numIter = 10;                 % iterations of the second part
A = pdf;                            
A = A(2:255,:);               % Remove end points to keep the best performance.
[numFun, numBin] = size(A);
numVar = numBin - 1;

% Phase 1: minimize maximum variance
x0 = ones(1,numVar);
f = @(x)objectiveFunction(x,A,numFun);                 % Objective function
options = optimoptions('fminimax','Display','iter','MaxFunctionEvaluations',numVar*1000,'UseParallel',true);
[x,fx,fx_max,exitflag,output] = fminimax(f,x0,[],[],ones(1,numVar),numVar,zeros(1,numVar),Inf*ones(1,numVar),[],options);

stabFunction = [0 cumsum(x)];
variance = f(x);
if plotFlag
    figure(1);plot(stabFunction);hold on;    
    figure(2);plot(variance);hold on
end

% Part 2: minimize absolute error iteratively
for ii = 1:numIter
% calculate the deriative
x0 = x;
deri = zeros(numFun,numVar);
step = 1e-10;
for jj = 1:numVar
    temp = x;
    temp(jj) = temp(jj) + step;
    deri(:,jj) = f(temp) - f(x);
end
deri = deri/step;
f = @(x)gapFunction(x,x0,A,deri,fx,fx_max,numFun);     % Objective function
options = optimoptions('fminimax','Display','iter','MaxFunctionEvaluations',numVar*1000,'UseParallel',true);
[x,fval,maxfval,exitflag,output] = fminimax(f,x0,[],[],ones(1,numVar),numVar,zeros(1,numVar),Inf*ones(1,numVar),[],options);
f = @(x)objectiveFunction(x,A,numFun);
fx = f(x);
end
stabFunction = [0 cumsum(x)];
variance = f(x);

if plotFlag
    figure(1);plot(stabFunction);hold on;    
    legend('Part 1 result', 'Part 2 result');
    figure(2);plot(variance);hold on
    legend('Part 1 variance', 'Part 2 variance');
end


function F = objectiveFunction(x,A,numFun)
    signal = [0 cumsum(x)];
    for ii = 1:numFun
        F(ii) = sum(A(ii,:).*(signal- signal(ii+1)).^2);
    end
end


function f = gapFunction(x,x0,A,deri,fx,fx_max,numFun)
    for ii = 1:numFun
        f(ii) = fx_max - fx(ii) - sum((x-x0).*deri(ii,:));
        f(ii+numFun) = sum(A(ii,:).*([0 cumsum(x)]-sum(x(1:ii))).^2) - fx_max;
    end
end


