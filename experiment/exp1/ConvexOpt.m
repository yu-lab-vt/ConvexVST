clc;clear;close all;
tic;
load dis1_3.mat
A = pdf;
var_ratio = ones(1,254);
numIter = 10;
A = A(2:255,:);
[numFun, numBin] = size(A);
numVar = numBin - 1;
result_c = zeros(1,numIter+1);
result_error = zeros(1,numIter+1);
% A = noise_distribution;

% % Convex Optimization
x0 = ones(1,numVar);
f = @(x)objectiveFunction(x,A,numFun,var_ratio);
f2 = @(x)validateFunction(x,pdf,numFun)
options = optimoptions('fminimax','Display','iter','MaxFunctionEvaluations',numVar*1000,'UseParallel',true);
% tic;
[x,fx,fx_max,exitflag,output] = fminimax(f,x0,[],[],ones(1,numVar),numVar,zeros(1,numVar),Inf*ones(1,numVar),[],options);
% [x,fx,fx_max,exitflag,output] = fminimax(f,x0,[],[],ones(1,numVar),numVar,0.2*ones(1,numVar),Inf*ones(1,numVar),[],options);
% toc
figure(1);plot(cumsum(x));hold on
figure(2);plot(f2(x));hold on
result_c(1) = median(f2(x));
result_error(1) = sum(abs(f2(x)-result_c(1)));


for ii = 1:numIter
% % calculate the deriative
x0 = x;
deri = zeros(numFun,numVar);
step = 1e-10;
for jj = 1:numVar
    temp = x;
    temp(jj) = temp(jj) + step;
    deri(:,jj) = f(temp) - f(x);
end
deri = deri/step;
f = @(x)gapFunction(x,x0,A,deri,fx,fx_max,numFun,var_ratio);
options = optimoptions('fminimax','Display','iter','MaxFunctionEvaluations',numVar*1000,'UseParallel',true);
% tic;
[x,fval,maxfval,exitflag,output] = fminimax(f,x0,[],[],ones(1,numVar),numVar,zeros(1,numVar),Inf*ones(1,numVar),[],options);
% toc
% figure(1);plot(cumsum(x));hold on;     
f = @(x)objectiveFunction(x,A,numFun,var_ratio);
fx = f(x);
figure(2);plot(f2(x));hold on
result_c(ii+1) = median(f2(x));
result_error(ii+1) = sum(abs(f2(x)-result_c(ii+1)));
end
variance = f2(x);
save('result1_3.mat','result_c','result_error','variance','x');

toc

function f = gapFunction(x,x0,A,deri,fx,fx_max,numFun,var_ratio)
    for ii = 1:numFun
        f(ii) = fx_max - fx(ii) - sum((x-x0).*deri(ii,:));
        f(ii+numFun) = sum(A(ii,:).*([0 cumsum(x)]-sum(x(1:ii))).^2) - fx_max;
    end
end

function F = objectiveFunction(x,A,numFun,var_ratio)
    signal = [0 cumsum(x)];
    for ii = 1:numFun
        F(ii) = sum(A(ii,:).*(signal- signal(ii+1)).^2);
    end

end

function F2 = validateFunction(x,pdf,numFun)
    signal = [0 cumsum(x)];
    for ii = 1:numFun+2
        F2(ii) = sum(pdf(ii,:).*(signal- signal(ii)).^2);
    end

end
