function [variance, stabilizeFunction] = ConvexVST(his,parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ConvexVST is a variance -stabilization algorithm for discrete random
% variables based on convex optimization. This code is a demo to show
% how ConvexVST works.
% Please gp to http://proceedings.mlr.press/v139/wang21p.html for more info

% ----- Input  -----
% hist:         The distribution of the family of random variables. (M x N)
%               M is the number of choices in the sample space and N is the
%               number of bins. 
% parameters.
%   histCenters:The choices in the sample space in ascending order.
%               (required = True, 1 x M) 
%   binCenters: Values of bins in ascending order.(required = True, 1 x N)
%   varRatio:   The weight of different histograms. Smaller varRatio means 
%               larger weight. (default = all ones, 1 x M)
%   algorithm:  Solver name you want to use. 'matlab' will use the matlab 
%               solver fminmax, 'mosek' will use the commerical Mosek
%               solver. It's faster. (default = 'matlab')
%   mosekPath:  The path of Mosek solver if you want o use.
%               (required = True if algorithm == 'mosek')
%   
%
% ----- Output -----
% variance:         Stabilized variance after transform. (1 x M)
% stabilizeFunction:Transform function. (1 x N)
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

if any(sum(his,2) < 1 - 1e-6) || any(sum(his,2) > 1 + 1e-6)
    warning('The sum of some histograms may be not 1.');
    his = his./sum(his,2);
end
if ~isfield(parameters,'histCenters')
error('parameters.histCenters is necessary.');
end
if ~isfield(parameters,'binCenters')
error('parameters.binCenters is necessary.');
end
histCenters = parameters.histCenters;
binCenters = parameters.binCenters; 
if min(histCenters) < min(binCenters) || max(histCenters) > max(binCenters)
    error('The range of histCenters should not exceed the range of binCenters.');
end
numOfHists = length(histCenters);
numOfBins = length(binCenters);
if ~isfield(parameters,'varRatio')
    varRatio = ones(1,numOfHists);
else
    varRatio = parameters.varRatio;
end
if ~isfield(parameters,'algorithm')
    algorithm = 'matlab';
else
    algorithm = 'mosek';
end

% optimization part first step
tic;
disp('Start First Step...');

if strcmp(algorithm, 'mosek')
    % MOSEK solver
    if ~isfield(parameters,'mosekPath')
        error('Please give the correct path of MOSEK solver.');
    end
    addpath(parameters.mosekPath);
    Hess = cell(1,numOfHists);
    for ii = 1:numOfHists
        Hess{ii} = zeros(numOfBins-1,numOfBins-1);
        left_ind = find(binCenters<histCenters(ii),1,'last'); % also the split variable
        right_ind = left_ind+1; 
        alpha = (histCenters(ii) - binCenters(left_ind))/(binCenters(right_ind) - binCenters(left_ind));
        beta = (binCenters(right_ind) - histCenters(ii))/(binCenters(right_ind) - binCenters(left_ind));
        if isempty(left_ind)
            left_ind = 0; right_ind = 1; alpha = 1; beta = 0;
        end
        for jj = 1:numOfBins
            if jj < left_ind
                Hess{ii}(jj:left_ind-1,jj:left_ind-1) = Hess{ii}(jj:left_ind-1,jj:left_ind-1) + his(ii,jj);
                Hess{ii}(left_ind,jj:left_ind-1) = Hess{ii}(left_ind,jj:left_ind-1) + alpha*his(ii,jj);
                Hess{ii}(jj:left_ind-1,left_ind) = Hess{ii}(jj:left_ind-1,left_ind) + alpha*his(ii,jj);
                Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + alpha^2*his(ii,jj);
            elseif jj == left_ind
                Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + alpha^2*his(ii,jj);
            elseif jj == right_ind
                if left_ind > 0
                    Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + beta^2*his(ii,jj);
                end
            elseif jj > right_ind
                Hess{ii}(right_ind:jj-1,right_ind:jj-1) = Hess{ii}(right_ind:jj-1,right_ind:jj-1) + his(ii,jj);
                if left_ind > 0 
                    Hess{ii}(left_ind,right_ind:jj-1) = Hess{ii}(left_ind,right_ind:jj-1) + beta*his(ii,jj);
                    Hess{ii}(right_ind:jj-1,left_ind) = Hess{ii}(right_ind:jj-1,left_ind) + beta*his(ii,jj);
                    Hess{ii}(left_ind,left_ind) = Hess{ii}(left_ind,left_ind) + beta^2*his(ii,jj);
                end
            end
        end
        Hess{ii} = tril(Hess{ii})*2/varRatio(ii);
        Hess{ii}(Hess{ii}<1e-16) = 0;
    end
    prob.qcsubk = [];
    prob.qcsubi = [];
    prob.qcsubj = [];
    prob.qcval = [];
    for ii = 1:numOfHists
        [row,column] = find(Hess{ii});
        prob.qcsubk = [prob.qcsubk; ii*ones(length(row),1)];
        prob.qcsubi = [prob.qcsubi; row];
        prob.qcsubj = [prob.qcsubj; column];
        prob.qcval = [prob.qcval; Hess{ii}(sub2ind(size(Hess{ii}),row,column))];
    end
    prob.c = [zeros(1,numOfBins-1) 1];
    prob.buc = [zeros(numOfHists,1); binCenters(end)-binCenters(1)];
    prob.blc = [-inf*ones(numOfHists,1); binCenters(end)-binCenters(1)];
    prob.blx = zeros(numOfBins,1);
    prob.a = repmat([zeros(1,numOfBins-1) -1],[numOfHists 1]);
    prob.a = [prob.a; [ones(1,numOfBins-1) 0]];
    param.MSK_IPAR_INTPNT_SOLVE_FORM = 'MSK_SOLVE_DUAL';
    [~,res] = mosekopt('minimize',prob,param);

    x = res.sol.itr.xx(1:end-1)';
    stabilizeFunction = cumsum([binCenters(1) x]);
    f = @(x)objectiveFunction(x,his,numOfHists,binCenters,histCenters);
    variance = f(x);
elseif strcmp(algorithm,'matlab')
    numVar = numOfBins - 1;
    x0 = diff(binCenters);
    f = @(x)objectiveFunction(x,his,numOfHists,varRatio);                 
    options = optimoptions('fminimax','Display','iter','MaxFunctionEvaluations',numVar*1000,'UseParallel',true);
    [x,fx,fx_max,exitflag,output] = fminimax(f,x0,[],[],ones(1,numVar),numVar,zeros(1,numVar),Inf*ones(1,numVar),[],options);
    stabilizeFunction = [0 cumsum(x)];
    variance = f(x);
else
    error('Please input the correct solver name. Support matlab and mosek only');
end
end

function F = objectiveFunction(x,A,numFun,varRatio)
    signal = [0 cumsum(x)];
    for ii = 1:numFun
        F(ii) = sum(A(ii,:).*(signal- signal(ii+1)).^2)/varRatio;
    end
end