clc;clear;close all;
% dbstop if error

addpath('C:\Program Files\Mosek\9.2\toolbox\R2015a');
result_error = zeros(4,3);
result_bar = zeros(4,3);
for dd = 1:3
    close all;
load(['simulation_' num2str(dd*5) '.mat'],'data_simulation','variance','realSignal');
data = data_simulation;
variance_original = variance;
realSignal_original = realSignal;

% our method
kerSize = [3 3 1];
numFrame = size(data,4);
histEdges = 0.5:254.5;
binEdges = -0.5:255.5;
options.histEdges = 0.5:254.5;
options.binEdges = -0.5:255.5;
options.histCenters = 1:254;
options.binCenters = 0:255;
options.sampleSize = 200;
options.ratio = 0.03;
parameters_original= options;

figure(2);
[his,variance,parameters] = HistogramCount(data,kerSize,options); 
his_original = his;
[variance, stabilizeFunction] = optimization(his,parameters);
figure(3);plot(0:255,stabilizeFunction,'LineWidth',2);hold on

realSignal_stab = realSignal_original+1;
realSignal_stab = stabilizeFunction(realSignal_stab);
data_stab = data+1;
data_stab = stabilizeFunction(data_stab);
figure(2);
[variance,his,parameters] = histogramCount(data_stab,realSignal_stab,binEdges,binEdges);  % validation
inv_stab = inverseFunction(stabilizeFunction);
figure(1);plot(inv_stab(parameters.histCenters+1),variance,'.','MarkerSize',10);hold on;
result_error(dd,1) = mean(abs(variance-median(variance)));
result_bar(dd,1) = max(variance) - min(variance);
result_error(dd,4) = mean(abs(variance(2:end-1)-median(variance(2:end-1))));
result_bar(dd,4) = max(variance(2:end-1)) - min(variance(2:end-1));

% generalized anscome transformation
start_ind = floor(length(variance_original)*0.1);
end_ind = ceil(length(variance_original)*0.5);
linear_fit = fit(parameters_original.histCenters(start_ind:end_ind)',variance_original(start_ind:end_ind)','poly1');
poisson = sqrt(linear_fit.p1); gaussian = sqrt(max(linear_fit.p2,0));
stabilizeFunction = 2/poisson*sqrt(poisson*parameters_original.binCenters+3*poisson^2/8+gaussian^2);
c = (max(stabilizeFunction) - min(stabilizeFunction));
stabilizeFunction = (stabilizeFunction - min(stabilizeFunction))/(max(stabilizeFunction) - min(stabilizeFunction))...
                        * (max(parameters_original.binCenters) - min(parameters_original.binCenters));
figure(3);plot(0:255,stabilizeFunction,'LineWidth',2);              
realSignal_stab = realSignal_original+1;
realSignal_stab = stabilizeFunction(realSignal_stab);
data_stab = data+1;
data_stab = stabilizeFunction(data_stab);
figure(2);
[variance,his,parameters] = histogramCount(data_stab,realSignal_stab,binEdges,binEdges);  % validation
inv_stab = inverseFunction(stabilizeFunction);
figure(1);plot(inv_stab(parameters.histCenters+1),variance,'.','MarkerSize',10);hold on;
result_error(dd,2) = mean(abs(variance-median(variance)));
result_bar(dd,2) = max(variance) - min(variance);


% foi
ou = 1.5; ol = 1.5; ru = 0.8; rl = 0.8;
ru1 = 0.2; rl1 = 0.2; ru2 = 0.5; rl2 = 0.5;

sigmaf = sqrt(variance_original);
fk = 0:255;
realSignal = zeros(size(data));
for ii = 1:numFrame
    realSignal(:,:,:,ii) = medfilt3(data(:,:,:,ii),kerSize);     % most time consuming part
end
figure(2);
[variance,his,parameters] = histogramCount(data,realSignal,binEdges,binEdges);   % generate variance curve

for kk = 1:500
ef = sigmaf - c;
ef_bar = max(-rl2,min(ru2,ef));
sigmaf_bar = max(c-rl2,min(c+ru2,sigmaf));
x_0lf = ef_bar<0;
x_0rg = ef_bar>0;
x_ru1lf = ef_bar<=ru1;
x_ru1rg = ef_bar>ru1;
x_rl1lf = ef_bar<=-rl1;
x_rl1rg = ef_bar>-rl1;
phi_ef = ru*x_0rg.*( (1-((ef-ru1)/ru1).^2).^(ou-1).*x_ru1lf + x_ru1rg )...
      +  rl*x_0lf.*( (1-((ef+rl1)/rl1).^2).^(ol-1).*x_rl1rg + x_rl1lf );
ifk = 1 - phi_ef.*ef_bar./sigmaf_bar;
increment = (ifk(1:end-1)+ifk(2:end))/2;
fk = [ 0 cumsum(increment.*diff(fk))];
fk = fk/max(fk)*255;
for ii = 1:256
    variance(ii) = sum(his(ii,:).*(fk- fk(ii)).^2);
end
sigmaf = sqrt(variance);
end

figure(3);plot(0:255,fk,'LineWidth',2);   
figure(3);plot(0:255,0:255,'--','LineWidth',2);   
realSignal_stab = realSignal+1;
realSignal_stab = stabilizeFunction(realSignal_stab);
data_stab = data+1;
data_stab = stabilizeFunction(data_stab);
figure(2);
[variance,his,parameters] = histogramCount(data_stab,realSignal_stab,binEdges,binEdges);  % validation
inv_stab = inverseFunction(stabilizeFunction);
figure(1);plot(inv_stab(parameters.histCenters+1),variance,'.','MarkerSize',10);hold on;
result_error(dd,3) = mean(abs(variance-median(variance)));
result_bar(dd,3) = max(variance) - min(variance);
figure(1);plot([0:255],variance_original,'.','MarkerSize',10);

figure(3);xlabel('\theta','FontWeight','bold','FontSize',15);
ylabel('$f(\theta)$','FontWeight','bold','Interpreter','latex','FontSize',15);
legend('ConvexOpt','GAT','Foi','Original','Location','northwest');

% export_fig(['..\figure\exp3\exp3_fig1' num2str(dd)], '-transparent', '-m2')
figure(1);xlabel('\theta','FontWeight','bold','FontSize',15);
ylabel('$\sigma^2_{\theta}$','FontWeight','bold','Interpreter','latex','FontSize',15);
legend('ConvexOpt','GAT','Foi','Original','Location','northwest');

% export_fig(['..\figure\exp2\exp2_fig2' num2str(dd)], '-transparent', '-m2')
end
result_error


function [variance, stabilizeFunction] = optimization(his,parameters)
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
if ~isfield(parameters,'var_ratio')
    var_ratio = ones(1,numOfHists);
else
    var_ratio = parameters.var_ratio;
end
if ~isfield(parameters,'step')
    step = 1;
else
    step = parameters.step;
end
if ~isfield(parameters,'algorithm')
    algorithm = 'aaaaa';
else
    algorithm = 'bbbbb';
end

% optimization part first step
tic;
disp('Start First Step...');

% MOSEK Part
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
    Hess{ii} = tril(Hess{ii})*2/var_ratio(ii);
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
end

function F = objectiveFunction(x,his,numOfHists,binCenters,histCenters)
    for ii = 1:numOfHists
        left_ind = find(binCenters<histCenters(ii),1,'last');
        if isempty(left_ind)
            F(ii) = sum(his(ii,:).*([0 cumsum(x)]).^2);
        else
            beta = (binCenters(left_ind+1) - histCenters(ii))/(binCenters(left_ind+1) - binCenters(left_ind));
            F(ii) = sum(his(ii,:).*([0 cumsum(x)]-(sum(x(1:left_ind))-beta*x(left_ind))).^2);
        end
    end
end

function [variance,his,parameters] = histogramCount(data,realSignal,histEdges,binEdges)
    data = data(:);                  
    realSignal = realSignal(:);
    
    cluster = discretize(realSignal,histEdges);
    index = label2idx(cluster); 

    numOfHists = length(histEdges) - 1;
    numOfBins = length(binEdges) - 1;
    histCenters = (histEdges(1:end-1) + histEdges(2:end))/2;
    binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;

    if min(cluster) > 1
        warning(['The smallest histogram does not contain any samples. Please'...
            ' double check options.histEdges.']);
    end
    if max(cluster) < length(histCenters)
        warning(['The largest histogram does not contain any samples. Please'...
            ' double check options.histEdges.']);
        histEdges = histEdges(1:length(index)+1);
        histCenters = histCenters(1:length(index));
    end

    % merge empty histograms
    index_label = zeros(1,length(index));
    for ii = 1:length(index)-1
        if length(index{ii}) < 200
            index_label(ii) = 1;
            if length(index{ii+1}) + length(index{ii}) > 0 
                histCenters(ii+1) = (histCenters(ii+1)*length(index{ii+1}) + ...
                    histCenters(ii)*length(index{ii}))/(length(index{ii+1})+length(index{ii}));
            end
            index{ii+1} = [index{ii+1};index{ii};];
        end
    end
    histEdges(logical([0 index_label])) = []; 
    histCenters(logical(index_label)) = [];
    index = index(logical(1-index_label));

    variance = zeros(1,length(index));
    his = zeros(length(index),numOfBins);
    for ii = 1:length(index)
        noise = data(index{ii}) - realSignal(index{ii});
        distribution = histogram(noise + histCenters(ii), binEdges);
        his(ii,:) = distribution.Values/length(noise);
        variance(ii) = sum(his(ii,:).*(binCenters-histCenters(ii)).^2);
    end
    parameters.histCenters = histCenters;
    parameters.histEdges = histEdges;
    parameters.binCenters = binCenters;
    parameters.binEdges = binEdges;
end


function inv_fx = inverseFunction(fx)

% For any monotonic increasing function f(x) in [0-255], this function will output
% f-1(x)

len = length(fx);
inv_fx = zeros(1,len);

for ii = 2:len
    if round(fx(ii)) >= ceil(fx(ii-1))
        for jj = ceil(fx(ii-1)):round(fx(ii))
            if jj > 0
                inv_fx(jj+1) = (jj - fx(ii-1))/(fx(ii) - fx(ii-1)) + ii - 2;
            end
        end
    end
end 
end