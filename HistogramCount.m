function [his,variance,parameters] = HistogramCount(data,kerSize,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to generate histogram from the input data.
%
%
% Input: data ----------- 2D/3D/4D imaging data.
%        kerSize -------- A one-row array with 2/3 elements represents 
%                         the size of a 2D/3D median kernel used to 
%                         calculate the real signal of one pixel by 
%                         calculating the median of its neighbors. 
%                         Elements must be odd.
%        
%        options: A structure contains some choices.
%        .histEdges ----- The edge of histograms. Pixels whose real 
%                         signals are between two adjacent edges will 
%                         be plotted in one histogram. Default value 
%                         is 0.5:254.5.
%        .binEdges ------ The edge of bins in each histograms. Default 
%                         value is -0.5:255.5.
%        .sampleSize ---- Default = 200. The smallest amount of samples 
%                         of each histograms. The histogram whose 
%                         number of samples is smaller than this 
%                         threshold will be merged to the next histogram.
%        .ratio --------- Default = 0.03. The proportion cutted in the
%                         truncated Gaussian fitting part.
%        .display ------- Default = True. Disaly the variacne cureve or
%                         not.
% Output: his ----------- An estimated histogram matrix. Each row is a
%                         histogram whose real signals are close to the
%                         corresponding entry in histCenters, and the
%                         boundary of real signals are the n and n+1
%                         entries in histEdges.
%         variance ------ Estimated variance of each histogram.
%         var_ratio ----- After truncated Gaussian fitting, the
%                         variance will be decreased. var_ratio is an
%                         estimated ratio of the variance of the
%                         truncated histogram to that of the real
%                         histogram.
%         histCenters --- New centers of histograms. Some histograms
%                         without enough samples will be merged. The
%                         new centers is the weighted mean of old
%                         centers.
%         histEdges ----- New edges of histograms. Some histograms
%                         without enough samples will be merged.
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

    if nargin == 2
        options.histEdges = 0.5:254.5;
        options.binEdges = -0.5:255.5;
        options.sampleSize = 200;
        options.ratio = 0.03;
        options.display = true;
    elseif nargin == 3
        if ~isfield(options,'histEdges')
            options.histEdges = 0.5:254.5;
        end
        if ~isfield(options,'binEdges')
            options.binEdges = -0.5:255.5;
        end
        if ~isfield(options,'ratio')
            options.ratio = 0.03;
        end
        if ~isfield(options,'sampleSize')
            options.sampleSize = 200;
        end
        if ~isfield(options,'display')
            options.display = true;
        end
    else
        error('Wrong number of inputs.');
    end

    % start counting
    tic;
    disp('Start Counting...');
    data = double(data);
    kerDim = length(kerSize);
    if kerDim == 2
        kerSize = [kerSize 1];
    end
    dataSize = size(data);
    dataDim = length(dataSize);
    if dataDim <= 3
        realSignal = medfilt3(data,kerSize);
    else
        numFrame = dataSize(4);
        realSignal = zeros(dataSize);
        for ii = 1:numFrame
            realSignal(:,:,:,ii) = medfilt3(data(:,:,:,ii),kerSize);     % most time consuming part
        end
    end

    kerRadi = floor(kerSize/2);
    if dataDim == 2
        data = data(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2));
        realSignal = realSignal(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2));
    elseif dataDim == 3
        data = data(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2),...
            kerRadi(3)+1:end-kerRadi(3));
        realSignal = realSignal(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2),...
            kerRadi(3)+1:end-kerRadi(3));
    elseif dataDim == 4
        data = data(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2),...
            kerRadi(3)+1:end-kerRadi(3),:);
        realSignal = realSignal(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2),...
            kerRadi(3)+1:end-kerRadi(3),:);
    end
    data = data(:);                  
    realSignal = realSignal(:);

    histEdges = options.histEdges;
    binEdges = options.binEdges;
    cluster = discretize(realSignal,histEdges);
    index = label2idx(cluster); 

    numOfHists = length(histEdges) - 1;
    numOfBins = length(binEdges) - 1;
    histCenters = (histEdges(1:end-1) + histEdges(2:end))/2;
    binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;
%     histCenters = options.histCenters;
%     binCenters = options.binCenters;

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
        if length(index{ii}) < options.sampleSize
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
    var_ratio = zeros(1,length(index));
    ratio = options.ratio;
    his = zeros(length(index),numOfBins);
    for ii = 1:length(index)
        noise = data(index{ii}) - realSignal(index{ii});
        distribution = histogram(noise + histCenters(ii), binEdges);
        his(ii,:) = distribution.Values/length(noise);

        % trauncated Gaussian fitting
        left_ind = 1; right_ind = numOfBins;
        ratio_residual = ratio;
        ratio_left = 0;
        ratio_right = 0;
        while ratio_residual > 0 && (left_ind < ii || right_ind > ii)
            while his(ii,left_ind) == 0 && left_ind < ii
                left_ind = left_ind + 1;
            end
            while his(ii,right_ind) == 0 && right_ind > ii
                right_ind = right_ind - 1;
            end
            if binCenters(ii) - binCenters(left_ind) <= binCenters(right_ind) - binCenters(ii)
                % truncate at right part
                if his(ii,right_ind) <= ratio_residual
                    ratio_residual = ratio_residual - his(ii,right_ind);
                    ratio_right = ratio_right + his(ii,right_ind);
                    his(ii,right_ind) = 0;
                else
                    ratio_right = ratio_right + ratio_residual;
                    his(ii,right_ind) = his(ii,right_ind) - ratio_residual;
                    ratio_residual = 0;
                end
            else
                % truncate at left part
                if his(ii,left_ind) <= ratio_residual
                    ratio_residual = ratio_residual - his(ii,left_ind);
                    ratio_left = ratio_left + his(ii,left_ind);
                    his(ii,left_ind) = 0;
                else
                    ratio_left = ratio_left + ratio_residual;
                    his(ii,left_ind) = his(ii,left_ind) - ratio_residual;
                    ratio_residual = 0;
                end
            end
        end
        if(ratio_left ~= 0 || ratio_right~=0)
            if(ratio_left == 0)
                beta = norminv(1 - ratio_right);
                phi_beta = normpdf(beta);
                var_ratio(ii) = 1 + (- beta*phi_beta)/(1 - ratio_right) ...
                - ((phi_beta)/(1 - ratio_right))^2;
            elseif(ratio_right == 0)
                alpha = norminv(ratio_left);
                phi_alpha =  normpdf(alpha);
                var_ratio(ii) = 1 + (alpha*phi_alpha)/(1 - ratio_left) ...
                - ((phi_alpha)/(1 - ratio_left))^2;
            else
                alpha = norminv(ratio_left);
                beta = norminv(1 - ratio_right);
                phi_alpha = normpdf(alpha);
                phi_beta = normpdf(beta);
                var_ratio(ii) = 1 + (alpha*phi_alpha - beta*phi_beta)/(1 - ratio_left - ratio_right) ...
                - ((phi_alpha - phi_beta)/(1 - ratio_left - ratio_right))^2;
            end
        else
            var_ratio(ii) = 1;
        end

        his(ii,:) = his(ii,:)/sum(his(ii,:));
        variance(ii) = sum(his(ii,:).*(binCenters-histCenters(ii)).^2)/var_ratio(ii);
    end
    disp(['Finished. Counting part running time is ' num2str(toc) ' seconds.']);
    parameters.var_ratio = var_ratio;
    parameters.histCenters = histCenters;
    parameters.histEdges = histEdges;
    parameters.binCenters = binCenters;
    parameters.binEdges = binEdges;
    if options.display == true
        plot(histCenters,variance,'o');
    end
