function anaPlot(ana,savePath)
% ANAPLOT produce several plots corresponds to the evaluation information contained
% in model analysis structure.
%
%   anaPlot(ana) would show analytic plots according to information in structure ANA.
%
%   anaPlot(ana,savePath) would produce analytic plots and save them into graphic
%   files under the folder specified by SAVEPATH.
%
% see also modelAnalysis, anaParamSetup.
%
% MooGu Z. <hzhu@case.edu>
% -----------------------------------------------------------------------------------
% Version 0.1 [June 22, 2014] - Start Version

swPrint = exist('savePath','var');

% when print plots into files only one invisible frame would be created
if swPrint
    if ~exist(savePath,'dir')
        system(['mkdir -p "',savePath,'"']);
    end

    f = figure('visible','off');
end

% draw SNR plot
if ~swPrint, f = figure(); end
plot(ana.summary.SNR,'o-k','MarkerSize',7); hold on; grid on;
plot(ana.summary.SNR,'.r','MarkerSize',13); hold off;
title('SNR Curve');
xlabel('# of iterations / 10,000');
ylabel('SNR / dB');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'SNRCurve.eps']);
end

% draw Probability Plots
if ~swPrint, f = figure(); end
plot(ana.summary.prob,'o-k','MarkerSize',7); hold on; grid on;
curveProb = plot(ana.summary.prob,'.g','MarkerSize',13);
plot(ana.summary.likelihood,'o-k','MarkerSize',7);
curveLh = plot(ana.summary.likelihood,'.r','MarkerSize',13);
plot(ana.summary.sparse,'o-k','MarkerSize',7);
curveSparse = plot(ana.summary.sparse,'.y','MarkerSize',13);
plot(ana.summary.slow,'o-k','MarkerSize',7);
curveSlow = plot(ana.summary.slow,'.b','MarkerSize',13); hold off;
legend([curveProb,curveLh,curveSparse,curveSlow],'Joint Probability', ...
       'Likelihood(Noise)','Sparseness','Slowness');
title('Probability Curves');
xlabel('# of iterations / 10,000');
ylabel('Probabilities');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'ProbCurve.eps']);
end

% independent likelihood plot 
if ~swPrint, f = figure(); end
plot(ana.summary.likelihood,'o-k','MarkerSize',7); hold on; grid on;
plot(ana.summary.likelihood,'.r','MarkerSize',13); hold off;
title('Likelihood (Noise) Curves');
xlabel('# of iterations / 10,000');
ylabel('Probabilities');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'LikelihoodCurve.eps']);
end

% independent sparseness plot
if ~swPrint, f = figure(); end
plot(ana.summary.sparse,'o-k','MarkerSize',7); hold on; grid on;
plot(ana.summary.sparse,'.y','MarkerSize',13); hold off;
title('Sparseness Curves');
xlabel('# of iterations / 10,000');
ylabel('Probabilities');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'SparseCurve.eps']);
end

% independent slowness plot 
if ~swPrint, f = figure(); end
plot(ana.summary.slow,'o-k','MarkerSize',7); hold on; grid on;
plot(ana.summary.slow,'.b','MarkerSize',13); hold off;
title('Slowness Curves');
xlabel('# of iterations / 10,000');
ylabel('Probabilities');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'SlowCurve.eps']);
end

[nmodel,nbox] = size(ana.resp.dist.damp);
% Decide model profile that would be plot
if nmodel == 1
    index = 1;
elseif nmodel == 2
    index = [1,2];
else
    index = [1,ceil(nmodel/2),nmodel];
end
% get iteration number of each model profile
iterNum  = cell(size(index));
distAmp  = zeros(numel(index),nbox);
distDAmp = zeros(numel(index),nbox);
for i = 1 : numel(index)
    iterNum{i} = sum(cellfun(@str2num,regexp(ana.model.nameList{index(i)}, ...
        '-?\d+\.?\d*|-?\d*\.?\d+','match')));
    distAmp(i,:) = ana.resp.dist.amp(index(i),:) / ...
        sum(ana.resp.dist.amp(index(i),:));
    distDAmp(i,:) = ana.resp.dist.damp(index(i),:) / ...
        sum(ana.resp.dist.damp(index(i),:));

end
% define a temporal function to compose legends string
lgdstr = @(num) (['iteration - ',num2str(num)]);
% distribution plots
cmap = cmapgen(numel(index));
axHandle = zeros(1,numel(index));
% distribution of 'a'
if ~swPrint, f = figure(); else clf; end
hold on
for i = 1 : numel(index)
    plot(ana.resp.dist.ampcrd,distAmp(i,:),'o-k','MarkerSize',7);
    axHandle(i) = plot(ana.resp.dist.ampcrd,distAmp(i,:),'.', ...
        'MarkerSize',13,'Color',cmap(i,:));
end
legend(axHandle,cellfun(lgdstr,iterNum,'UniformOutput',false)); 
grid on; hold off;
title('Distribution of a');
xlabel('Value of a');
ylabel('Frequency');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'aDist.eps']);
end
% distribution of 'delta a'
if ~swPrint, f = figure(); else clf; end
hold on
for i = 1 : numel(index)
    plot(ana.resp.dist.dampcrd,distDAmp(i,:),'o-k','MarkerSize',7);
    axHandle(i) = plot(ana.resp.dist.dampcrd,distDAmp(i,:),'.', ...
        'MarkerSize',13,'Color',cmap(i,:));
end
legend(axHandle,cellfun(lgdstr,iterNum,'UniformOutput',false)); 
grid on; hold off;
title('Distribution of \Delta a');
xlabel('Value of \Delta a');
ylabel('Frequency');
if swPrint
    print(f,'-depsc2','-r300',[savePath,'daDist.eps']);
end

end

