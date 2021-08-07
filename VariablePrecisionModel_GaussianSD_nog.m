% VARIABLEPRECISIONMODEL_GAUSSIANSD returns a structure for a variable precision mixture model
% in which the standard deviations of observers' reports are assumed to be
% themselves distributed as a normal distribution.

function model = VariablePrecisionModel_GaussianSD_nog()
  model.name = 'Variable precision model (gaussian over sd)';
	model.paramNames = {'mnSTD', 'stdSTD'};
	model.lowerbound = [0 0]; % Lower bounds for the parameters
	model.upperbound = [100 100]; % Upper bounds for the parameters
	model.movestd = [0.15, 0.07];
	model.pdf = @vp_pdf;
  model.modelPlot = @model_plot;
	model.start = [15, 5;
                   20, 10;
                   10, 2;
                   30, 3];

  % To specify a prior probability distribution, change and uncomment
  % the following line, where p is a vector of parameter values, arranged
  % in the same order that they appear in model.paramNames:
  % model.prior = @(p) (1);

  % For speed, calculate these all out here
  stdsSumOver = linspace(0.5, 100, 500);
  kValues = deg2k(stdsSumOver)';
  baseK = log(besseli(0, kValues, 1)) + kValues;
  lastX = [];

  function y = vp_pdf(data,mnSTD,stdSTD)
      g=0;
    % Probability of each of these
    probEachSD = normpdf(stdsSumOver, mnSTD, stdSTD);
    probEachSD = probEachSD./sum(probEachSD);

    if length(data.errors)~=length(lastX) || any(data.errors~=lastX)
      % Calculate pdf for each STD; only do if the data is different than
      % last time
      x = repmat(data.errors(:), [1 length(kValues)]);
      k = repmat(kValues', [length(data.errors) 1]);
      newBaseK = repmat(baseK', [length(data.errors) 1]);
      model.v = exp((k.*cos((pi/180)*x)) - (log(360) + newBaseK));
      lastX = data.errors;
    end

    % Make final model prediction and sum
%     probDataUnderThisNormal = (1-g).*model.v + (g).*1/360;
    probDataUnderThisNormal = model.v;
    probEachSDBig = repmat(probEachSD, [size(probDataUnderThisNormal,1), 1]);
    y = sum(probDataUnderThisNormal.*probEachSDBig,2);
  end

  % Use our custom modelPlot to make a higher-order distribution plot
  function figHand = model_plot(data, params, varargin)
    figHand = figure();
    if isstruct(params) && isfield(params, 'vals')
      maxParams = MCMCSummarize(params, 'maxPosterior');
      which = randsample(size(params.vals,1), 100);
      likeVals = params.like(which);
      params = params.vals(which,:);
    else
      maxParams = params;
      likeVals = 1;
    end
    set(gcf, 'Color', [1 1 1]);
    x = stdsSumOver;
    for i=1:size(params,1)
      y = normpdf(stdsSumOver, params(i,1), params(i,2));
      colorOfLine = fade([0.54, 0.61, 0.06], ...
        exp(likeVals(i) - max(likeVals)));
      plot(x, y, 'Color', colorOfLine); hold on;
    end
    y = normpdf(stdsSumOver, maxParams(1), maxParams(2));
    plot(x, y, 'Color', [0.54, 0.61, 0.06], 'LineWidth', 3);
    title('Higher-order distribution', 'FontSize', 14);
    xlabel('Standard dev. (degrees)', 'FontSize', 14);
    ylabel('Probability', 'FontSize', 14);
  end
end
