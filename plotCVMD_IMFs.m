function plotCVMD_IMFs(t, x, IMFs, plotTitle)
% plotCVMD_IMFs  Plot original signal and CVMD IMFs
% with identical y-limits for all IMF subplots.
%
% Inputs:
%   t         : time axis (1×N or N×1)
%   x         : original signal (1×N or N×1)
%   IMFs      : IMFs from CVMD, size K×N
%   plotTitle : (optional) title string for sgtitle
%
% Example:
%   plotCVMD_IMFs(t, xn, IMF_best, 'CVMD');

    if nargin < 4
        plotTitle = 'CVMD decomposition';
    end

    % Ensure shapes
    t = t(:).';
    x = x(:).';
    if size(IMFs,2) ~= numel(t)
        error('IMFs must have size K×N with N == numel(t).');
    end

    [K, N] = size(IMFs); %#ok<ASGLU>

    % ----- Compute common y-limits for all IMFs -----
    allIMFvalues = IMFs(:);
    ymin = min(allIMFvalues);
    ymax = max(allIMFvalues);
    if ymin == ymax
        ymin = ymin - 1; ymax = ymax + 1;
    end
    yrange = ymax - ymin;
    ylimCommon = [ymin - 0.05*yrange, ymax + 0.05*yrange];

    % ----- Create figure and subplots -----
    figure('Color','w');
    ax = gobjects(K+1,1);

    % (1) Original signal
    ax(1) = subplot(K+1,1,1);
    plot(t, x, 'k', 'LineWidth',1);
    ylabel('Signal');
    set(ax(1), 'FontName','Times New Roman','FontSize',10, ...
        'XLim',[t(1) t(end)], 'Box','on');
    title('Original signal');

    % (2…K+1) Each IMF
    for k = 1:K
        ax(k+1) = subplot(K+1,1,k+1);
        plot(t, IMFs(k,:), 'b', 'LineWidth',1);
        ylabel(sprintf('IMF_%d',k));
        set(ax(k+1), 'FontName','Times New Roman','FontSize',10, ...
            'XLim',[t(1) t(end)], ...
            'YLim', ylimCommon, ...        % same y-range for all IMFs
            'Box','on');
        if k == K
            xlabel('Time (s)');
        end
    end

    % Link x-axes for easier zoom/pan
    linkaxes(ax,'x');

    % Global title
    sgtitle(plotTitle, 'FontName','Times New Roman','FontSize',12);
end