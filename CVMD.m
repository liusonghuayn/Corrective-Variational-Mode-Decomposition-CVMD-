function [Kpro, IMF_best, omega_best] = CVMD(xn, Fs, kmax)
% CVMD  Kurtosis–based adaptive selection of the number of VMD modes.
%
% Usage:
%   [Kpro, IMF_best, omega_best] = CVMD(xn, Fs, kmax);
%
% Inputs:
%   xn   : input signal (1×N or N×1)
%   Fs   : sampling frequency
%   kmax : maximum number of candidate modes (optional, default 10)
%
% Outputs:
%   Kpro      : optimal number of modes (K)
%   IMF_best  : IMFs obtained by VMD with Kpro (Kpro×N)
%   omega_best: center frequencies associated with IMF_best

    if nargin < 3
        kmax = 10;           % default maximum K
    end

    % Ensure row vector
    xn   = xn(:).';
    len  = length(xn);
    aim  = round(log2(len)); %#ok<NASGU>

    % If amplitude normalization is needed, enable the following line:
    % xnn = xn/std(xn);
    xnn  = xn;

    % VMD parameters
    alpha = 1000;    % moderate bandwidth constraint
    tau   = 0;       % noise-tolerance (no strict fidelity enforcement)
    DC    = 0;
    init  = 1;       % initialize omegas uniformly
    tol1  = 1e-7;

    % Containers for evaluation values
    KU  = zeros(1, kmax-1);
    PCC = zeros(1, kmax-1);

    % “Current best” variables to avoid running VMD once more at the end
    bestKU      = inf;
    bestK       = 2;
    IMF_best    = [];
    omega_best  = [];

    for ii = 2:kmax
        K = ii;

        % -------- VMD decomposition --------
        [u, ~, omega] = VMD(xnn, alpha, tau, K, DC, init, tol1);  % u: K×N

        % Correlation between reconstructed signal and original signal
        % (overall reconstruction quality, optional)
        PCC(ii-1) = corr(xn', sum(u,1)');

        % -------- Power analysis & spectrum‑weighted kurtosis --------
        [~, ~, realy] = localPowerAnalysis(u, Fs);

        % Compute signal power and set threshold
        signal_power = sum(xn.^2) / len;
        Fthreshold   = signal_power / 100;

        newrealy = realy;
        newrealy(newrealy < Fthreshold) = 0;

        modesize = size(u);
        ku2      = zeros(modesize(1),1);

        for i = 1:modesize(1)
            miszero = (newrealy(i,:) ~= 0);
            ku      = kurtosis(u(i,:));
            norm_0  = sum(miszero);
            ku2(i)  = norm_0 * ku;
        end

        KU(ii-1) = sum(ku2);

        % -------- Update “current best” solution if this K is better --------
        if KU(ii-1) < bestKU
            bestKU     = KU(ii-1);
            bestK      = K;
            IMF_best   = u;        % save IMFs for this K
            omega_best = omega;    % save corresponding center frequencies
        end
    end

    % Final optimal K
    Kpro = bestK;
end


%% ================== Local power analysis function ==================
function [FY, realf, realy] = localPowerAnalysis(xn, Fs)
% Local power spectrum analysis.
% Intended for internal use by CVMD only.

    [a,~] = size(xn);

    for i = 1:a
        y    = xn(i,:);
        nfft = length(y);

        % --- Force all lengths/orders to be integers to avoid hanning warnings ---
        winLen   = max(1, round(nfft/4));  % window length
        window   = hanning(winLen);
        noverlap = floor(winLen/2);        % 50% overlap
        nfft_psd = 2*nfft;

        [P11, ~] = pwelch(y, window, noverlap, nfft_psd, Fs);
        realy(i,:)   = P11;
    end

    % Power spectrum peak of the original (summed) signal
    if a==1
        y1 = xn;
    else
        y1 = sum(xn);
    end

    nfft    = length(y1);
    winLen  = max(1, round(nfft/4));
    window  = hanning(winLen);
    noverlap = floor(winLen/2);
    nfft_psd = 2*nfft;

    [P11, realf] = pwelch(y1, window, noverlap, nfft_psd, Fs);
    realy1 = P11;
    FY     = max(realy1);
end