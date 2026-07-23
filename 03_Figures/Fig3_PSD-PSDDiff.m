%% =========================================================
% Fig4 grouped plot
% Left panel:
%   T01 / T02 / ZSD02 whole-period PSD
%   plotted as median + IQR (25%–75%)
%
% Right panel:
%   PSDdiff = PSD(T02) - PSD(T01)
%   computed from whole-period hourly PSD
%   plotted as median + IQR + zero line
%% =========================================================
clc;
clear;
close all;

% Paths
rd4PSD_MAIN = '/home/lyang/Disk2/DATA_HDB/Test_SUSTech/Result/PSD/PSD_WT';
rd4PSD_ZSD  = '/home/lyang/Disk2/DATA_HDB/Test_SUSTech/Result/PSD/PSD_WT/PSD_ZSD/PSD_WT/data';
rd4RESULT   = '/home/lyang/Disk2/DATA_HDB/Test_SUSTech/Result/Fig/Fig4';
if ~exist(rd4RESULT, 'dir')
    mkdir(rd4RESULT);
end

stations = {'T01','T02','ZSD02'};
Comps    = {'BHE','BHN','BHZ'};

st_1 = '20250314';
et_1 = '20250331';

% Colors
color_T01 = [0.0000 0.4470 0.7410];
color_T02 = [0.8500 0.3250 0.0980];
color_ZSD = [0.4660 0.6740 0.1880];

median_lw  = 2.5;
noise_lw   = 1.5;
face_alpha = 0.30;

% Noise model (NLNM / NHNM)
PL = [0.10 0.17 0.40 0.80 1.24 2.40 4.30 5.00 6.00 10.00 12.00 15.60 ...
      21.90 31.60 45.00 70.00 101.00 154.00 328.00 600.00 10000.00 100000.00];
AL = [-162.36 -166.70 -170.00 -166.40 -168.60 -159.98 -141.10 -72.36 -97.26 ...
      -132.18 -205.27 -37.65 -114.37 -160.58 -187.50 -216.47 -185.00 -168.34 ...
      -217.43 -258.28 -346.88 -346.88];
BL = [5.64 0.00 -8.30 28.90 52.48 29.81 0.00 -99.77 -66.49 -31.57 36.16 ...
      -104.33 -47.10 -16.28 0.00 15.70 0.00 -7.61 11.90 26.60 48.75 48.75];

PH = [0.10 0.22 0.32 0.80 3.80 4.60 6.30 7.90 15.40 20.00 345.80 100000.00];
AH = [-108.73 -150.34 -122.31 -116.85 -108.48 -74.66 0.66 -93.37 73.54 ...
      -151.52 -206.66 -206.66];
BH = [-17.23 -80.50 -23.87 32.51 18.08 -32.95 -127.18 -22.42 -162.98 10.01 ...
      34.63 31.63];

NLNM = zeros(size(PL));
fL    = zeros(size(PL));
for i = 1:length(PL)
    NLNM(i) = AL(i) + BL(i)*log10(PL(i));
    fL(i)   = 1 / PL(i);
end

NHNM = zeros(size(PH));
fH    = zeros(size(PH));
for i = 1:length(PH)
    NHNM(i) = AH(i) + BH(i)*log10(PH(i));
    fH(i)   = 1 / PH(i);
end

%% Main loop
for c = 1:numel(Comps)
    comp = Comps{c};
    fprintf('\nProcessing component: %s\n', comp);


    % Left panel statistics: whole-period PSD median + IQR
    PSDstat = struct();
    for s = 1:numel(stations)
        station = stations{s};
        switch station
            case 'T01'
                mainColor = color_T01;
                comp_use  = comp;
                psdDir    = rd4PSD_MAIN;
                p_name    = fullfile(psdDir, ...
                            sprintf('PSD_%s_%s_%s_%s.PSD', station, comp_use, st_1, et_1));
            case 'T02'
                mainColor = color_T02;
                comp_use  = comp;
                psdDir    = rd4PSD_MAIN;
                p_name    = fullfile(psdDir, ...
                            sprintf('PSD_%s_%s_%s_%s.PSD', station, comp_use, st_1, et_1));
            case 'ZSD02'
                mainColor = color_ZSD;
                comp_use  = strrep(comp, 'BH', 'HH');
                psdDir    = rd4PSD_ZSD;
                p_name    = fullfile(psdDir, ...
                            sprintf('PSD_ZSD02_%s_%s_%s.PSD', comp_use, st_1, et_1));
        end
        fprintf('Reading left-panel file: %s\n', p_name);
        if ~exist(p_name, 'file')
            warning('Missing file: %s', p_name);
            continue;
        end

        P_data = load(p_name);
        F_ref  = P_data(2:end, 1);
        PP     = P_data(2:end, 2:end);
        valid_cols = ~all(isnan(PP), 1);
        PP = PP(:, valid_cols);
        if isempty(PP)
            warning('No valid hourly PSD columns found in %s', p_name);
            continue;
        end

        PSDstat.(station).F      = F_ref;
        PSDstat.(station).PP     = PP;
        PSDstat.(station).median = median(PP, 2, 'omitnan');
        PSDstat.(station).p25    = prctile(PP', 25)';
        PSDstat.(station).p75    = prctile(PP', 75)';
        PSDstat.(station).color  = mainColor;
    end

    % Right panel: PSDdiff = PSD(T02) - PSD(T01)
    % median + IQR + zero line
    file1 = fullfile(rd4PSD_MAIN, sprintf('PSD_%s_%s_%s_%s.PSD', 'T01', comp, st_1, et_1));
    file2 = fullfile(rd4PSD_MAIN, sprintf('PSD_%s_%s_%s_%s.PSD', 'T02', comp, st_1, et_1));

    if isfile(file1) && isfile(file2)
        Pdata1 = load(file1);
        Pdata2 = load(file2);
        F_diff_full = Pdata1(2:end,1);
        PSD1 = Pdata1(2:end,2:end);
        PSD2 = Pdata2(2:end,2:end);
        ncol = min(size(PSD1,2), size(PSD2,2));
        PSD1 = PSD1(:,1:ncol);
        PSD2 = PSD2(:,1:ncol);
        valid_pair = ~(all(isnan(PSD1),1) | all(isnan(PSD2),1));
        PSD1 = PSD1(:, valid_pair);
        PSD2 = PSD2(:, valid_pair);

        if ~isempty(PSD1) && ~isempty(PSD2)
            PSDdiff = PSD2 - PSD1;

            medianPSD = median(PSDdiff, 2, 'omitnan');
            Q1        = prctile(PSDdiff, 25, 2);
            Q3        = prctile(PSDdiff, 75, 2);
            minPSD    = min(PSDdiff, [], 2, 'omitnan');
            maxPSD    = max(PSDdiff, [], 2, 'omitnan');

            % reduction probability
            valid_count  = sum(~isnan(PSDdiff), 2);
            reduce_count = sum(PSDdiff < 0, 2, 'omitnan');
            reduce_prob  = 100 * reduce_count ./ valid_count;
            reduce_prob(valid_count == 0) = NaN;

            step = 1;
            idx  = 1:step:length(F_diff_full);

            Fp_diff = F_diff_full(idx);
            medianP = medianPSD(idx);
            Q1p     = Q1(idx);
            Q3p     = Q3(idx);
            minP       = minPSD(idx);
            maxP       = maxPSD(idx);
            ProbReduce = reduce_prob(idx);
        end
    else
        warning('Whole-time PSD files missing for component %s. Right panel skipped.', comp);
    end


    %% Figure
    fig = figure('Visible','off');
    set(fig, 'Color', 'w', 'Position', [100 100 1800 460]);
    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % Left panel: median + IQR
    ax1 = nexttile(tl,1);
    hold(ax1,'on');
    order_left = {'T01','T02','ZSD02'};
    hMedian = gobjects(numel(order_left),1);

    % IQR patches first
    for k = 1:numel(order_left)
        stn = order_left{k};
        if ~isfield(PSDstat, stn)
            continue;
        end

        F_left = PSDstat.(stn).F;
        p25    = PSDstat.(stn).p25;
        p75    = PSDstat.(stn).p75;
        clr    = PSDstat.(stn).color;

        x_patch = [F_left; flipud(F_left)];
        y_patch = [p25; flipud(p75)];

        patch(ax1, x_patch, y_patch, clr, ...
            'FaceAlpha', face_alpha, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    % Noise model
    semilogx(ax1, fL, NLNM, '--', 'Color', [0.2 0.2 0.2], ...
        'LineWidth', noise_lw, 'HandleVisibility','off');
    semilogx(ax1, fH, NHNM, '--', 'Color', [0.2 0.2 0.2], ...
        'LineWidth', noise_lw, 'HandleVisibility','off');

    text(ax1, 0.012, -124, 'NHNM', 'FontSize', 17, ...
        'Color', [0.2 0.2 0.2], 'FontAngle', 'italic', 'FontName', 'Times New Roman');
    text(ax1, 0.012, -180, 'NLNM', 'FontSize', 17, ...
        'Color', [0.2 0.2 0.2], 'FontAngle', 'italic', 'FontName', 'Times New Roman');

    % Median curves last
    for k = 1:numel(order_left)
        stn = order_left{k};
        if ~isfield(PSDstat, stn)
            continue;
        end

        F_mean = PSDstat.(stn).F;
        mpsd   = PSDstat.(stn).median;
        clr    = PSDstat.(stn).color;

        hMedian(k) = semilogx(ax1, F_mean, mpsd, '-', ...
            'Color', clr, ...
            'LineWidth', median_lw, ...
            'DisplayName', stn);
    end

    for k = 1:numel(hMedian)
        if isgraphics(hMedian(k))
            uistack(hMedian(k), 'top');
        end
    end

    set(ax1, 'XScale', 'log', ...
        'XLim', [0.01 40], ...
        'YLim', [-190 -60], ...
        'XTick', [0.01 0.05 0.1 0.5 1 5 10 40], ...
        'XTickLabel', {'0.01','0.05','0.1','0.5','1','5','10','40'}, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 17, ...
        'LineWidth', 1.2, ...
        'TickDir', 'in', ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on');

    grid(ax1,'on');
    ax1.GridAlpha = 0.12;
    ax1.MinorGridAlpha = 0.08;
    xlabel(ax1, 'Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 20);
    ylabel(ax1, 'PSD 10log_{10}(m/s^2)^2/Hz dB', 'FontName', 'Times New Roman', 'FontSize', 20);

    valid_idx = isgraphics(hMedian);
    valid_handles = hMedian(valid_idx);
    valid_names = order_left(valid_idx);

    if strcmp(comp,'BHZ')
        legend(ax1, valid_handles, valid_names,'Location', 'northeast','Box', 'off', ...
            'FontName', 'Times New Roman', 'FontSize', 13);
    end

    % Right panel: median PSDdiff + IQR + zero line
    ax2 = nexttile(tl,2);
    hold(ax2,'on');

    % min-max error bars
    step_err = 2;
    idx_err  = 1:step_err:length(Fp_diff);
    lowerErr = medianP(idx_err) - minP(idx_err);
    upperErr = maxP(idx_err)    - medianP(idx_err);
    h_err = errorbar(ax2, Fp_diff(idx_err), medianP(idx_err), ...
        lowerErr, upperErr, 'LineStyle', 'none', ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.9, 'CapSize', 3, ...
        'DisplayName', 'Min–max range');

    % IQR patch
    x_patch2 = [Fp_diff; flipud(Fp_diff)];
    y_patch2 = [Q1p; flipud(Q3p)];
    h_iqr = patch(ax2, x_patch2, y_patch2, [0.55 0.55 0.55], ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none','DisplayName', 'IQR (25%–75%)');

    % Zero line
    yline(ax2, 0, '-', 'Color', [0.55 0.55 0.55],'LineWidth', 1.0, 'HandleVisibility', 'off');

    % Median line
    h_median_diff = semilogx(ax2, Fp_diff, medianP, '-', ...
        'Color',[0.80 0.05 0.05],'LineWidth',median_lw,'DisplayName','Median PSDdiff');

    if strcmp(comp,'BHZ')
        legend(ax2, h_median_diff, {'Median PSDdiff'}, ...
            'Location', 'northeast', 'Box', 'off', ...
            'FontName', 'Times New Roman', 'FontSize', 13);
    end
    set(ax2, 'XScale', 'log','XLim', [0.01 40],'XTick', [0.01 0.05 0.1 0.2 0.5 1 5 10 40], ...
        'XTickLabel', {'0.01','0.05','0.1','0.2','0.5','1','5','10','40'}, ...
        'FontName', 'Times New Roman','FontSize', 17, 'LineWidth', 1.2, ...
        'TickDir', 'in', 'Box', 'on', 'Layer', 'top','XMinorTick', 'off','YMinorTick', 'off');

    switch comp
        case {'BHE','BHN'}
            ylim(ax2, [-50 30]);
            yticks(ax2, -40:10:20);
        case 'BHZ'
            ylim(ax2, [-90 50]);
            yticks(ax2, -80:20:40);
    end

    % grid(ax2, 'on');
    % box(ax2, 'on');
    % ax2.GridLineStyle = '--';
    grid(ax2,'on');
    ax2.GridAlpha = 0.12;
    ax2.MinorGridAlpha = 0.08;

    xlabel(ax2, 'Frequency (Hz)', 'FontSize', 20, 'FontName', 'Times New Roman');
    ylabel(ax2, 'PSD Difference (dB)', 'FontSize', 20, 'FontName', 'Times New Roman');
    set(ax1, 'Position', [0.07 0.16 0.40 0.76]);
    set(ax2, 'Position', [0.56 0.16 0.40 0.76]);

    % Save
    output_png = fullfile(rd4RESULT, sprintf('Fig4_%s.tiff', comp));
    exportgraphics(fig, output_png, 'Resolution', 600);
    close(fig);

    %% Separate figure: reduction probability: black fraction curve
    have_reduction = true;
    if have_reduction
        fig_prob = figure('Visible','off');
        set(fig_prob, 'Color', 'w', 'Position', [100 100 800 300]);
        hold on;

        semilogx(Fp_diff, ProbReduce, '-', 'Color', 'r', 'LineWidth', 2.0);
        set(gca, 'XScale', 'log', 'XLim', [0.01 40], 'YLim', [0 100], ...
            'XTick', [0.05 0.1 0.5 1 5 10 40], ...
            'XTickLabel', {'0.05','0.1','0.5','1','5','10','40'},'YTick', 0:20:100, ...
            'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.0, ...
            'TickDir', 'in', 'TickLength', [0.018 0.018], ...
            'Box', 'on', 'Layer', 'top', 'XMinorTick', 'on', 'YMinorTick', 'on');
        grid on;
        axp = gca;
        axp.GridAlpha = 0.10;
        axp.MinorGridAlpha = 0.05;
        axp.GridLineStyle = '-';
        axp.MinorGridLineStyle = '-';
        xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 12);
        ylabel('Reduction probability (%)','FontName', 'Times New Roman', 'FontSize', 12);
        out_tif_prob = fullfile(rd4RESULT, sprintf('Reduction_probability_%s.tiff', comp));
        exportgraphics(fig_prob, out_tif_prob, 'Resolution', 600);
        close(fig_prob);
        fprintf('Saved reduction probability: %s\n', out_tif_prob);
    end
end

disp('Finished.');