% Figure 3: horizontal-total (HH) and vertical (HHZ) PSD
% Left column: whole-period PSD median and IQR for T01, T02, and ZSD02.
% Right column: PSDdiff = PSD(T02) - PSD(T01), median and IQR.
%
% HH is the rotation-invariant horizontal total power. Because the input
% files contain PSD in dB, the two horizontal PSDs are converted to linear
% power, summed, and converted back to dB:
%   PSD_HH = 10*log10(10.^(PSD_E/10) + 10.^(PSD_N/10)).

clc;
clear;
close all;

% Repository-relative paths
script_file = mfilename('fullpath');
script_dir = fileparts(script_file);
repo_root = fileparts(script_dir);
project_root = fileparts(repo_root);
data_root = fullfile(project_root,'ESS_OBS_Shield_Data_v1.1.0');
rd4PSD_MAIN = fullfile(data_root,'Processed_PSD_v1.0');
rd4PSD_ZSD  = rd4PSD_MAIN;
rd4RESULT   = fullfile(repo_root,'outputs','Figure_3');
if ~exist(rd4RESULT, 'dir')
    mkdir(rd4RESULT);
end

stations = {'T01','T02','ZSD02'};
plotComps = {'HH','HHZ'};
st_1 = '20250314';
et_1 = '20250331';

% Colors
color_T01 = [0.0000 0.4470 0.7410];
color_T02 = [0.8500 0.3250 0.0980];
color_ZSD = [0.4660 0.6740 0.1880];

median_lw  = 2.5;
noise_lw   = 1.5;
face_alpha = 0.30;

% Noise models (NLNM and NHNM)
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

NLNM = AL + BL .* log10(PL);
NHNM = AH + BH .* log10(PH);
fL = 1 ./ PL;
fH = 1 ./ PH;

for c = 1:numel(plotComps)
    plotComp = plotComps{c};
    fprintf('\nProcessing component: %s\n', plotComp);

    % Load and summarize whole-period PSD
    PSDstat = struct();
    for s = 1:numel(stations)
        station = stations{s};
        switch station
            case 'T01'
                mainColor = color_T01;
            case 'T02'
                mainColor = color_T02;
            case 'ZSD02'
                mainColor = color_ZSD;
        end

        try
            [F_ref, time_ref, PP] = load_plot_psd(station, plotComp, rd4PSD_MAIN, ...
                rd4PSD_ZSD, st_1, et_1);
        catch ME
            warning('%s', ME.message);
            continue;
        end

        valid_cols = ~all(isnan(PP), 1);
        PP = PP(:, valid_cols);
        time_ref = time_ref(valid_cols);
        if isempty(PP)
            warning('No valid hourly PSD columns found for %s %s.', station, plotComp);
            continue;
        end

        PSDstat.(station).F = F_ref;
        PSDstat.(station).time = time_ref;
        PSDstat.(station).PP = PP;
        PSDstat.(station).median = median(PP, 2, 'omitnan');
        PSDstat.(station).p25 = prctile(PP', 25)';
        PSDstat.(station).p75 = prctile(PP', 75)';
        PSDstat.(station).color = mainColor;
    end

    % T01 and T02 have already been combined to HH before differencing.
    have_diff = isfield(PSDstat, 'T01') && isfield(PSDstat, 'T02');
    if have_diff
        assert_same_frequency(PSDstat.T01.F, PSDstat.T02.F, 'T01 and T02');
        [common_time, idx1, idx2] = intersect(PSDstat.T01.time, PSDstat.T02.time, 'stable');
        if isempty(common_time)
            warning('No common hourly columns for T01 and T02 %s.', plotComp);
            have_diff = false;
        else
            PSD1 = PSDstat.T01.PP(:, idx1);
            PSD2 = PSDstat.T02.PP(:, idx2);
            valid_pair = ~(all(isnan(PSD1), 1) | all(isnan(PSD2), 1));
            PSD1 = PSD1(:, valid_pair);
            PSD2 = PSD2(:, valid_pair);
            if isempty(PSD1)
                warning('No valid paired hourly PSD columns for %s.', plotComp);
                have_diff = false;
            end
        end
    end

    if have_diff
        Fp_diff = PSDstat.T01.F;
        PSDdiff = PSD2 - PSD1;
        medianP = median(PSDdiff, 2, 'omitnan');
        Q1p = prctile(PSDdiff', 25)';
        Q3p = prctile(PSDdiff', 75)';
        minP = min(PSDdiff, [], 2, 'omitnan');
        maxP = max(PSDdiff, [], 2, 'omitnan');

        valid_count = sum(~isnan(PSDdiff), 2);
        reduce_count = sum(PSDdiff < 0, 2);
        ProbReduce = 100 * reduce_count ./ valid_count;
        ProbReduce(valid_count == 0) = NaN;
    end

    % Figure: absolute PSD and PSD difference
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1800 460]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl, 1);
    hold(ax1, 'on');
    order_left = {'T01','T02','ZSD02'};
    hMedian = gobjects(numel(order_left), 1);

    for k = 1:numel(order_left)
        stn = order_left{k};
        if ~isfield(PSDstat, stn)
            continue;
        end
        F_left = PSDstat.(stn).F;
        p25 = PSDstat.(stn).p25;
        p75 = PSDstat.(stn).p75;
        clr = PSDstat.(stn).color;
        patch(ax1, [F_left; flipud(F_left)], [p25; flipud(p75)], clr, ...
            'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    semilogx(ax1, fL, NLNM, '--', 'Color', [0.2 0.2 0.2], ...
        'LineWidth', noise_lw, 'HandleVisibility', 'off');
    semilogx(ax1, fH, NHNM, '--', 'Color', [0.2 0.2 0.2], ...
        'LineWidth', noise_lw, 'HandleVisibility', 'off');
    text(ax1, 0.012, -124, 'NHNM', 'FontSize', 17, 'Color', [0.2 0.2 0.2], ...
        'FontAngle', 'italic', 'FontName', 'Times New Roman');
    text(ax1, 0.012, -180, 'NLNM', 'FontSize', 17, 'Color', [0.2 0.2 0.2], ...
        'FontAngle', 'italic', 'FontName', 'Times New Roman');

    for k = 1:numel(order_left)
        stn = order_left{k};
        if ~isfield(PSDstat, stn)
            continue;
        end
        hMedian(k) = semilogx(ax1, PSDstat.(stn).F, PSDstat.(stn).median, '-', ...
            'Color', PSDstat.(stn).color, 'LineWidth', median_lw, 'DisplayName', stn);
    end
    for k = 1:numel(hMedian)
        if isgraphics(hMedian(k))
            uistack(hMedian(k), 'top');
        end
    end

    set(ax1, 'XScale', 'log', 'XLim', [0.01 40], 'YLim', [-190 -60], ...
        'XTick', [0.01 0.05 0.1 0.5 1 5 10 40], ...
        'XTickLabel', {'0.01','0.05','0.1','0.5','1','5','10','40'}, ...
        'FontName', 'Times New Roman', 'FontSize', 17, 'LineWidth', 1.2, ...
        'TickDir', 'in', 'Box', 'on', 'Layer', 'top', ...
        'XMinorTick', 'on', 'YMinorTick', 'on');
    grid(ax1, 'on');
    ax1.GridAlpha = 0.12;
    ax1.MinorGridAlpha = 0.08;
    xlabel(ax1, 'Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 20);
    ylabel(ax1, 'PSD 10log_{10}(m/s^2)^2/Hz dB', ...
        'FontName', 'Times New Roman', 'FontSize', 20);

    valid_idx = isgraphics(hMedian);
    if strcmp(plotComp, 'HHZ') && any(valid_idx)
        legend(ax1, hMedian(valid_idx), order_left(valid_idx), 'Location', 'northeast', ...
            'Box', 'off', 'FontName', 'Times New Roman', 'FontSize', 13);
    end

    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    if have_diff
        step_err = 2;
        idx_err = 1:step_err:numel(Fp_diff);
        lowerErr = medianP(idx_err) - minP(idx_err);
        upperErr = maxP(idx_err) - medianP(idx_err);
        errorbar(ax2, Fp_diff(idx_err), medianP(idx_err), lowerErr, upperErr, ...
            'LineStyle', 'none', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.9, ...
            'CapSize', 3, 'HandleVisibility', 'off');
        patch(ax2, [Fp_diff; flipud(Fp_diff)], [Q1p; flipud(Q3p)], [0.55 0.55 0.55], ...
            'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        yline(ax2, 0, '-', 'Color', [0.55 0.55 0.55], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
        h_median_diff = semilogx(ax2, Fp_diff, medianP, '-', ...
            'Color', [0.80 0.05 0.05], 'LineWidth', median_lw, ...
            'DisplayName', 'Median PSDdiff');
        if strcmp(plotComp, 'HHZ')
            legend(ax2, h_median_diff, {'Median PSDdiff'}, 'Location', 'northeast', ...
                'Box', 'off', 'FontName', 'Times New Roman', 'FontSize', 13);
        end
    else
        text(ax2, 0.5, 0.5, 'No paired T01/T02 PSD data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontName', 'Times New Roman');
    end

    set(ax2, 'XScale', 'log', 'XLim', [0.01 40], ...
        'XTick', [0.01 0.05 0.1 0.2 0.5 1 5 10 40], ...
        'XTickLabel', {'0.01','0.05','0.1','0.2','0.5','1','5','10','40'}, ...
        'FontName', 'Times New Roman', 'FontSize', 17, 'LineWidth', 1.2, ...
        'TickDir', 'in', 'Box', 'on', 'Layer', 'top', ...
        'XMinorTick', 'off', 'YMinorTick', 'off');
    if strcmp(plotComp, 'HH')
        ylim(ax2, [-50 30]);
        yticks(ax2, -40:10:20);
    else
        ylim(ax2, [-90 50]);
        yticks(ax2, -80:20:40);
    end
    grid(ax2, 'on');
    ax2.GridAlpha = 0.12;
    ax2.MinorGridAlpha = 0.08;
    xlabel(ax2, 'Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 20);
    ylabel(ax2, 'PSD Difference (dB)', 'FontName', 'Times New Roman', 'FontSize', 20);
    set(ax1, 'Position', [0.07 0.16 0.40 0.76]);
    set(ax2, 'Position', [0.56 0.16 0.40 0.76]);

    output_tif = fullfile(rd4RESULT, sprintf('Figure_3_%s.tiff', plotComp));
    exportgraphics(fig, output_tif, 'Resolution', 600);
    close(fig);
    fprintf('Saved: %s\n', output_tif);

    if have_diff
        fig_prob = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 300]);
        semilogx(Fp_diff, ProbReduce, '-', 'Color', 'r', 'LineWidth', 2.0);
        set(gca, 'XScale', 'log', 'XLim', [0.01 40], 'YLim', [0 100], ...
            'XTick', [0.05 0.1 0.5 1 5 10 40], ...
            'XTickLabel', {'0.05','0.1','0.5','1','5','10','40'}, ...
            'YTick', 0:20:100, 'FontName', 'Times New Roman', 'FontSize', 10, ...
            'LineWidth', 1.0, 'TickDir', 'in', 'TickLength', [0.018 0.018], ...
            'Box', 'on', 'Layer', 'top', 'XMinorTick', 'on', 'YMinorTick', 'on');
        grid on;
        axp = gca;
        axp.GridAlpha = 0.10;
        axp.MinorGridAlpha = 0.05;
        xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 12);
        ylabel('Reduction probability (%)', 'FontName', 'Times New Roman', 'FontSize', 12);
        out_prob = fullfile(rd4RESULT, sprintf('Reduction_probability_%s.tiff', plotComp));
        exportgraphics(fig_prob, out_prob, 'Resolution', 600);
        close(fig_prob);
        fprintf('Saved reduction probability: %s\n', out_prob);
    end
end

disp('Finished.');


function [F, time, PSDdB] = load_plot_psd(station, plotComp, mainDir, zsdDir, st, et)
% Load one plotted component. HH is formed in linear power.
    if strcmp(plotComp, 'HH')
        if strcmp(station, 'ZSD02')
            channels = {'HHE','HHN'};
        else
            channels = {'BHE','BHN'};
        end
    elseif strcmp(plotComp, 'HHZ')
        if strcmp(station, 'ZSD02')
            channels = {'HHZ'};
        else
            channels = {'BHZ'};
        end
    else
        error('Unsupported plotted component: %s', plotComp);
    end

    files = cell(size(channels));
    for i = 1:numel(channels)
        if strcmp(station, 'ZSD02')
            files{i} = fullfile(zsdDir, sprintf('PSD_ZSD02_%s_%s_%s.PSD', ...
                channels{i}, st, et));
        else
            files{i} = fullfile(mainDir, sprintf('PSD_%s_%s_%s_%s.PSD', ...
                station, channels{i}, st, et));
        end
        if ~isfile(files{i})
            error('Missing file: %s', files{i});
        end
        fprintf('Reading: %s\n', files{i});
    end

    [F, time, PSDdB] = read_psd_file(files{1});
    if numel(files) == 1
        return;
    end

    [F2, time2, PSD2dB] = read_psd_file(files{2});
    assert_same_frequency(F, F2, sprintf('%s horizontal components', station));
    [time, idx1, idx2] = intersect(time, time2, 'stable');
    if isempty(time)
        error('No common hourly columns between the two horizontal files for %s.', station);
    end
    PSDdB = horizontal_total_db(PSDdB(:, idx1), PSD2dB(:, idx2));
end


function [F, time, PSDdB] = read_psd_file(filename)
    data = load(filename);
    if size(data, 1) < 2 || size(data, 2) < 2
        error('Invalid PSD matrix: %s', filename);
    end
    F = data(2:end, 1);
    time = data(1, 2:end);
    PSDdB = data(2:end, 2:end);
end


function PSD_HH_dB = horizontal_total_db(PSD1_dB, PSD2_dB)
% Sum orthogonal horizontal PSDs in linear power, not in dB.
    if ~isequal(size(PSD1_dB), size(PSD2_dB))
        error('Horizontal PSD matrices have different sizes after time alignment.');
    end
    PSD_HH_dB = nan(size(PSD1_dB));
    valid = isfinite(PSD1_dB) & isfinite(PSD2_dB);
    PSD_HH_dB(valid) = 10 .* log10(10.^(PSD1_dB(valid) ./ 10) + ...
        10.^(PSD2_dB(valid) ./ 10));
end


function assert_same_frequency(F1, F2, label)
    if numel(F1) ~= numel(F2)
        error('Frequency-grid length mismatch for %s.', label);
    end
    tolerance = max(1e-12, 1e-8 * max(abs(F1)));
    if any(abs(F1 - F2) > tolerance)
        error('Frequency-grid mismatch for %s.', label);
    end
end

