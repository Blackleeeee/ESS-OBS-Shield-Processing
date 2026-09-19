% Plot PSD and DiffPSD spectrograms with tide-level curve
% Yang Li

clc; clear;

% Repository-relative paths
script_file = mfilename('fullpath');
script_dir = fileparts(script_file);
repo_root = fileparts(script_dir);
project_root = fileparts(repo_root);
data_root = fullfile(project_root,'ESS_OBS_Shield_Data_v1.1.0');
rd4PSD = fullfile(data_root,'Processed_PSD_v1.0');
rd4DIFF = fullfile(data_root,'Processed_PSDDiff_v1.0');
rd4OUT = fullfile(repo_root,'outputs','Figure_S1');
if ~exist(rd4OUT,'dir'); mkdir(rd4OUT); end

components = ["BHE"; "BHN"; "BHZ"];
stations = ["T01", "T02"];

% Time range
st = '20250314';
et = '20250331';

% Figure style
FIG_W_CM = 18.0;
FIG_H_CM = 5.0;

FONT_NAME = 'Times New Roman';
FS_TICK = 6.5;
FS_LABEL = 7.0;
FS_CBAR = 7.0;
AX_LW = 0.75;

% Tide
tide_file = fullfile(data_root,'Raw_Current_Tide_v1.1','Tide_UTC.txt');
opts = detectImportOptions(tide_file);
opts = setvartype(opts, {'Time'}, 'char');
tide_raw = readtable(tide_file, opts);
Tide_Time = datenum(tide_raw.Time, 'yyyy-mm-dd HH:MM:SS');
Tide_Level = tide_raw.waterlevel / 100;

%% Plot raw PSD
for s = 1:length(stations)
    for c = 1:length(components)
        station = stations(s);
        component = components(c);

        cd(rd4PSD)
        PSD_file = strcat('PSD_', station, '_', component, '_', st, '_', et, '.PSD');
        PSD = load(PSD_file);
        PSD_time = PSD(1, 2:end);

        [Fig1, ax1, ch] = TFsurfaceDF(PSD, 'PSD', FIG_W_CM, FIG_H_CM, FONT_NAME, FS_TICK, FS_LABEL, FS_CBAR, AX_LW);

        hold(ax1, 'on')
        yyaxis(ax1, 'right')
        plot(ax1, Tide_Time, Tide_Level, 'k-', 'LineWidth', 1.0)
        ylabel(ax1, 'Tide level (m)', 'FontSize', FS_LABEL, 'FontName', FONT_NAME, 'Color', 'k')
        ylim(ax1, [0 4.5])
        xlim(ax1, [PSD_time(1), PSD_time(end)])

        yticks(ax1, 0:1:4)
        ax1.YAxis(2).Color = 'k';
        ax1.YAxis(2).FontSize = FS_TICK;
        ax1.YAxis(2).LineWidth = AX_LW;
        ax1.YAxis(2).MinorTick = 'on';
        try
            ax1.YAxis(2).MinorTickValues = 0.5:0.5:4.5;
        catch
        end

        yyaxis(ax1, 'left')
        ax1.YAxis(1).Color = 'k';

        cd(rd4OUT)
        outbase = strcat('TFS2TIDE_', station, '_', component, '_', st, '_', et);
        print(Fig1, strcat(outbase, '.tiff'), '-dtiff', '-r300');
        close(Fig1)
    end
end

%% Plot DiffPSD
for c = 1:length(components)
    component = components(c);

    cd(rd4DIFF)
    PSD_file = strcat('DiffPSD_', component, '_', st, '_', et, '.PSD');
    PSD = load(PSD_file);
    PSD_time = PSD(1, 2:end);

    [Fig2, ax2, ch] = TFsurfaceDF(PSD, 'DIFF', FIG_W_CM, FIG_H_CM, FONT_NAME, FS_TICK, FS_LABEL, FS_CBAR, AX_LW);

    hold(ax2, 'on')
    yyaxis(ax2, 'right')
    plot(ax2, Tide_Time, Tide_Level, 'k-', 'LineWidth', 1.0)
    ylabel(ax2, 'Tide level (m)', 'FontSize', FS_LABEL, 'FontName', FONT_NAME, 'Color', 'k')
    ylim(ax2, [0 4.5])
    xlim(ax2, [PSD_time(1), PSD_time(end)])

    yticks(ax2, 0:1:4)
    ax2.YAxis(2).Color = 'k';
    ax2.YAxis(2).FontSize = FS_TICK;
    ax2.YAxis(2).LineWidth = AX_LW;
    ax2.YAxis(2).MinorTick = 'on';
    try
        ax2.YAxis(2).MinorTickValues = 0.5:0.5:4.5;
    catch
    end

    yyaxis(ax2, 'left')
    ax2.YAxis(1).Color = 'k';

    cd(rd4OUT)
    outbase = strcat('Figure_S1_', component, '_', st, '_', et);
    print(Fig2, strcat(outbase, '.tiff'), '-dtiff', '-r300');
    close(Fig2)
end

%% Function
function [Fig, ax, ch] = TFsurfaceDF(data, plot_type, fig_w_cm, fig_h_cm, font_name, fs_tick, fs_label, fs_cbar, ax_lw)

    T = data(1, 2:end);
    Fa = data(2:end, 1);

    DFl = 0.01;
    DFh = 40.0;
    [~, FLI] = min(abs(Fa - DFl));
    [~, FHI] = min(abs(Fa - DFh));

    F = data(FLI+1:FHI+1, 1);
    S = data(FLI+1:FHI+1, 2:end);

    if length(T) > 1000
        n = length(T);
    else
        n = 1000;
    end

    % Dense frequency
    f = 1:length(F);
    ff = linspace(1, length(F), n);
    FF = interp1(f, F, ff, 'linear');
    FFF = repmat(FF(:), 1, n);

    % Dense time
    t = 1:length(T);
    tt = linspace(1, length(T), n);
    TT = interp1(t, T, tt, 'linear');
    TTT = repmat(TT, n, 1);

    % Dense S
    IM = 'linear';
    if size(data, 1) > size(data, 2)
        TF = linspace(min(T), max(T), size(S, 1));
        ss = zeros(length(F), length(TF));
        for i = 1:length(F)
            ss(i, :) = interp1(T, S(i, :), TF, IM);
        end
        TR = repmat(TF, length(F), 1);
        FR = repmat(F, 1, length(TF));
        SS = interp2(TR, FR, ss, TTT, FFF);
    elseif size(data, 1) < size(data, 2)
        FT = linspace(min(F), max(F), size(S, 2));
        ss = zeros(length(FT), length(T));
        for i = 1:length(T)
            ss(:, i) = interp1(F, S(:, i), FT, IM);
        end
        TR = repmat(T, length(FT), 1);
        FR = repmat(FT(:), 1, length(T));
        SS = interp2(TR, FR, ss, TTT, FFF);
    else
        SS = S;
        TT = T;
        FF = F;
    end

    % Figure
    Fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2, 2, fig_w_cm, fig_h_cm]);

    ax = axes('Parent', Fig, ...
        'Position', [0.05, 0.08, 0.84, 0.88], ...
        'YScale', 'log', ...
        'YMinorTick', 'on', ...
        'YTick', [0.01 0.05 0.1 0.5 1.0 5.0 10.0 40.0], ...
        'YTickLabel', {'0.01', '0.05', '0.1', '0.5', '1.0', '5.0', '10.0', '40.0'}, ...
        'FontSize', fs_tick, ...
        'FontName', font_name, ...
        'LineWidth', ax_lw, ...
        'TickDir', 'in', ...
        'Box', 'on', ...
        'Layer', 'top');

    surface(ax, TT, FF, SS, 'LineStyle', 'none');
    view(ax, 2)
    axis(ax, [min(T) max(T) DFl DFh])

    ylabel(ax, 'Frequency (Hz)', 'FontSize', fs_label, 'FontName', font_name)
    xlabel(ax, '')

    ax.TickLength = [0.010 0.010];
    ax.XColor = 'k';
    ax.YColor = 'k';

    % X axis
    x0 = ceil(min(T));
    x1 = floor(max(T));
    xt_major = x0+1:2:x1;
    ax.XTick = xt_major;
    ax.XTickLabel = cellstr(datestr(xt_major, 'mmm.dd'));
    ax.XTickLabelRotation = 0;
    ax.XMinorTick = 'on';
    try
        ax.XRuler.MinorTickValues = x0:0.5:x1;
    catch
    end

    % Colorbar
    ch = colorbar(ax, 'eastoutside');
    ch.FontSize = fs_tick;
    ch.FontName = font_name;
    ch.LineWidth = ax_lw;
    ch.TickDirection = 'in';
    ch.TickLength = 0.04;
    ch.Position = [0.924, 0.08, 0.01, 0.88];

    if strcmpi(plot_type, 'PSD')
        colormap(Fig, jet)
        caxis(ax, [-180 -50])
        ylabel(ch, 'PSD 10log_{10}(m/s^2)^2/Hz dB', 'FontSize', fs_cbar, 'FontName', font_name)

    elseif strcmpi(plot_type, 'DIFF')
        cmap = bluewhitered(256);
        colormap(Fig, cmap)

        % maxAbs = max(abs(SS(:)), [], 'omitnan');
        % if ~isfinite(maxAbs) || maxAbs == 0
        %     maxAbs = 1;
        % end
        % caxis(ax, [-maxAbs maxAbs])
        caxis(ax, [-30 30])

        ylabel(ch, 'PSD difference (dB)', 'FontSize', fs_cbar, 'FontName', font_name)
    end
end

function cmap = bluewhitered(n)
    if mod(n, 2) ~= 0
        n = n + 1;
    end
    m = n / 2;
    bluePart = [linspace(0,1,m)' linspace(0,1,m)' ones(m,1)];
    redPart  = [ones(m,1) linspace(1,0,m)' linspace(1,0,m)'];
    cmap = [bluePart; redPart];
end
