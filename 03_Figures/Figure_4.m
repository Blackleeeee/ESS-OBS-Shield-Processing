clc; clear; close all;

% PSD difference versus current speed for HHZ and HH.
% A single 2 x 3 figure is produced: HHZ on the first row and HH on the second.
%
% HH is calculated by summing the two horizontal PSDs in linear power:
%   PSD_HH = 10*log10(10.^(PSD_BHE/10) + 10.^(PSD_BHN/10))
% The configuration difference is then calculated in dB:
%   PSDdiff = PSD_T02 - PSD_T01

% Repository-relative paths
script_file = mfilename('fullpath');
script_dir = fileparts(script_file);
repo_root = fileparts(script_dir);
project_root = fileparts(repo_root);
data_root = fullfile(project_root,'ESS_OBS_Shield_Data_v1.1.0');
rd4PSD = fullfile(data_root,'Processed_PSD_v1.0');
rd4Curr = fullfile(data_root,'Raw_Current_Tide_v1.1');
floodFile = fullfile(rd4Curr,'FloodEvents.txt');
ebbFile = fullfile(rd4Curr,'EbbEvents.txt');
rd4Out = fullfile(repo_root,'outputs','Figure_4');
if ~exist(rd4Out,'dir'); mkdir(rd4Out); end

% Components and period
plotComps = {'HHZ','HH'};
st = '20250314';
et = '20250331';

% Frequency bands
Bands = {
    '0.01–0.1 Hz', 0.01, 0.1;
    '0.1–10 Hz',   0.1,  10;
    '10–40 Hz',    10,   40
};

% Current-speed bins
speed_bin_width = 0.05;
speed_edges = 0:speed_bin_width:0.85;
speed_centers = (speed_edges(1:end-1)+speed_edges(2:end))/2;
minN = 30;

% Plot settings
col_flood = [0.8500 0.3250 0.0980];
col_ebb = [0.4940 0.1840 0.5560];
line_width = 1.0;
marker_size = 3.2;
shade_alpha = 0.14;

fontname = 'Times New Roman';
fontsize = 7.0;
label_fontsize = 7.0;
panel_fontsize = 7.0;
legend_fontsize = 6.5;

xlim_use = [0 0.8];
xticks_use = 0:0.2:0.8;
xminor_ticks = 0:0.1:0.8;
ylim_use = [-40 20];
yticks_use = -40:10:20;

% Current speed
speed_file = fullfile(rd4Curr,'Flow_Speed_UTC.txt');
opts = detectImportOptions(speed_file,'Delimiter','\t');
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];

flow_raw = readtable(speed_file,opts);
Time_Flow = datenum(datetime(flow_raw{:,1},'InputFormat','yyyy-MM-dd HH:mm:ss'));
Speed_Flow = flow_raw{:,2};

[Time_Flow,ia] = unique(Time_Flow);
Speed_Flow = Speed_Flow(ia);

% Flood and ebb intervals
[FloodStart,FloodEnd] = read_tidal_events(floodFile);
[EbbStart,EbbEnd] = read_tidal_events(ebbFile);

Results = struct();
Stats = struct();

% Figure: slightly larger outer margins and tighter horizontal spacing
fig = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[2 2 18 7.5],'PaperPositionMode','auto');
axs = gobjects(2,3);

left_margin = 0.090;
right_margin = 0.070;
horizontal_gap = 0.02;
ax_width = (1-left_margin-right_margin-2*horizontal_gap)/3;
ax_height = 0.36;
row_bottom = [0.52 0.10];

for c = 1:numel(plotComps)
    plotComp = plotComps{c};
    fprintf('\nProcessing component: %s\n',plotComp);

    [PSD_freq,PSD_time,PSDdiff] = build_psddiff(plotComp,rd4PSD,st,et);
    Speed_psd = interp1(Time_Flow,Speed_Flow,PSD_time,'linear',NaN);
    isFlood = in_periods(PSD_time,FloodStart,FloodEnd);
    isEbb = in_periods(PSD_time,EbbStart,EbbEnd);

    Results.(plotComp).PSD_time = PSD_time;
    Results.(plotComp).PSD_freq = PSD_freq;
    Results.(plotComp).PSDdiff = PSDdiff;
    Results.(plotComp).Speed_psd = Speed_psd;
    Results.(plotComp).isFlood = isFlood;
    Results.(plotComp).isEbb = isEbb;

    for b = 1:size(Bands,1)
        ax_left = left_margin+(b-1)*(ax_width+horizontal_gap);
        axs(c,b) = axes(fig,'Units','normalized', ...
            'Position',[ax_left row_bottom(c) ax_width ax_height]);
        ax = axs(c,b);
        hold(ax,'on');

        f1 = Bands{b,2};
        f2 = Bands{b,3};
        [xFlood,yFlood] = extract_freqtime_samples(PSD_freq,PSDdiff,Speed_psd,isFlood,f1,f2);
        [xEbb,yEbb] = extract_freqtime_samples(PSD_freq,PSDdiff,Speed_psd,isEbb,f1,f2);

        Stats.(plotComp)(b).Flood = plot_binned_median_iqr_xy( ...
            ax,xFlood,yFlood,speed_edges,col_flood,minN,line_width,marker_size,shade_alpha);
        Stats.(plotComp)(b).Ebb = plot_binned_median_iqr_xy( ...
            ax,xEbb,yEbb,speed_edges,col_ebb,minN,line_width,marker_size,shade_alpha);
        Stats.(plotComp)(b).BandName = Bands{b,1};
        Stats.(plotComp)(b).f1 = f1;
        Stats.(plotComp)(b).f2 = f2;

        yline(ax,0,'--','Color','k','LineWidth',0.8,'HandleVisibility','off');
        set(ax,'XLim',xlim_use,'YLim',ylim_use,'XTick',xticks_use,'YTick',yticks_use, ...
            'FontName',fontname,'FontSize',fontsize,'LineWidth',0.7,'TickDir','in', ...
            'TickLength',[0.014 0.014],'Box','on','Layer','top', ...
            'XMinorTick','on','YMinorTick','off','XColor','k','YColor','k');

        try
            ax.XAxis.MinorTickValues = xminor_ticks;
        catch
        end
        ax.XGrid = 'on';
        ax.XMinorGrid = 'off';
        ax.YGrid = 'off';
        ax.YMinorGrid = 'off';
        ax.GridLineStyle = '--';
        ax.GridColor = [0.65 0.65 0.65];
        ax.GridAlpha = 0.45;

        if c == 1
            % title(ax,Bands{b,1},'FontName',fontname,'FontSize',title_fontsize,'FontWeight','normal');
            set(ax,'XTickLabel',[]);
        else
            xlabel(ax,'Current speed (m/s)','FontName',fontname,'FontSize',label_fontsize);
        end

        if b == 1
            ylabel(ax,[plotComp ' PSD_{\rm diff} (dB)'], ...
                'Interpreter','tex', ...
                'FontName',fontname, ...
                'FontSize',label_fontsize);
        else
            set(ax,'YTickLabel',[]);
        end
    end

end

% Shared legend in the first panel
h1 = plot(axs(1,1),nan,nan,'o-','Color',col_flood, ...
    'MarkerFaceColor',col_flood,'MarkerEdgeColor','k', ...
    'MarkerSize',marker_size,'LineWidth',line_width);
h2 = plot(axs(1,1),nan,nan,'o-','Color',col_ebb, ...
    'MarkerFaceColor',col_ebb,'MarkerEdgeColor','k', ...
    'MarkerSize',marker_size,'LineWidth',line_width);
lgd = legend(axs(1,1),[h1 h2],{'Flood tide','Ebb tide'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontName',fontname, ...
    'FontSize',legend_fontsize);
lgd.ItemTokenSize = [10 7];

out_tif = fullfile(rd4Out,'Figure_4.tiff');
exportgraphics(fig,out_tif,'Resolution',600);
close(fig);
fprintf('Saved figure: %s\n',out_tif);

out_mat = fullfile(rd4Out,sprintf( 'FreqTimeScatter_DPSD_vs_Current_HHZ_HH_%s_%s.mat',st,et));
save(out_mat,'Results','Bands','Stats','speed_edges','speed_centers', 'speed_bin_width','minN');
fprintf('\nSaved statistics: %s\n',out_mat);
disp('Finished.');


function [F,time,PSDdiff] = build_psddiff(plotComp,rd4PSD,st,et)

[F1,time1,PSD1] = load_station_psd('T01',plotComp,rd4PSD,st,et);
[F2,time2,PSD2] = load_station_psd('T02',plotComp,rd4PSD,st,et);
assert_same_frequency(F1,F2,sprintf('T01 and T02 %s',plotComp));

[time,idx1,idx2] = intersect(time1,time2,'stable');
if isempty(time)
    error('No common hourly PSD columns for T01 and T02 %s.',plotComp);
end

PSD1 = PSD1(:,idx1);
PSD2 = PSD2(:,idx2);
valid_pair = ~(all(~isfinite(PSD1),1) | all(~isfinite(PSD2),1));

time = time(valid_pair);
PSD1 = PSD1(:,valid_pair);
PSD2 = PSD2(:,valid_pair);
if isempty(time)
    error('No valid paired hourly PSD columns for %s.',plotComp);
end

F = F1;
PSDdiff = PSD2-PSD1;

end


function [F,time,PSDdB] = load_station_psd(station,plotComp,rd4PSD,st,et)

if strcmp(plotComp,'HH')
    channels = {'BHE','BHN'};
elseif strcmp(plotComp,'HHZ')
    channels = {'BHZ'};
else
    error('Unsupported component: %s',plotComp);
end

file1 = fullfile(rd4PSD,sprintf('PSD_%s_%s_%s_%s.PSD', ...
    station,channels{1},st,et));
if ~isfile(file1); error('Missing file: %s',file1); end
fprintf('Reading: %s\n',file1);
[F,time,PSDdB] = read_psd_file(file1);

if numel(channels) == 1
    return
end

file2 = fullfile(rd4PSD,sprintf('PSD_%s_%s_%s_%s.PSD', ...
    station,channels{2},st,et));
if ~isfile(file2); error('Missing file: %s',file2); end
fprintf('Reading: %s\n',file2);
[F2,time2,PSD2dB] = read_psd_file(file2);

assert_same_frequency(F,F2,sprintf('%s horizontal components',station));
[time,idx1,idx2] = intersect(time,time2,'stable');
if isempty(time)
    error('No common hourly PSD columns between BHE and BHN for %s.',station);
end

PSDdB = horizontal_total_db(PSDdB(:,idx1),PSD2dB(:,idx2));

end


function [F,time,PSDdB] = read_psd_file(filename)

data = load(filename);
if size(data,1)<2 || size(data,2)<2
    error('Invalid PSD matrix: %s',filename);
end

F = data(2:end,1);
time = data(1,2:end);
PSDdB = data(2:end,2:end);

end


function PSD_HH_dB = horizontal_total_db(PSD1_dB,PSD2_dB)

if ~isequal(size(PSD1_dB),size(PSD2_dB))
    error('Horizontal PSD matrices have different sizes after time alignment.');
end

PSD_HH_dB = nan(size(PSD1_dB));
valid = isfinite(PSD1_dB) & isfinite(PSD2_dB);
PSD_HH_dB(valid) = 10*log10( ...
    10.^(PSD1_dB(valid)/10) + 10.^(PSD2_dB(valid)/10));

end


function assert_same_frequency(F1,F2,label)

if numel(F1) ~= numel(F2)
    error('Frequency-grid length mismatch for %s.',label);
end

tolerance = max(1e-12,1e-8*max(abs(F1)));
if any(abs(F1-F2)>tolerance)
    error('Frequency-grid mismatch for %s.',label);
end

end


function [StartTime,EndTime] = read_tidal_events(eventFile)

opts = detectImportOptions(eventFile,'Delimiter','\t');
opts.VariableTypes = {'string','datetime','datetime'};
opts = setvaropts(opts,2,'InputFormat','yyyy-MM-dd HH:mm:ss');
opts = setvaropts(opts,3,'InputFormat','yyyy-MM-dd HH:mm:ss');

Tab = readtable(eventFile,opts);
StartTime = datenum(Tab{:,2});
EndTime = datenum(Tab{:,3});

end


function mask = in_periods(T,StartTime,EndTime)

mask = false(size(T));
for i = 1:length(StartTime)
    mask = mask | (T>=StartTime(i) & T<EndTime(i));
end

end


function [x,y] = extract_freqtime_samples(PSD_freq,PSDdiff,Speed,idxStage,f1,f2)

fidx = PSD_freq>=f1 & PSD_freq<=f2;
if ~any(fidx)
    x = [];
    y = [];
    return
end

tidx = idxStage & isfinite(Speed);
if ~any(tidx)
    x = [];
    y = [];
    return
end

D = PSDdiff(fidx,tidx);
S = Speed(tidx);
S_mat = repmat(S(:).',size(D,1),1);

x = S_mat(:);
y = D(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

end


function Stat = plot_binned_median_iqr_xy( ...
    ax,X,Y,edges,color_rgb,minN,line_width,marker_size,shade_alpha)

nB = length(edges)-1;
xc = nan(1,nB);
medv = nan(1,nB);
q1v = nan(1,nB);
q3v = nan(1,nB);
nv = zeros(1,nB);

for i = 1:nB
    idb = X>=edges(i) & X<edges(i+1) & isfinite(X) & isfinite(Y);
    nv(i) = sum(idb);
    if nv(i)<minN
        continue
    end

    vals = Y(idb);
    xc(i) = (edges(i)+edges(i+1))/2;
    medv(i) = median(vals,'omitnan');
    q1v(i) = prctile(vals,25);
    q3v(i) = prctile(vals,75);
end

Stat.xc = xc;
Stat.median = medv;
Stat.q1 = q1v;
Stat.q3 = q3v;
Stat.N = nv;

valid = isfinite(xc) & isfinite(medv);
if sum(valid)<2
    return
end

x = xc(valid);
med = medv(valid);
q1 = q1v(valid);
q3 = q3v(valid);

patch(ax,[x fliplr(x)],[q1 fliplr(q3)],color_rgb, ...
    'FaceAlpha',shade_alpha,'EdgeColor','none','HandleVisibility','off');
plot(ax,x,med,'-','Color',color_rgb,'LineWidth',line_width, ...
    'HandleVisibility','off');
plot(ax,x,med,'o','MarkerFaceColor',color_rgb,'MarkerEdgeColor','k', ...
    'MarkerSize',marker_size,'LineWidth',0.4,'HandleVisibility','off');

end
