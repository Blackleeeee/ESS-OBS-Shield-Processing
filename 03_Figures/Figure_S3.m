% Figure S3: T01/T02 PSD for joint current-direction and tidal-stage groups.
% Pool frequency-time points in each speed bin; median and IQR (not a CI).
% HH sums horizontal PSDs in linear power. Hour counts set the sample gate.
% Separate HHZ and HH figures: rows T01/T02, columns three frequency bands.
clc; clear; close all;

% Repository-relative paths
script_file = mfilename('fullpath');
script_dir = fileparts(script_file);
repo_root = fileparts(script_dir);
project_root = fileparts(repo_root);
data_root = fullfile(project_root,'ESS_OBS_Shield_Data_v1.1.0');
rdPSD = fullfile(data_root,'Processed_PSD_v1.0');
rdCurrent = fullfile(data_root,'Raw_Current_Tide_v1.1');
speedFile = fullfile(rdCurrent,'Flow_Speed_UTC.txt');
directionFile = fullfile(rdCurrent,'Flow_Direction_UTC.txt');
floodFile = fullfile(rdCurrent,'FloodEvents.txt');
ebbFile = fullfile(rdCurrent,'EbbEvents.txt');
rdOut = fullfile(repo_root,'outputs','Figure_S3');
st = '20250314'; et = '20250331';

% Time alignment: confirm against the script that generated the hourly PSD.
% 'start': [timestamp,timestamp+1h); 'center': +/-30min; 'end': preceding hour.
cfg.psdTimestamp = 'start';
cfg.windowSeconds = 3600;
cfg.currentSampleSeconds = 10;  % 0.1 Hz, used to assess hourly coverage
cfg.minCoverage = 0.80;
cfg.directionCenters = [120 300];
cfg.directionHalfWidth = 45;    % degrees about each observed direction
cfg.minDirectionFraction = 0.80; % fraction of valid samples within one sector
cfg.slackSpeed = 0.05;          % ignore unreliable direction at weak currents
cfg.minStageFraction = 0.80;    % mixed/unclassified hours are excluded from joint groups
cfg.speedEdges = 0:0.05:0.85;
cfg.minHoursPerDirection = 5;   % HOURS, unlike original minN frequency points
cfg.bands = [0.01 0.1; 0.1 10; 10 40];
cfg.bandNames = {'0.01-0.1 Hz','0.1-10 Hz','10-40 Hz'};
cfg.groupNames = ["Flood tide","Ebb tide"];
cfg.groupDirections = [300 120];
cfg.colors = [0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560]; % flood, ebb
cfg.figureSizeCm = [18 7.5];
cfg.xLim = [0 0.8];

% Read environmental records; match timestamps, never interpolate through gaps.
[ts,vs] = read_environment(speedFile);
[td,vd] = read_environment(directionFile);
[tc,is,id] = intersect(ts,td);
if isempty(tc); error('Speed and direction have no common timestamps.'); end
speed = vs(is); direction = vd(id);
direction(direction<0 | direction>360) = NaN;
direction = mod(direction,360);
speed(speed<0) = NaN;
[fs,fe] = read_events(floodFile);
[es,ee] = read_events(ebbFile);
if ~exist(rdOut,'dir'); mkdir(rdOut); end
fprintf('PSD timestamps treated as %s of %.0f-second windows.\n',cfg.psdTimestamp,cfg.windowSeconds);
fprintf('Using %d paired speed/direction records; no environmental interpolation.\n',numel(tc));

allHourly = table(); allBins = table();
components = {'HHZ','HH'};
for c = 1:2
    comp = components{c};
    [f,t1,p1] = load_psd(rdPSD,'T01',comp,st,et);
    [f2,t2,p2] = load_psd(rdPSD,'T02',comp,st,et);
    same_frequency(f,f2);
    [t,i1,i2] = intersect(t1,t2);
    if isempty(t); error('No paired PSD times for %s.',comp); end
    p1 = p1(:,i1); p2 = p2(:,i2);

    % Analysis period: 15–30 March 2025
    ta = round(datenum(datetime(2025,3,15,0,0,0))*86400);
    tb = round(datenum(datetime(2025,3,31,0,0,0))*86400);  % half-open end
    keep = t>=ta & t<tb;
    t = t(keep);
    p1 = p1(:,keep);
    p2 = p2(:,keep);

    E = hourly_environment(t,tc,speed,direction,fs,fe,es,ee,cfg);
    E.JointGroup = repmat("Excluded",height(E),1);
    for g = 1:2
        j = E.EnvironmentOK & E.DirectionGroup_deg==cfg.groupDirections(g) & E.Stage==cfg.groupNames(g);
        E.JointGroup(j) = cfg.groupNames(g);
    end
    E.JointSelectionOK = E.JointGroup~="Excluded";
    H = table(); B = table();
    for b = 1:3
        fi = f>=cfg.bands(b,1) & f<=cfg.bands(b,2);
        if ~any(fi); error('No frequency points in %s.',cfg.bandNames{b}); end
        Q1 = p1(fi,:); Q2 = p2(fi,:);
        complete = all(isfinite(Q1) & isfinite(Q2),1);
        Q1(:,~complete) = NaN; Q2(:,~complete) = NaN;
        T = E;
        T.Component = repmat(string(comp),height(T),1);
        T.Band = repmat(string(cfg.bandNames{b}),height(T),1);
        T.NFrequencyPoints = repmat(sum(fi),height(T),1);
        T.T01_HourlyMedianPSD_dB = median(Q1,1,'omitnan')';
        T.T02_HourlyMedianPSD_dB = median(Q2,1,'omitnan')';
        T.UseForComparison = T.JointSelectionOK & complete(:);
        S = bin_statistics(T,{Q1,Q2},cfg);
        H = [H; T]; B = [B; S];
    end
    allHourly = [allHourly; H]; allBins = [allBins; B];
    fprintf('%s: paired %d h; joint flood %d h; joint ebb %d h; excluded %d h before PSD QC.\n', ...
        comp,numel(t),sum(E.JointGroup=="Flood tide"),sum(E.JointGroup=="Ebb tide"),sum(~E.JointSelectionOK));
end
for c = 1:2
    make_figure(allBins,components{c},cfg,rdOut);
end
writetable(allHourly,fullfile(rdOut,'Figure_S3_Hourly_QC.csv'));
writetable(allBins,fullfile(rdOut,'Figure_S3_PooledPSD_Statistics.csv'));
save(fullfile(rdOut,'Figure_S3_Statistics.mat'),'cfg','allHourly','allBins');
fprintf('Saved Figure_S3_HHZ.tiff, Figure_S3_HH.tiff and statistics to:\n%s\n',rdOut);

function [t,v] = read_environment(file)
    if ~isfile(file); error('Missing environmental file: %s',file); end
    o = detectImportOptions(file,'Delimiter','\t');
    o.VariableNamesLine = 1; o.DataLines = [2 Inf];
    o = setvartype(o,o.VariableNames{1},'string');
    A = readtable(file,o);
    t = datetime(string(A{:,1}),'InputFormat','yyyy-MM-dd HH:mm:ss');
    v = A{:,2};
    if ~isnumeric(v); v = str2double(string(v)); end
    t = round(datenum(t)*86400); v = double(v(:));
    good = isfinite(t); t = t(good); v = v(good);
    [tu,~,group] = unique(t);
    nDuplicate = numel(t)-numel(tu);
    if nDuplicate>0
        % Each timestamp contributes once. Preserve an available finite value;
        % inconsistent finite values are missing, never averaged arbitrarily.
        finite = isfinite(v);
        vmin = accumarray(group(finite),v(finite),[numel(tu) 1],@min,NaN);
        vmax = accumarray(group(finite),v(finite),[numel(tu) 1],@max,NaN);
        conflict = isfinite(vmin) & isfinite(vmax) & abs(vmax-vmin)>1e-10;
        v = vmin; v(conflict) = NaN;
        fprintf('%s: merged %d duplicate rows; %d conflicting timestamps set to NaN.\n', ...
            file,nDuplicate,sum(conflict));
    else
        [~,j] = sort(t); v = v(j);
    end
    t = tu;
end

function [s,e] = read_events(file)
    o = detectImportOptions(file,'Delimiter','\t');
    o = setvartype(o,o.VariableNames(2:3),'string');
    A = readtable(file,o);
    s = round(datenum(datetime(string(A{:,2}),'InputFormat','yyyy-MM-dd HH:mm:ss'))*86400);
    e = round(datenum(datetime(string(A{:,3}),'InputFormat','yyyy-MM-dd HH:mm:ss'))*86400);
    if any(e<=s); error('Invalid tidal intervals in %s.',file); end
    [s,j] = sort(s); e = e(j);
    if any(s(2:end)<e(1:end-1)); error('Overlapping intervals in %s.',file); end
end

function [f,t,p] = load_psd(folder,station,comp,st,et)
    if strcmp(comp,'HHZ'); ch = {'BHZ'}; else; ch = {'BHE','BHN'}; end
    [f,t,p] = read_psd(fullfile(folder,sprintf('PSD_%s_%s_%s_%s.PSD',station,ch{1},st,et)));
    if numel(ch)==2
        [f2,t2,p2] = read_psd(fullfile(folder,sprintf('PSD_%s_%s_%s_%s.PSD',station,ch{2},st,et)));
        same_frequency(f,f2); [t,i,j] = intersect(t,t2);
        a = p(:,i); b = p2(:,j); m = max(a,b);
        p = m+10*log10(10.^((a-m)/10)+10.^((b-m)/10));
        p(~isfinite(a)|~isfinite(b)) = NaN;
    end
end

function [f,t,p] = read_psd(file)
    A = load(file);
    if size(A,1)<2 || size(A,2)<2; error('Invalid PSD matrix: %s',file); end
    f = A(2:end,1); t = A(1,2:end)'; p = A(2:end,2:end);
    if any(~isfinite(t)) || any(t<datenum(2025,1,1) | t>datenum(2026,1,1))
        error('Expected MATLAB datenums for 2025 in PSD row 1: %s',file);
    end
    if any(~isfinite(f)) || any(diff(f)<=0); error('Frequency grid must increase: %s',file); end
    t = round(t*86400);
    if numel(unique(t))~=numel(t); error('Duplicate PSD timestamps: %s',file); end
    [t,j] = sort(t); p = p(:,j);
end

function same_frequency(f,g)
    if numel(f)~=numel(g) || any(abs(f-g)>max(1e-12,1e-8*max(abs(f))))
        error('PSD frequency grids differ. No automatic resampling is performed.');
    end
end

function E = hourly_environment(t,tc,u,d,fs,fe,es,ee,cfg)
    switch cfg.psdTimestamp
        case 'start'; a = t;
        case 'center'; a = t-cfg.windowSeconds/2;
        case 'end'; a = t-cfg.windowSeconds;
        otherwise; error('psdTimestamp must be start, center or end.');
    end
    z = a+cfg.windowSeconds; n = numel(t);
    meanSpeed = nan(n,1); meanDir = nan(n,1); coverage = zeros(n,1);
    frac = zeros(n,2); group = nan(n,1); stage = strings(n,1);
    rise = zeros(n,1); fall = zeros(n,1); validN = zeros(n,1);
    for k = 1:n
        j = tc>=a(k) & tc<z(k) & isfinite(u) & isfinite(d);
        validN(k) = sum(j); coverage(k) = min(1,validN(k)*cfg.currentSampleSeconds/cfg.windowSeconds);
        rise(k) = sum(max(0,min(z(k),fe)-max(a(k),fs)))/cfg.windowSeconds;
        fall(k) = sum(max(0,min(z(k),ee)-max(a(k),es)))/cfg.windowSeconds;
        if rise(k)+fall(k)>1+1e-6; error('Flood/ebb event files overlap in a PSD window.'); end
        stage(k) = "Other";
        if rise(k)>=cfg.minStageFraction; stage(k) = "Flood tide";
        elseif fall(k)>=cfg.minStageFraction; stage(k) = "Ebb tide"; end
        if ~any(j); continue; end
        uk = u(j); dk = d(j); meanSpeed(k) = mean(uk);
        moving = uk>=cfg.slackSpeed;
        if any(moving)
            meanDir(k) = mod(atan2d(mean(sind(dk(moving))),mean(cosd(dk(moving)))),360);
        end
        for g = 1:2
            angle = abs(mod(dk-cfg.directionCenters(g)+180,360)-180);
            frac(k,g) = sum(moving & angle<=cfg.directionHalfWidth)/numel(uk);
        end
        [v,g] = max(frac(k,:));
        if coverage(k)>=cfg.minCoverage && v>=cfg.minDirectionFraction
            group(k) = cfg.directionCenters(g);
        end
    end
    E = table(datetime(t/86400,'ConvertFrom','datenum'),datetime(a/86400,'ConvertFrom','datenum'), ...
        datetime(z/86400,'ConvertFrom','datenum'),meanSpeed,meanDir,group,coverage,validN, ...
        frac(:,1),frac(:,2),stage,rise,fall,isfinite(group), ...
        'VariableNames',{'PSDTime_UTC','WindowStart_UTC','WindowEnd_UTC','Speed_ms', ...
        'MeanDirection_deg','DirectionGroup_deg','Coverage','NCurrentSamples', ...
        'Sector1Fraction','Sector2Fraction','Stage','RisingFraction','FallingFraction','EnvironmentOK'});
end

function S = bin_statistics(T,Q,cfg)
    rows = {}; edges = cfg.speedEdges; stations = ["T01","T02"];
    for k = 1:numel(edges)-1
        in = T.UseForComparison & T.Speed_ms>=edges(k) & T.Speed_ms<edges(k+1);
        masks = {in & T.JointGroup==cfg.groupNames(1),in & T.JointGroup==cfg.groupNames(2)};
        counts = [sum(masks{1}) sum(masks{2})];
        common = all(counts>=cfg.minHoursPerDirection);
        for m = 1:2
            P = Q{m};
            for g = 1:2
                j = masks{g}; values = P(:,j); values = values(isfinite(values));
                med = NaN; q = [NaN NaN];
                if ~isempty(values); med = median(values); q = prctile(values,[25 75]); end
                rows(end+1,:) = {T.Component(1),T.Band(1),stations(m),cfg.groupNames(g), ...
                    cfg.groupDirections(g),edges(k),edges(k+1),counts(g),numel(values), ...
                    mean(T.Speed_ms(j)),med,q(1),q(2),common}; %#ok<AGROW>
            end
        end
    end
    S = cell2table(rows,'VariableNames',{'Component','Band','Station','TidalGroup', ...
        'Direction_deg','SpeedLow','SpeedHigh','NHours','NSpectralPoints', ...
        'MeanSpeed_ms','Median_dB','Q25_dB','Q75_dB','CommonBin'});
end

function make_figure(B,comp,cfg,out)
    stations = ["T01","T02"];
    marks = {'^','v'};
    styles = {'-','--'};
    fig = figure('Visible','off','Color','w','Units','centimeters', 'Position',[2 2 cfg.figureSizeCm],'PaperPositionMode','auto');
    ax = gobjects(2,3);

    % Layout control
    left   = 0.095;
    right  = 0.075;
    bottom = 0.10;
    top    = 0.08;

    hgap = 0.015;
    vgap = 0.040;

    w  = (1-left-right-2*hgap)/3;
    ht = (1-bottom-top-vgap)/2;
    rowBottom = [bottom+ht+vgap, bottom];
    B = B(B.Component==string(comp),:);

    for m = 1:2
        for b = 1:3
            xpos = left + (b-1)*(w+hgap);
            ypos = rowBottom(m);

            ax(m,b) = axes(fig,'Units','normalized', 'Position',[xpos ypos w ht]);
            a = ax(m,b);
            hold(a,'on');

            for g = 1:2
                R = B(B.Station==stations(m) & ...
                      B.Band==string(cfg.bandNames{b}) & ...
                      B.TidalGroup==cfg.groupNames(g),:);

                x  = (R.SpeedLow + R.SpeedHigh)/2;
                y  = R.Median_dB;
                lo = R.Q25_dB;
                hi = R.Q75_dB;

                y(~R.CommonBin)  = NaN;
                lo(~R.CommonBin) = NaN;
                hi(~R.CommonBin) = NaN;

                col = cfg.colors(g,:);

                errorbar(a,x,y,y-lo,hi-y, ...
                    'LineStyle',styles{g}, ...
                    'Marker',marks{g}, ...
                    'Color',col, ...
                    'MarkerFaceColor',col, ...
                    'MarkerEdgeColor',col, ...
                    'MarkerSize',3.2, ...
                    'LineWidth',1, ...
                    'CapSize',2, ...
                    'HandleVisibility','off');
            end

            set(a,'XLim',cfg.xLim, ...
                'XTick',0.2:0.2:0.6, ...
                'FontName','Times New Roman', ...
                'FontSize',6.5, ...
                'LineWidth',0.5, ...
                'TickDir','in', ...
                'TickLength',[0.014 0.014], ...
                'Box','on', ...
                'Layer','top', ...
                'XMinorTick','on', ...
                'YMinorTick','off', ...
                'XColor','k', ...
                'YColor','k');
            a.XGrid = 'on';
            a.YGrid = 'off';
            a.YMinorGrid = 'off';
            a.GridLineStyle = '--';
            a.GridColor = [0.65 0.65 0.65];
            a.GridAlpha = 0.45;

            if m == 1
                a.XTickLabel = [];
            else
                xlabel(a,'Current speed (m/s)','FontSize',7);
            end

            if b == 1
                ylabel(a,sprintf('%s PSD (dB)',char(stations(m))),'FontSize',7);
            else
                a.YTickLabel = [];
            end

            label = sprintf('(%c)','a'+(m-1)*3+b-1);
            if m == 1
                label = [label ' ' cfg.bandNames{b}];
            end

            text(a,0.035,0.95,label,'Units','normalized', ...
                'VerticalAlignment','top', ...
                'FontName','Times New Roman', ...
                'FontSize',6.5);
        end
    end

    % Shared y limits across all panels, including the full IQR
    xc = (B.SpeedLow + B.SpeedHigh)/2;
    R = B(B.CommonBin & xc>=cfg.xLim(1) & xc<=cfg.xLim(2),:);
    v = [R.Q25_dB; R.Q75_dB];
    v = v(isfinite(v));

    if ~isempty(v)
        ymin = min(v);
        ymax = max(v);
        pad = 0.05*(ymax-ymin);
        if pad == 0
            pad = 1;
        end
        set(ax(:),'YLim',[ymin-pad ymax+pad],'YTickMode','auto');
    end

    linkaxes(ax(:),'x');
    h = gobjects(1,2);
    for g = 1:2
        h(g) = plot(ax(1,1),nan,nan, ...
            'LineStyle',styles{g}, ...
            'Marker',marks{g}, ...
            'Color',cfg.colors(g,:), ...
            'MarkerFaceColor',cfg.colors(g,:), ...
            'MarkerSize',3.0, ...
            'LineWidth',1);
    end
    lg = legend(ax(1,1),h,{'Flood tide','Ebb tide'}, ...
        'Location','southeast', ...
        'NumColumns',1, ...
        'Box','off', ...
        'FontName','Times New Roman', ...
        'FontSize',6, ...
        'AutoUpdate','off');
    lg.ItemTokenSize = [8 6];
    exportgraphics(fig,fullfile(out,['Figure_S3_' comp '.tiff']),'Resolution',600);
    close(fig);
end
