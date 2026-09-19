clc; clear; close all;

% Repository-relative paths
script_file = mfilename('fullpath');
script_dir = fileparts(script_file);
repo_root = fileparts(script_dir);
project_root = fileparts(repo_root);
data_root = fullfile(project_root,'ESS_OBS_Shield_Data_v1.1.0');
rd4data = fullfile(data_root,'Pro_Seismic_acceleration_v1.0');
rd4psd = fullfile(data_root,'Processed_PSDDiff_v1.0');
rd4fig = fullfile(repo_root,'outputs','Figure_S2');
if ~exist(rd4fig,'dir'); mkdir(rd4fig); end

outfile = fullfile(rd4fig,'Figure_S2.tiff');

% Parameters
Fs = 100;
WL_hour = 3600;
WL_sub = 300;
fmin = 0.01;
fmax = 40;

freq_ticks = [0.01 0.05 0.1 0.2 0.5 1 5 10 20 40];
freq_labels = {'0.01','0.05','0.1','0.2','0.5','1','5','10','20','40'};

fontname = 'Times New Roman';
fontsize = 7.5;
panel_fontsize = 7.5;
legend_fontsize = 7;

colors = [0, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.9290, 0.6940, 0.1250];

c_zz = [0.4940 0.1840 0.5560];
c_zero = [0 0 0];

% Read SAC
disp('Reading SAC waveforms ...')

H11 = read_sac_series(rd4data,'ZS.T01.LC.BHN*.SAC');
H12 = read_sac_series(rd4data,'ZS.T01.LC.BHE*.SAC');
Z1 = read_sac_series(rd4data,'ZS.T01.LC.BHZ*.SAC');

H21 = read_sac_series(rd4data,'ZS.T02.LC.BHN*.SAC');
H22 = read_sac_series(rd4data,'ZS.T02.LC.BHE*.SAC');
Z2 = read_sac_series(rd4data,'ZS.T02.LC.BHZ*.SAC');

nTotal = min([numel(H11) numel(H12) numel(Z1) numel(H21) numel(H22) numel(Z2)]);

H11 = H11(1:nTotal);
H12 = H12(1:nTotal);
Z1 = Z1(1:nTotal);
H21 = H21(1:nTotal);
H22 = H22(1:nTotal);
Z2 = Z2(1:nTotal);

% Frequency
nfft = 2^nextpow2(WL_sub*Fs);
freq_all = (0:nfft/2)'/nfft*Fs;
fmask = freq_all>=fmin & freq_all<=fmax;
F = freq_all(fmask);

% Hourly windows
nHourSample = WL_hour*Fs;
nWin = floor(nTotal/nHourSample);
nFreq = numel(F);

fprintf('Hourly windows: %d\n',nWin);

% Coherence matrices
C_H1H2_T01 = nan(nFreq,nWin);
C_H1Z_T01 = nan(nFreq,nWin);
C_H2Z_T01 = nan(nFreq,nWin);

C_H1H2_T02 = nan(nFreq,nWin);
C_H1Z_T02 = nan(nFreq,nWin);
C_H2Z_T02 = nan(nFreq,nWin);

C_ZZ = nan(nFreq,nWin);

% Hourly coherence
disp('Computing hourly coherence ...')

for iw = 1:nWin
    id = (1:nHourSample)+(iw-1)*nHourSample;

    c = coh_hour(H11(id),H12(id),Fs,WL_sub);
    C_H1H2_T01(:,iw) = c(fmask);

    c = coh_hour(H11(id),Z1(id),Fs,WL_sub);
    C_H1Z_T01(:,iw) = c(fmask);

    c = coh_hour(H12(id),Z1(id),Fs,WL_sub);
    C_H2Z_T01(:,iw) = c(fmask);

    c = coh_hour(H21(id),H22(id),Fs,WL_sub);
    C_H1H2_T02(:,iw) = c(fmask);

    c = coh_hour(H21(id),Z2(id),Fs,WL_sub);
    C_H1Z_T02(:,iw) = c(fmask);

    c = coh_hour(H22(id),Z2(id),Fs,WL_sub);
    C_H2Z_T02(:,iw) = c(fmask);

    c = coh_hour(Z1(id),Z2(id),Fs,WL_sub);
    C_ZZ(:,iw) = c(fmask);

    if mod(iw,24)==0 || iw==nWin
        fprintf('  %d / %d hours finished\n',iw,nWin);
    end
end

% Coherence median
M_H1H2_T01 = median(C_H1H2_T01,2,'omitnan');
M_H1Z_T01 = median(C_H1Z_T01,2,'omitnan');
M_H2Z_T01 = median(C_H2Z_T01,2,'omitnan');

M_H1H2_T02 = median(C_H1H2_T02,2,'omitnan');
M_H1Z_T02 = median(C_H1Z_T02,2,'omitnan');
M_H2Z_T02 = median(C_H2Z_T02,2,'omitnan');

M_ZZ = median(C_ZZ,2,'omitnan');

% Read BHZ PSDdiff
disp('Reading vertical PSDdiff ...')

psdfile = fullfile(rd4psd,'DiffPSD_BHZ_20250314_20250331.PSD');
[fPSD,PSDdiff] = read_diffpsd(psdfile);

S_PSD = curve_stats(PSDdiff);

psd_med = interp1(fPSD,S_PSD.med,F,'linear',nan);
psd_q25 = interp1(fPSD,S_PSD.q25,F,'linear',nan);
psd_q75 = interp1(fPSD,S_PSD.q75,F,'linear',nan);

fprintf('PSDdiff loaded: %d frequencies x %d hourly records\n',size(PSDdiff,1),size(PSDdiff,2));

% Common PSDdiff axis
v = [psd_q25; psd_q75];
v = v(isfinite(v));

if isempty(v)
    psd_lim = 10;
else
    psd_lim = max(abs(v))*1.05;
    psd_lim = max(5,ceil(psd_lim/5)*5);
end

% Figure
set(groot,'defaultAxesFontName',fontname);
set(groot,'defaultTextFontName',fontname);

fig = figure('Units','centimeters','Position',[2 2 18 12.5],'Color','w');
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% (a) T01
ax1 = nexttile(tl,1);
hold(ax1,'on');

h11 = plot(ax1,F,M_H1H2_T01,'-','Color',colors(1,:),'LineWidth',1.1);
h12 = plot(ax1,F,M_H1Z_T01,'-','Color',colors(2,:),'LineWidth',1.1);
h13 = plot(ax1,F,M_H2Z_T01,'-','Color',colors(3,:),'LineWidth',1.1);

format_axis(ax1,freq_ticks,freq_labels,fmin,fmax,fontsize);
ylabel(ax1,'Coherence');
text(ax1,0.015,0.91,'(a) T01','Units','normalized', ...
    'FontSize',panel_fontsize,'FontWeight','normal');

lg1 = legend(ax1,[h11 h12 h13], ...
    {'HH1-HH2','HH1-HHZ','HH2-HHZ'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontSize',legend_fontsize, ...
    'NumColumns',1);
if isprop(lg1,'ItemTokenSize')
    lg1.ItemTokenSize = [10 6];
end

% (b) T02
ax2 = nexttile(tl,2);
hold(ax2,'on');

h21 = plot(ax2,F,M_H1H2_T02,'-','Color',colors(1,:),'LineWidth',1.1);
h22 = plot(ax2,F,M_H1Z_T02,'-','Color',colors(2,:),'LineWidth',1.1);
h23 = plot(ax2,F,M_H2Z_T02,'-','Color',colors(3,:),'LineWidth',1.1);

format_axis(ax2,freq_ticks,freq_labels,fmin,fmax,fontsize);
ylabel(ax2,'Coherence');
text(ax2,0.015,0.91,'(b) T02','Units','normalized', ...
    'FontSize',panel_fontsize,'FontWeight','normal');

% (c) T01-T02
ax3 = nexttile(tl,3);
hold(ax3,'on');

yyaxis(ax3,'right')
hp3 = plot_psd(F,psd_med,psd_q25,psd_q75);
yline(ax3,0,'--','Color',c_zero,'LineWidth',0.8,'HandleVisibility','off');
ylim(ax3,[-psd_lim psd_lim]);
ylabel(ax3,'PSD difference (dB)');

yyaxis(ax3,'left')
h33 = plot(ax3,F,M_ZZ,'-','Color',c_zz,'LineWidth',1.1);

format_axis(ax3,freq_ticks,freq_labels,fmin,fmax,fontsize);
ylabel(ax3,'Coherence');
xlabel(ax3,'Frequency (Hz)');
set_yyaxis_black(ax3);

text(ax3,0.015,0.91,'(c) T01-T02','Units','normalized', ...
    'FontSize',panel_fontsize,'FontWeight','normal');

lg3 = legend(ax3,[h33 hp3], ...
    {'HHZ-HHZ','PSD_{diff}'}, ...
    'Location','northeast', ...
    'Box','off', ...
    'FontSize',legend_fontsize, ...
    'NumColumns',1);
if isprop(lg3,'ItemTokenSize')
    lg3.ItemTokenSize = [10 6];
end


% Save TIFF
disp('Saving TIFF ...')
save_tiff_lzw(fig,outfile,600);
fprintf('Saved: %s\n',outfile);


% Functions

function data = read_sac_series(rd,pattern)

files = dir(fullfile(rd,pattern));
if isempty(files); error('No SAC files found: %s',fullfile(rd,pattern)); end

[~,idx] = sort({files.name});
files = files(idx);

data = [];
for i = 1:numel(files)
    [d,~] = rdSac(fullfile(rd,files(i).name));
    data = [data; d(:)];
end

end


function coh = coh_hour(x,y,Fs,WL_sub)

nSubSample = WL_sub*Fs;
nSub = floor(min(numel(x),numel(y))/nSubSample);
nfft = 2^nextpow2(nSubSample);
nf = nfft/2+1;
win = hann(nSubSample);

Sxx = zeros(nf,1);
Syy = zeros(nf,1);
Sxy = zeros(nf,1);
nValid = 0;

for k = 1:nSub
    id = (1:nSubSample)+(k-1)*nSubSample;
    xx = x(id);
    yy = y(id);

    if any(~isfinite(xx)) || any(~isfinite(yy)) || var(xx)==0 || var(yy)==0
        continue
    end

    X = fft(detrend(xx(:)).*win,nfft);
    Y = fft(detrend(yy(:)).*win,nfft);

    X = X(1:nf);
    Y = Y(1:nf);

    Sxx = Sxx+X.*conj(X);
    Syy = Syy+Y.*conj(Y);
    Sxy = Sxy+X.*conj(Y);
    nValid = nValid+1;
end

coh = nan(nf,1);
if nValid==0; return; end

Sxx = Sxx/nValid;
Syy = Syy/nValid;
Sxy = Sxy/nValid;

den = real(Sxx).*real(Syy);
valid = isfinite(den) & den>0;

coh(valid) = abs(Sxy(valid)).^2./den(valid);
coh(coh<0) = 0;
coh(coh>1) = 1;

end


function [freq,D] = read_diffpsd(file)

% PSD format:
% row 1    = hourly time
% column 1 = frequency
% body     = hourly PSDdiff

A = importdata(file);
if isstruct(A); A = A.data; end

if ~isnumeric(A) || isempty(A)
    error('Unable to read PSD file: %s',file);
end

if size(A,1)<2 || size(A,2)<2
    error('Invalid PSD matrix: %s',file);
end

freq = A(2:end,1);
D = A(2:end,2:end);

valid = isfinite(freq) & freq>0;
freq = freq(valid);
D = D(valid,:);

[freq,idx] = sort(freq);
D = D(idx,:);

[freq,idx] = unique(freq,'stable');
D = D(idx,:);

end


function S = curve_stats(M)

n = size(M,1);
S.med = nan(n,1);
S.q25 = nan(n,1);
S.q75 = nan(n,1);

for i = 1:n
    v = M(i,:);
    v = v(isfinite(v));

    if isempty(v); continue; end

    S.med(i) = median(v);
    q = prctile(v,[25 75]);
    S.q25(i) = q(1);
    S.q75(i) = q(2);
end

end


function h = plot_psd(F,med,q25,q75)

c_line = [0.80 0.05 0.05];
c_fill = [0.55 0.55 0.55];

valid = isfinite(F) & isfinite(med) & isfinite(q25) & isfinite(q75);
x = F(valid);
m = med(valid);
q1 = q25(valid);
q3 = q75(valid);

fill([x;flipud(x)],[q1;flipud(q3)],c_fill,'FaceAlpha',0.20,'EdgeColor','none','HandleVisibility','off');
h = plot(x,m,'-','Color',c_line,'LineWidth',1.1);

end


function format_axis(ax,ticks,labels,fmin,fmax,fontsize)
set(ax,'XScale','log','XLim',[fmin fmax],'YLim',[0 1]);
set(ax,'XTick',ticks,'XTickLabel',labels,'YTick',0:0.2:1);
set(ax,'FontName','Times New Roman','FontSize',fontsize,'LineWidth',0.7);
set(ax,'TickDir','in','TickLength',[0.012 0.012],'Box','on');
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.12;
grid(ax,'on');
end


function set_yyaxis_black(ax)

ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

end


function save_tiff_lzw(fig,outfile,dpi)

tmp = [tempname '.tif'];

set(fig,'PaperPositionMode','auto');
print(fig,tmp,'-dtiff',['-r' num2str(dpi)]);

img = imread(tmp);
imwrite(img,outfile,'tif','Compression','lzw','Resolution',dpi);

if exist(tmp,'file'); delete(tmp); end

end


function [data,hd] = rdSac(sacFile)

fid = fopen(sacFile,'r');
if fid<0; error('Cannot open SAC file: %s',sacFile); end

hd = fread(fid,70,'single');
hd(71:158) = fread(fid,88,'int');
data = fread(fid,inf,'single');

fclose(fid);

end
