close all
clear all

% code also uses:
% cbrewer (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
% m_map (https://www.eoas.ubc.ca/~rich/map.html)
% SEM_calc.m (https://github.com/thertz/HertzLabCode/blob/master/Matlab/General/SEM_calc.m)
% cmocean (https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)

% define results file (i.e., where results files are stored)
fid  = '/Users/riwoolway/Google Drive/Manuscripts/ice_thickness/results';

% custom colors for figures
ireds = cbrewer('seq','Reds',9);
set1 = cbrewer('qual','Set1',9);
ispec = cbrewer('div','RdYlBu',11);
rcp_cols = {set1(2,:),set1(5,:),set1(1,:)};

% continued
CT = cbrewer('div', 'RdYlBu', 11);
CT2 = cbrewer('div', 'RdYlGn', 11);
CT3 = cbrewer('qual', 'Set2', 8);
CT4 = cbrewer('div', 'PuOr', 11);
CT6 = cbrewer('seq', 'Purples', 12);
major_cols = [CT(2,:);CT(4,:);CT2(9,:);CT4(9,:);[0 0 0]];
iminor_str = {'Af','Am','As','Aw',...
    'BWh','BWk','BSh','BSk',...
    'Csa','Csb','Csc','Cwa','Cwb','Cwc','Cfa','Cfb','Cfc',...
    'Dsa','Dsb','Dsc','Dsd','Dwa','Dwb','Dwc','Dwd','Dfa','Dfb','Dfc','Dfd',...
    'ET','EF'};
wwint = flipud(winter(9));
aaut = autumn(8);
ckoepp = [aaut(1:4,:);...
    aaut(5:8,:);...
    wwint;...
    CT6;...
    [.6 .6 .6];[.7 .7 .7]];

% population count data
myncfile = [fid '/' 'global_population_at2020.nc'];
ilat = ncread(myncfile,'lat');
ilon = ncread(myncfile,'lon');
iilon = mod((ilon+180),360)-180;
pop_count = ncread(myncfile,'Population');

% population density
myncfile = [fid '/' 'global_population_density_at2020_remapcon2.nc'];
pop_dens = ncread(myncfile,'Population_Density');

% koppen classification
myncfile = [fid '/' 'koppen_climate_classification_clm5grid.nc'];
ikopp = ncread(myncfile,'type');
ikopp2 = ikopp;
ikopp2(ikopp <= 31) = 5;
ikopp2(ikopp <= 29) = 4;
ikopp2(ikopp <= 17) = 3;
ikopp2(ikopp <= 8) = 2;
ikopp2(ikopp <= 4) = 1;

% surface area of each grid
myncfile = [fid '/' 'grid_cell_area.nc'];
cell_area = ncread(myncfile,'area');

% sort lon -180 to 180
[ilon2,idx_lon] = sort(iilon);
pop_count2 = pop_count(idx_lon,:);
pop_dens2 = pop_dens(idx_lon,:);
cell_area2 = cell_area(idx_lon,:);
ikopp2 = ikopp2(idx_lon,:);
 
% only select grids above 30N
idx_lat = find(ilat >= 30);
ilat2 = ilat(idx_lat);
pop_count2 = pop_count2(:,idx_lat);
pop_dens2 = pop_dens2(:,idx_lat);
cell_area2 = cell_area2(:,idx_lat);
ikopp2b = ikopp2(:,idx_lat);

% population density bins
pop_bins = [0 1; 1 5; 5 10; 10 20; 20 50; 50 100;...
    100 200; 200 500; 500 1e3; 1e3 1e99];

% global air temperature
ifid = fopen([fid '/' 'trefht_glob_mean2.csv']);
dd = textscan(ifid,'%f%f%f%f%f%f%f','delimiter',',','headerlines',1);
fclose(ifid);
iyear = dd{1};
idx_year = find(iyear >= 1900 & iyear <= 1929);
itas = dd{2};
itas = itas - nanmean(itas(idx_year));
itas_arctic = nanmean([dd{6} dd{7}],2);
itas_arctic = itas_arctic - nanmean(itas_arctic(idx_year));

% rank air temperature
[rank_itas,idx_itas] = sort(itas);
    
% indices for future warming levels
idx_1 = find(itas > 1.5,1,'first');
idx_2 = find(itas > 2,1,'first');
idx_3 = find(itas > 3,1,'first');

% 15-year window
n_pre = idx_year;
n1 = (idx_1 - 7):(idx_1 + 7);
n2 = (idx_2 - 7):(idx_2 + 7);
n3 = (idx_3 - 7):(idx_3 + 7);

%
k = 3; % thickness == 4 inches
  
    % ice duration time series
    myncfile = [fid '/' 'Iceduration_global_average.nc'];
    ice_dur_ts = ncread(myncfile,'iceduration',[1 1 k],[Inf Inf 1]);
        
    % calculate anomalies
    ice_dur_ts2 = ice_dur_ts;
    ice_dur_ts_anom2 = ice_dur_ts;
    for ii = 1:size(ice_dur_ts,2)
        
        %
        iice_dur_ts = ice_dur_ts(:,ii);
        
        % moving average
        ice_dur_ts2(:,ii) = movmean(iice_dur_ts(idx_itas),31);
            
        % anomaly
        ice_dur_ts_anom = iice_dur_ts - nanmean(iice_dur_ts(idx_year));
        
        % moving average
        ice_dur_ts_anom2(:,ii) = movmean(ice_dur_ts_anom(idx_itas),31);
    end

    %
    av_ice_dur_ts2 = nanmean(ice_dur_ts2,2);
    
    %
    av_ice_dur_ts_anom2 = nanmean(ice_dur_ts_anom2,2);    
                
    % load spatially resolved ice metrics
    myncfile = [fid '/' 'Iceduration_ensmean_2.nc'];
    ice_dur = ncread(myncfile,'iceduration',[1 1 1 k],[Inf Inf Inf 1]);
    
    myncfile = [fid '/' 'Iceduration_ensstd.nc'];
    ice_dur_std = ncread(myncfile,'iceduration',[1 1 1 k],[Inf Inf Inf 1]);
    
    %
    ice_dur2 = ice_dur(idx_lon,:,:,:);
    ice_dur2 = ice_dur2(:,idx_lat,:,:);
    
    ice_dur2_std = ice_dur_std(idx_lon,:,:,:);
    ice_dur2_std = ice_dur2_std(:,idx_lat,:,:);

    % calculate pre-industrial average
    ice_dur2_pre = nanmean(ice_dur2(:,:,n_pre),3);
    ice_dur2_std_pre = nanmean(ice_dur2_std(:,:,n_pre),3);

    % if no ice during pre-industrial, set to NA    
    ice_dur3 = ice_dur2;
    ice_dur3_std = ice_dur2_std;
    for ii = 1:size(ice_dur2,3)        
        
        %
        ice_dur2b = ice_dur2(:,:,ii);
        ice_dur2b(ice_dur2_pre == 0) = nan;
        ice_dur3(:,:,ii) = ice_dur2b;
        
        %
        ice_dur2b_std = ice_dur2_std(:,:,ii);
        ice_dur2b_std(ice_dur2_std_pre == 0) = nan;
        ice_dur3_std(:,:,ii) = ice_dur2b_std;
    end
    
    % averages
    ice_dur_pre = nanmean(ice_dur3(:,:,n_pre),3);
    ice_dur_fut1 = nanmean(ice_dur3(:,:,n1),3);
    ice_dur_fut2 = nanmean(ice_dur3(:,:,n2),3);
    ice_dur_fut3 = nanmean(ice_dur3(:,:,n3),3);   
    
    ice_dur_std_pre = nanmean(ice_dur3_std(:,:,n_pre),3);
    ice_dur_std_fut1 = nanmean(ice_dur3_std(:,:,n1),3);
    ice_dur_std_fut2 = nanmean(ice_dur3_std(:,:,n2),3);
    ice_dur_std_fut3 = nanmean(ice_dur3_std(:,:,n3),3);
    
    %
    ice_dur_fut1(ice_dur_fut1 > ice_dur_pre) = nan;
    ice_dur_fut2(ice_dur_fut2 > ice_dur_pre) = nan;
    ice_dur_fut3(ice_dur_fut3 > ice_dur_pre) = nan;
    
    ice_dur_std_fut1(ice_dur_fut1 > ice_dur_pre) = nan;
    ice_dur_std_fut2(ice_dur_fut2 > ice_dur_pre) = nan;
    ice_dur_std_fut3(ice_dur_fut3 > ice_dur_pre) = nan;
    
    % figure 1 for paper
    close all
    fh = figure(5);
    set(fh,'color','white','Units', 'Inches', 'Position', [0,0,7.2,4],...
        'PaperUnits', 'Inches', 'PaperSize', [7.2,4]);
        ax1 = subplot(121);
            m_proj('stereographic','lat',90,'long',30,'radius',60);
            m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
            m_grid('linewi',0.01,...
                'fontname','Arial','fontsize',7,'linest',':','col','k',...
                'xtick',[],'ytick',[],'linestyle','-',...
                'box','on','backcolor','w');
            hold on;
            m_pcolor(ilon2,ilat2,ice_dur_pre');
            c1 = colorbar('Location','southoutside');
            caxis(ax1,[0 360]);
            colormap(ax1,cmocean('-ice',18));
            ylabel(c1,'Number of days with safe ice (1900-1929)',...
                'fontname','Arial','fontsize',7);        
            set(ax1,'position',[0.02    0.14    0.45    0.85]);
            set(c1,'position',[0.08   0.1    0.34    0.04],...
                'fontname','Arial','fontsize',7,...
                'YTick',0:80:320);
    
        ax2 = subplot(122);
            plot(rank_itas,ice_dur_ts_anom2,'color',[.8 .8 .8],'linewidth',0.7)
            hold on;
            plot(rank_itas,av_ice_dur_ts_anom2,'color','k','linewidth',1.5)
            xlim(ax2,[-0.25 3.5]);
            yy2 = ylabel(ax2,'Change in the duration of safe ice per year (days)',...
                'fontname','Arial','fontsize',7);
            xx2 = xlabel(ax2,['Change in global air temperature (' char(176) 'C)'],...
                'fontname','Arial','fontsize',7);
            set(ax2,'position',[0.55    0.14    0.4    0.8150],...
                'XMinorTick','on','YMinorTick','on',...
                'fontname','Arial','fontsize',7);
            hold on;
            plot([1.5 1.5],[-35 5],':k');
            hold on;
            plot([2 2],[-35 5],':k');
            hold on;
            plot([3 3],[-35 5],':k');

            h = zeros(2,1);
            ccols = {[.8 .8 .8],[0 0 0]};
            for ii = 1:2
                h(ii) = plot(nan,nan,'-','color',ccols{ii},...
                    'linewidth',1.5,'visible','on');
            end
            ll = legend(h,{'Model ensemble','Average'},...
                'Location','northeast','FontName','Arial',...
                'fontsize',7,'NumColumns',1,'AutoUpdate','off');
            legend('boxoff');   
            
        xx2.Position(2) = -38;
        yy2.Position(1) = -0.6;
        
        annotation(fh,'textbox',[0.03 0.96 0.03 0.04],...
            'String',{'a'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.47 0.96 0.03 0.04],...
            'String',{'b'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
                           
    % absolute difference between pre industrial and future
    ice_dur_diff1 = ice_dur_fut1 - ice_dur_pre;
    ice_dur_diff2 = ice_dur_fut2 - ice_dur_pre;
    ice_dur_diff3 = ice_dur_fut3 - ice_dur_pre;   

    % percent difference between pre industrial and future
    ice_dur_diff1p = -(100.*((ice_dur_pre - ice_dur_fut1)./ice_dur_pre));
    ice_dur_diff2p = -(100.*((ice_dur_pre - ice_dur_fut2)./ice_dur_pre));
    ice_dur_diff3p = -(100.*((ice_dur_pre - ice_dur_fut3)./ice_dur_pre));    
           
    % load latitudinal ensemble
    myncfile = [fid '/' 'Iceduration_4inch_zonalmean_NH.nc'];
    ice_lat = ncread(myncfile,'iceduration');
    iilat = ncread(myncfile,'lat');
    idx_lat2 = find(iilat >= 30);
    
    ice_lat_pre = squeeze(nanmean(ice_lat(idx_lat2,n_pre,:),2));
    ice_lat_fut1 = squeeze(nanmean(ice_lat(idx_lat2,n1,:),2));
    ice_lat_fut2 = squeeze(nanmean(ice_lat(idx_lat2,n2,:),2));
    ice_lat_fut3 = squeeze(nanmean(ice_lat(idx_lat2,n3,:),2)); 
    
    %
    ice_lat_diff1 = ice_lat_fut1 - ice_lat_pre;
    ice_lat_diff2 = ice_lat_fut2 - ice_lat_pre;
    ice_lat_diff3 = ice_lat_fut3 - ice_lat_pre;

    % percent difference between pre industrial and future
    ice_lat_diff1p = -(100.*((ice_lat_pre - ice_lat_fut1)./ice_lat_pre));
    ice_lat_diff2p = -(100.*((ice_lat_pre - ice_lat_fut2)./ice_lat_pre));
    ice_lat_diff3p = -(100.*((ice_lat_pre - ice_lat_fut3)./ice_lat_pre));

    %
    ice_lat_diff1_mean = nanmean(ice_lat_diff1,2);
    ice_lat_diff2_mean = nanmean(ice_lat_diff2,2);
    ice_lat_diff3_mean = nanmean(ice_lat_diff3,2);
    
    ice_lat_diff1_sd = nanstd(ice_lat_diff1,[],2);
    ice_lat_diff2_sd = nanstd(ice_lat_diff2,[],2);
    ice_lat_diff3_sd = nanstd(ice_lat_diff3,[],2);
    
    %
    ice_lat_diff1p_mean = nanmean(ice_lat_diff1p,2);
    ice_lat_diff2p_mean = nanmean(ice_lat_diff2p,2);
    ice_lat_diff3p_mean = nanmean(ice_lat_diff3p,2);
    
    ice_lat_diff1p_sd = nanstd(ice_lat_diff1p,[],2);
    ice_lat_diff2p_sd = nanstd(ice_lat_diff2p,[],2);
    ice_lat_diff3p_sd = nanstd(ice_lat_diff3p,[],2);
    
    % figure 2 for paper
    close all
    fh = figure(8);
    set(fh,'color','white','Units', 'Inches', 'Position', [0,0,7.2,6],...
        'PaperUnits', 'Inches', 'PaperSize', [7.2,6]);
        
        ax1 = subplot(231);
            m_proj('stereographic','lat',90,'long',30,'radius',60);
            m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
            m_grid('linewi',0.01,...
                'fontname','Arial','fontsize',7,'linest',':','col','k',...
                'xtick',[],'ytick',[],'linestyle','-',...
                'box','on','backcolor','w');
            hold on;
            m_pcolor(ilon2,ilat2,ice_dur_diff1');
            caxis(ax1,[-50 10]);
            colormap(ax1,cmocean('-balance','pivot',0,12)); 
            title(ax1,['1.5' char(176) 'C global warming'],'fontname','Arial','fontsize',7);

        ax2 = subplot(232);
            m_proj('stereographic','lat',90,'long',30,'radius',60);
            m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
            m_grid('linewi',0.01,...
                'fontname','Arial','fontsize',7,'linest',':','col','k',...
                'xtick',[],'ytick',[],'linestyle','-',...
                'box','on','backcolor','w');
            hold on;
            m_pcolor(ilon2,ilat2,ice_dur_diff3');
            caxis(ax2,[-50 10]);
            colormap(ax2,cmocean('-balance','pivot',0,12));
            c1 = colorbar('Location','southoutside',...
                'fontname','Arial','fontsize',7);
            ylabel(c1,'Change in the number of days with safe ice per year (days)',...
                'fontname','Arial','fontsize',7);
            title(ax2,['3' char(176) 'C global warming'],'fontname','Arial','fontsize',7);
           
        ax3 = subplot(233);
            plot(ice_lat_diff1_mean,ilat2,'-','color',rcp_cols{1},...
                'linewidth',1.2);
            hold on;
            plot(ice_lat_diff1_mean-ice_lat_diff1_sd,ilat2,'-',...
                'color',[rcp_cols{1} .5],'linewidth',1.2);
            hold on;
            plot(ice_lat_diff1_mean+ice_lat_diff1_sd,ilat2,'-',...
                'color',[rcp_cols{1} .5],'linewidth',1.2);
            
            %
            plot(ice_lat_diff3_mean,ilat2,'-','color',rcp_cols{3},...
                'linewidth',1.2);
            hold on;
            plot(ice_lat_diff3_mean-ice_lat_diff3_sd,ilat2,'-',...
                'color',[rcp_cols{3} .5],'linewidth',1.2);
            hold on;
            plot(ice_lat_diff3_mean+ice_lat_diff3_sd,ilat2,'-',...
                'color',[rcp_cols{3} .5],'linewidth',1.2);
            ylim(ax3,[28 82]);
            yy1 = ylabel(ax3,['Latitude (' char(176) ')'],...
                'fontname','Arial','fontsize',7);   
            xx1 = xlabel(ax3,{'Change in the duration';'of safe ice (days)'},...
                'fontname','Arial','fontsize',7);
            box off;
                    
        ax4 = subplot(234);
            m_proj('stereographic','lat',90,'long',30,'radius',60);
            m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
            m_grid('linewi',0.01,...
                'fontname','Arial','fontsize',7,'linest',':','col','k',...
                'xtick',[],'ytick',[],'linestyle','-',...
                'box','on','backcolor','w');
            hold on;
            m_pcolor(ilon2,ilat2,ice_dur_diff1p');
            caxis(ax4,[-100 10]);
            colormap(ax4,cmocean('-balance','pivot',0,11)); 
            title(ax4,['1.5' char(176) 'C global warming'],...
                'fontname','Arial','fontsize',7);

        ax5 = subplot(235);
            m_proj('stereographic','lat',90,'long',30,'radius',60);
            m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
            m_grid('linewi',0.01,...
                'fontname','Arial','fontsize',7,'linest',':','col','k',...
                'xtick',[],'ytick',[],'linestyle','-',...
                'box','on','backcolor','w');
            hold on;
            m_pcolor(ilon2,ilat2,ice_dur_diff3p');
            caxis(ax5,[-100 10]);
            colormap(ax5,cmocean('-balance','pivot',0,11));
            c2 = colorbar('Location','southoutside',...
                'fontname','Arial','fontsize',7);
            ylabel(c2,'Percent change in the duration of safe ice per year (%)',...
                'fontname','Arial','fontsize',7);
            title(ax5,['3' char(176) 'C global warming'],...
                'fontname','Arial','fontsize',7);
            
        ax6 = subplot(236);
            plot(ice_lat_diff1p_mean,ilat2,'-','color',rcp_cols{1},...
                'linewidth',1.2);
            hold on;
            plot(ice_lat_diff1p_mean-ice_lat_diff1p_sd,ilat2,'-',...
                'color',[rcp_cols{1} .5],'linewidth',1.2);
            hold on;
            plot(ice_lat_diff1p_mean+ice_lat_diff1p_sd,ilat2,'-',...
                'color',[rcp_cols{1} .5],'linewidth',1.2);
           
            %
            plot(ice_lat_diff3p_mean,ilat2,'-','color',rcp_cols{3},...
                'linewidth',1.2);
            hold on;
            plot(ice_lat_diff3p_mean-ice_lat_diff3p_sd,ilat2,'-',...
                'color',[rcp_cols{3} .5],'linewidth',1.2);
            hold on;
            plot(ice_lat_diff3p_mean+ice_lat_diff3p_sd,ilat2,'-',...
                'color',[rcp_cols{3} .5],'linewidth',1.2);
            ylim(ax6,[28 82]);
            yy2 = ylabel(ax6,['Latitude (' char(176) ')'],...
                'fontname','Arial','fontsize',7);   
            xx2 = xlabel(ax6,{'Percent change in the duration';'of safe ice (%)'},...
                'fontname','Arial','fontsize',7);
            box off;
            
            %
            h3 = zeros(2,1);
            rcp_cols2 = {rcp_cols{1},rcp_cols{3}};
            for ii = 1:2
                h3(ii) = plot(nan,nan,'o','color',rcp_cols2{ii},...
                    'MarkerFaceColor',rcp_cols2{ii},'markers',3,...
                    'visible','on');
            end
            ll3 = legend(h3,{['1.5' char(176) 'C warming'],...
                ['3' char(176) 'C warming']},...
                'Location','northwest','FontName','Arial',...
                'fontsize',7,'NumColumns',1,'AutoUpdate','off');
            legend('boxoff');
            
            set(ax1,'position',[0.02    0.58    0.3    0.4],...
                'fontname','Arial','fontsize',7);
            set(ax2,'position',[0.34    0.58    0.3    0.4],...
                'fontname','Arial','fontsize',7);
            set(ax3,'position',[0.71    0.6    0.27    0.36],...
                'fontname','Arial','fontsize',7,'XMinorTick','on','YMinorTick','on');
            
            set(ax4,'position',[0.02    0.08    0.3    0.4],...
                'fontname','Arial','fontsize',7);
            set(ax5,'position',[0.34    0.08    0.3    0.4],...
                'fontname','Arial','fontsize',7);
            set(ax6,'position',[0.71    0.1    0.27    0.36],...
                'fontname','Arial','fontsize',7,...
                'XTick',-100:20:20,'XMinorTick','on','YMinorTick','on');
              
            set(c1,'position',[0.15    0.56    0.35    0.0209],...
                'fontsize',7,'fontname','Arial');
            set(c2,'position',[0.15    0.06    0.35    0.0208],...
                'fontsize',7,'fontname','Arial');
            
            set(ll3,'position',[0.7    0.3777    0.1680    0.0671],...
                'fontname','Arial','fontsize',7);
        
            yy1.Position(1) = -56;
            yy2.Position(1) = -91;
            
            xx1.Position(2) = 22.5;
            xx2.Position(2) = 22.5;
            
        annotation(fh,'textbox',[0.02 0.95 0.03 0.04],...
            'String',{'a'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.33 0.95 0.03 0.04],...
            'String',{'b'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.64 0.95 0.03 0.04],...
            'String',{'c'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.02 0.45 0.03 0.04],...
            'String',{'d'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.33 0.45 0.03 0.04],...
            'String',{'e'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.64 0.45 0.03 0.04],...
            'String',{'f'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        
        %
        ikopp3 = ikopp2b;
        ikopp3(isnan(ice_dur_diff1)) = nan;

        %
        num_kopp = nan(1,5);
        for ii = 1:length(num_kopp)
            ikopp3b = ikopp3;
            ikopp3b(ikopp3b ~= ii) = nan;
            num_kopp(ii) = sum(~isnan(ikopp3b(:)));
        end

        % fig 3 for paper       
        close all
        fh = figure(8);
        set(fh,'color','white','Units', 'Inches', 'Position', [0,0,7.2,5.5],...
            'PaperUnits', 'Inches', 'PaperSize', [7.2,5.5]);

            ax1 = subplot(2,3,1:2);
                m_proj('robinson','lat',[-90 90]);
                m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
                m_grid('linewi',0.01,'fontname','Arial','fontsize',7,...
                    'linest',':','col','k','linestyle','-',...
                    'box','on','backcolor','w');
                hold on;
                m_pcolor(ilon2,ilat,ikopp2');
                colormap(ax1,major_cols); 
                c1 = colorbar('Location','southoutside');
                set(c1,'YTick',[1.5 2.2 3 3.8 4.7],'YTickLabel',{'Tropical','Dry',...
                    'Mild Temperate','Snow','Polar'},...
                    'fontname','Arial','fontsize',7);

            ax2 = subplot(2,3,3);
                scatter(1:5,num_kopp,30,'filled',...
                    'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5])
                xlim(ax2,[0 6]);
                ylim(ax2,[-10 2500]);
                ylabel(ax2,'Number of lakes','fontname','Arial',...
                    'fontsize',7);
                set(ax2,'XTick',1:5,'XTickLabel',{'Tropical','Dry',...
                    'Mild Temperate','Snow','Polar'},...
                    'fontname','Arial','fontsize',7);

        %
        num_kopp2 = num_kopp(num_kopp > 0);    
        cc = [2,3,4,5];

        %
        aa1 = nan(length(ilon2)*length(ilat2),length(num_kopp2));
        aa2 = aa1;
        aa3 = aa1;
        aa4 = aa1;
        for ii = 1:length(cc)

            %
            ice_dur_diff1_kopp = ice_dur_diff1;
            ice_dur_diff1_kopp(ikopp2b ~= cc(ii)) = nan;

            %
            ice_dur_diff3_kopp = ice_dur_diff3;
            ice_dur_diff3_kopp(ikopp2b ~= cc(ii)) = nan;

            %
            ice_dur_diff1p_kopp = ice_dur_diff1p;
            ice_dur_diff1p_kopp(ikopp2b ~= cc(ii)) = nan;

            %
            ice_dur_diff3p_kopp = ice_dur_diff3p;
            ice_dur_diff3p_kopp(ikopp2b ~= cc(ii)) = nan;


            %
            aa1(:,ii) = ice_dur_diff1_kopp(:);
            aa2(:,ii) = ice_dur_diff3_kopp(:);
            aa3(:,ii) = ice_dur_diff1p_kopp(:);
            aa4(:,ii) = ice_dur_diff3p_kopp(:);
        end

        %  
        ax3 = subplot(223);
            errorbar(0.9:1:3.9,nanmean(aa1),nanstd(aa1),'o',...
                'MarkerFaceColor',CT(9,:),'color',CT(9,:))
            hold on;
            errorbar(1.1:1:4.1,nanmean(aa2),nanstd(aa2),'o',...
                'MarkerFaceColor',CT(2,:),'color',CT(2,:));
            xlim(ax3,[0.5 4.5]);
            set(ax3,'XTick',1:4,'XTickLabel',{'Dry',...
                'Mild Temperate','Snow','Polar'},...
                'fontname','Arial','fontsize',7);
            ylabel(ax3,{'Change in the duration';'of safe ice (days)'},...
                'fontname','Arial','fontsize',7);

        %
        ax4 = subplot(224);
            errorbar(0.9:1:3.9,nanmean(aa3),nanstd(aa3),'o',...
                'MarkerFaceColor',CT(9,:),'color',CT(9,:))
            hold on;
            errorbar(1.1:1:4.1,nanmean(aa4),nanstd(aa4),'o',...
                'MarkerFaceColor',CT(2,:),'color',CT(2,:))   
            xlim(ax4,[0.5 4.5]);
            set(ax4,'XTick',1:4,'XTickLabel',{'Dry',...
                'Mild Temperate','Snow','Polar'},...
                'fontname','Arial','fontsize',7);
            ylabel(ax4,{'Percent change in the duration';'of safe ice (%)'},...
                'fontname','Arial','fontsize',7);
            
        h = zeros(2,1);
        cc2 = [CT(9,:);CT(2,:)];
        for ii = 1:2
            h(ii) = plot(nan,nan,'o','color',cc2(ii,:),...
                'MarkerFaceColor',cc2(ii,:),'visible','on',...
                'markers',4);
        end
        ll = legend(h,{['1.5 ' char(176) 'C warming'],['3 ' char(176) 'C warming']},...
            'Location','northoutside','FontName','Arial',...
            'fontsize',7,'NumColumns',1,'AutoUpdate','off');        

        set(ax1,'position',[0.06    0.54    0.53    0.47]);
        set(c1,'position',[0.13    0.54    0.38    0.02]);
        set(ax2,'position',[0.66    0.6    0.3    0.32],...
            'YMinorTick','on','XMinorTick','on');
        set(ax3,'position',[0.12    0.09    0.35    0.3412],...
            'YMinorTick','on','XMinorTick','on');
        set(ax4,'position',[0.57    0.09    0.35    0.3412],...
            'YMinorTick','on','XMinorTick','on');
        set(ll,'position',[0.75    0.11    0.1622    0.0521]);
        
        annotation(fh,'textbox',[0.04 0.95 0.03 0.04],...
            'String',{'a'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.58 0.95 0.03 0.04],...
            'String',{'b'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.04 0.45 0.03 0.04],...
            'String',{'c'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.52 0.45 0.03 0.04],...
            'String',{'d'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
                     
    % relate the percent change to population
    pop_dens2b = pop_dens2;
    pop_dens3 = pop_dens2(:);    
    ice_dur_diff1p2 = ice_dur_diff1p(:);
    ice_dur_diff2p2 = ice_dur_diff2p(:);
    ice_dur_diff3p2 = ice_dur_diff3p(:);
    pop_bin_avg1 = nan(size(pop_bins,1),1);
    pop_bin_avg1b = nan(size(pop_bins,1),1);
    pop_bin_avg2 = nan(size(pop_bins,1),1);
    pop_bin_avg2b = nan(size(pop_bins,1),1);
    pop_bin_avg3 = nan(size(pop_bins,1),1);       
    pop_bin_avg3b = nan(size(pop_bins,1),1);
    for ii = 1:size(pop_bins,1)
        
        %
        pop_dens2b(pop_dens2 >= pop_bins(ii,1) & pop_dens2 < pop_bins(ii,2)) = ii;
        
        %
        idx_bin = find(pop_dens3 >= pop_bins(ii,1) & pop_dens3 < pop_bins(ii,2));
        
        %
        pop_bin_avg1(ii) = nanmean(ice_dur_diff1p2(idx_bin));
        pop_bin_avg1b(ii) = SEM_calc(ice_dur_diff1p2(idx_bin));
        
        %
        pop_bin_avg2(ii) = nanmean(ice_dur_diff2p2(idx_bin));
        pop_bin_avg2b(ii) = SEM_calc(ice_dur_diff2p2(idx_bin));
        
        %
        pop_bin_avg3(ii) = nanmean(ice_dur_diff3p2(idx_bin));
        pop_bin_avg3b(ii) = SEM_calc(ice_dur_diff3p2(idx_bin));
    end
    
    % remove grids without ice
    pop_dens2c = pop_dens2b;
    pop_dens2c(isnan(ice_dur_pre)) = nan;
    
    % figure 4
    close all
    fh = figure(4);
    set(fh,'color','white','Units', 'Inches', 'Position', [0,0,7.2,4],...
        'PaperUnits', 'Inches', 'PaperSize', [7.2,4]);
        ax1 = subplot(121);
            m_proj('stereographic','lat',90,'long',30,'radius',60);
            m_coast('patch',[.6 .6 .6],'edgecolor',[.4 .4 .4]);
            m_grid('linewi',0.01,'fontname','Arial','fontsize',7,...
                'linest',':','col','k','xtick',[],'ytick',[],'linestyle','-',...
                'box','on','backcolor','w');
            hold on;
            m_pcolor(ilon2,ilat2,pop_dens2c');
            caxis(ax1,[0 10]);
            colormap(ax1,cmocean('dense',10))
            c1 = colorbar('Location','southoutside',...
                'fontname','Arial','fontsize',7);
            set(c1,'YTick',0.5:1:9.5,'YTickLabel',...
                {'0-1','1-5','5-10','10-20','20-50','50-100',...
                '100-200','200-500','500-1000','>1000'},...
                'fontname','Arial','fontsize',7);
            ylabel(c1,'Population density (People/km^{2})',...
                'fontname','Arial','fontsize',7);
    
        ax2 = subplot(122);
            X = categorical({'0-1','1-5','5-10','10-20','20-50','50-100',...
                '100-200','200-500','500-1000','>1000'});
            X = reordercats(X,{'0-1','1-5','5-10','10-20','20-50','50-100',...
                '100-200','200-500','500-1000','>1000'});
            bb = bar(X,[-pop_bin_avg1 -pop_bin_avg2 -pop_bin_avg3]');
            box off;
            bb(1).FaceColor = ispec(9,:);
            bb(2).FaceColor = ispec(5,:);
            bb(3).FaceColor = ispec(2,:);
            bb(1).EdgeColor = 'none';
            bb(2).EdgeColor = 'none';
            bb(3).EdgeColor = 'none';
            ylabel(ax2,'Percent change in the duration of safe ice (%)',...
                'fontname','Arial','fontsize',7);
            xlabel(ax2,'Population density (People/km^{2})',...
                'fontname','Arial','fontsize',7);
            ll = legend({['1.5 ' char(176) 'C'],...
                ['2 ' char(176) 'C'],['3 ' char(176) 'C']},...
                'Location','northoutside','FontName','Arial',...
                'fontsize',7,'NumColumns',3,'AutoUpdate','off');
            legend('boxoff');

            set(ax1,'position',[-0.04    0.2    0.55    0.77]);
            set(c1,'position',[0.05   0.145    0.38    0.04]);
            
            set(ax2,'position',[0.55    0.17    0.43    0.74],...
                'fontname','Arial','fontsize',7,...
                'YTick',0:10:100,'YTickLabel',{'0','-10','-20','-30',...
                '-40','-50','-60','-70','-80','-90','-100'});
            
            set(ll,'position',[0.54    0.9427    0.4942    0.0417]);

        annotation(fh,'textbox',[0.02 0.95 0.03 0.04],...
            'String',{'a'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
        annotation(fh,'textbox',[0.45 0.95 0.03 0.04],...
            'String',{'b'},'fontsize',7,'fontname','Arial',...
            'fontweight','bold','FitBoxToText','off',...
            'EdgeColor','none');
   

% maps of viable ice season for transportation
ice_cols = [set1(3,:);ireds(7,:)];
myncfile = [fid '/' 'iceduration_grid_42inch.nc'];

% loop through ensemble
ts_area_unviable = nan(250,90);
ts_area_viable = nan(250,90);
for i = 1:90
    
    %
    ice_dur = ncread(myncfile,'iceduration',[1 1 1 i],[Inf Inf Inf 1]);
    
    %
    ice_dur2 = ice_dur(idx_lon,:,:,:);
    ice_dur2 = ice_dur2(:,idx_lat,:,:);

    % calculate pre-industrial average
    ice_dur2_pre = nanmean(ice_dur2(:,:,n_pre),3);

    % if no ice during pre-industrial, set to NA
    ice_dur3 = ice_dur2;
    for ii = 1:size(ice_dur2,3)        
        ice_dur2b = ice_dur2(:,:,ii);
        ice_dur2b(ice_dur2_pre == 0) = nan;
        ice_dur3(:,:,ii) = ice_dur2b;        
    end

    % define categories
    ice_dur2b = ice_dur3;
    ice_dur2b(ice_dur2b <= 7) = 2; % unviable
    ice_dur2b(ice_dur2b >= 7) = 1; % viable

    % total area per grid
    cell_area3 = repmat(cell_area2,1,1,size(ice_dur2b,3));

    %
    area_unviable = cell_area3;
    area_viable = cell_area3;

    %
    area_unviable(ice_dur2b ~= 2) = nan;
    area_viable(ice_dur2b ~= 1) = nan;

    %
    ts_area_unviable(:,i) = movmean(squeeze(sum(area_unviable,[1 2],'omitnan')),31);
    ts_area_viable(:,i) = movmean(squeeze(sum(area_viable,[1 2],'omitnan')),31);
end

%
ts_area_viable2a = nanmean(ts_area_viable,2)./1e6;
ts_area_unviable2a = nanmean(ts_area_unviable,2)./1e6;

%
ts_area_viable2b = nanstd(ts_area_viable,[],2)./1e6;
ts_area_unviable2b = nanstd(ts_area_unviable,[],2)./1e6;
 
%
close all
fh = figure(1);
set(fh,'color','white','Units', 'Inches', 'Position', [0,0,7.2,3.5],...
    'PaperUnits', 'Inches', 'PaperSize', [7.2,3.5]);
    
    ax1 = subplot(1,2,1);
        plot(itas,itas_arctic,'o','color',[.6 .6 .6],...
            'MarkerFaceColor',[.6 .6 .6],'markers',3)
        xlim(ax1,[-0.25 3.5]);
        box off;
        yy1 = ylabel(ax1,['Change in average air temperature above 60' char(176) 'N (' char(176) 'C)'],...
            'fontname','Arial','fontsize',7);
        xx1 = xlabel(ax1,['Change in global air temperature (' char(176) 'C)'],...
            'fontname','Arial','fontsize',7);
        set(ax1,'XMinorTick','on','YMinorTick','on',...
            'fontname','Arial','fontsize',7);
                
    ax2 = subplot(1,2,2);
        yyaxis left
        plot(rank_itas,ts_area_viable2a,'-','color',ice_cols(1,:),...
            'linewidth',1.5)
        hold on; 
        plot(rank_itas,ts_area_viable2a-ts_area_viable2b,...
            '-','color',[ice_cols(1,:) .5],...
            'linewidth',1)
        hold on; 
        plot(rank_itas,ts_area_viable2a+ts_area_viable2b,...
            '-','color',[ice_cols(1,:) .5],...
            'linewidth',1)
        xlim(ax2,[-0.25 3.5]);
        hold on;
        plot([1.5 1.5],[0 5],':k');
        hold on;
        plot([2 2],[0 5],':k');
        hold on;
        plot([3 3],[0 5],':k');
        yy2 = ylabel(ax2,'Total surface area of a viable ice season (km^{2} x 10^{6})',...
            'fontname','Arial','fontsize',7);
        xx2 = xlabel(ax2,['Change in global air temperature (' char(176) 'C)'],...
            'fontname','Arial','fontsize',7);
        set(ax2,'XMinorTick','on','YMinorTick','on',...
            'fontname','Arial','fontsize',7);
        ax2.YAxis(1).Color = 'k';
        ax2.YAxis(2).Color = 'k';
        set(ax2,'XMinorTick','on','YMinorTick','on',...
            'fontname','Arial','fontsize',7);
        
        yyaxis right
        plot(rank_itas,ts_area_unviable2a,'-','color',ice_cols(2,:),...
            'linewidth',1.5)
        hold on;
        plot(rank_itas,ts_area_unviable2a-ts_area_unviable2b,...
            '-','color',[ice_cols(2,:) .6],'linewidth',1)
        hold on;
        plot(rank_itas,ts_area_unviable2a+ts_area_unviable2b,...
            '-','color',[ice_cols(2,:) .5],'linewidth',1)
        yy3 = ylabel(ax2,'Total surface area of an unviable ice season (km^{2} x 10^{6})',...
            'fontname','Arial','fontsize',7);
        
        h = zeros(2,1);
        for ii = 1:2
            h(ii) = plot(nan,nan,'-','color',ice_cols(ii,:),...
                'linewidth',1.2,'visible','on');
        end
        ll = legend(h,{'Viable','Unviable'},...
            'Location','northoutside','FontName','Arial',...
            'fontsize',7,'NumColumns',3,'AutoUpdate','off');
        legend('boxoff');
            
	set(ax1,'position',[0.06    0.1100    0.4    0.8]);
    set(ax2,'position',[0.53    0.1100    0.4    0.8]);
    
    yy1.Position(1) = -0.5;
    xx1.Position(2) = -1.6;
    
    yy2.Position(1) = -0.5;
    xx2.Position(2) = -0.3;
    
    yy3.Position(1) = 3.7;
    
    annotation(fh,'textbox',[0.0 0.96 0.03 0.04],...
        'String',{'a'},'fontsize',7,'fontname','Arial',...
        'fontweight','bold','FitBoxToText','off',...
        'EdgeColor','none');    
    annotation(fh,'textbox',[0.46 0.96 0.03 0.04],...
        'String',{'b'},'fontsize',7,'fontname','Arial',...
        'fontweight','bold','FitBoxToText','off',...
        'EdgeColor','none');

    set(ll,'position',[0.6052    0.93    0.2490    0.0476]);
