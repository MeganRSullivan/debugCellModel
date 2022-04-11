%% Cell Model Plots

outPath = '/Users/megansullivan/Documents/UC Irvine/GitHub/CellTesting/defaultCellModel/';

DATEstring = datestr(now,'yyyy-mm-dd');

%% figure settings
set(groot,'defaultAxesFontName','Times',...
    'defaultAxesFontSize',14,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','on',...
    'defaultAxesYMinorTick','on');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Times',...
    'defaultTextInterpreter','latex');

%% run default cell model
%x0=[];
%[CellOut, parBIO] = CellCNPcopy(par,x0, P0,N0,T0,Irr0);
%%
lat = grid.yt;
lon = grid.xt;

% gridded output
par.CellOut.C2P = M3d*0;
par.CellOut.N2P = M3d*0;
par.CellOut.C2N = M3d*0;

par.CellOut.C2P(iprod) = CellOut.CP;
par.CellOut.N2P(iprod) = CellOut.NP;
par.CellOut.C2N(iprod) = CellOut.CN;

par.CellOut.C2P(isnan(par.CellOut.C2P)) = 0; %remove NaNs

% save cell model allocations for analysis
par.CellOut.LimType = M3d(:,:,1:2)*NaN;
par.CellOut.r       = M3d(:,:,1:2)*NaN;
par.CellOut.mu      = M3d(:,:,1:2)*NaN;
par.CellOut.E       = M3d(:,:,1:2)*NaN;
par.CellOut.L       = M3d(:,:,1:2)*NaN;
par.CellOut.A       = M3d(:,:,1:2)*NaN;
par.CellOut.PLip    = M3d(:,:,1:2)*NaN;
par.CellOut.PStor   = M3d(:,:,1:2)*NaN;
par.CellOut.QP      = M3d(:,:,1:2)*NaN;
par.CellOut.QC      = M3d(:,:,1:2)*NaN;

par.CellOut.LimType(iprod) = CellOut.LimType;
par.CellOut.r(iprod)       = CellOut.r;
par.CellOut.mu(iprod)      = CellOut.mu;
par.CellOut.E(iprod)       = CellOut.E;
par.CellOut.L(iprod)       = CellOut.L;
par.CellOut.A(iprod)       = CellOut.A;
par.CellOut.PLip(iprod)    = CellOut.PLip;
par.CellOut.PStor(iprod)   = CellOut.PStor;
par.CellOut.QP(iprod)      = CellOut.QP;
par.CellOut.QC(iprod)      = CellOut.QC;

model_cell.CellOut = par.CellOut;

%% unpack cell model 
C2P     = model_cell.CellOut.C2P;
N2P     = model_cell.CellOut.N2P;
C2N     = model_cell.CellOut.C2N;
radius  = model_cell.CellOut.r;
LimType = model_cell.CellOut.LimType;
mu      = model_cell.CellOut.mu;
L       = model_cell.CellOut.L;
E       = model_cell.CellOut.E;
A       = model_cell.CellOut.A;

% set land points to NaN
C2P((C2P==0)) = NaN;
C2N((C2N==0)) = NaN;
N2P((N2P==0)) = NaN;
radius((radius==0)) = NaN;

L(~M3d(:,:,1:2))= NaN;
E(~M3d(:,:,1:2))= NaN;
A(~M3d(:,:,1:2))= NaN;


gammaDNA = par.BIO.gammaDNA; %.016;        % DNA fraction of cell
gammaLipid = .173;       % Structural lipid fraction of cell
gammaS = gammaDNA + gammaLipid;
%M=A;
S = 2*A + gammaS;

%%
xhat = par.BIO
%dim = [0.1199 0.695 0.1 0.2];
dim = [0.152042857142856 0.535476190476201 0.220089285714286 0.297619047619048];
%parstr = 'Cell Model Parameters';
parstr = {'Cell Model Parameters'};
kk =1;
if isfield(xhat,'Q10Photo')
	kk=kk+1;
    parstr{kk,1} = ['Q10Photo=' num2str(xhat.Q10Photo)];
end
if isfield(xhat,'fStorage')
	kk=kk+1;
    parstr{kk,1} = ['fStorage=' num2str(xhat.fStorage)] ;
end
if isfield(xhat,'fRibE')
	kk=kk+1;
    parstr{kk,1} = ['fRibE=' num2str(xhat.fRibE)] ;
end
if isfield(xhat,'kST0')
	kk=kk+1;
    parstr{kk,1} = ['kST0=' num2str(xhat.kST0)] ;
end
if isfield(xhat,'PLip_PCutoff')
	kk=kk+1;
    parstr{kk,1} = ['PLip PCutoff=' num2str(xhat.PLip_PCutoff)] ;
end
if isfield(xhat,'PLip_scale')
	kk=kk+1;
    parstr{kk,1} = ['PLip scale=' num2str(xhat.PLip_scale)] ;
end
if isfield(xhat,'PStor_rCutoff')
	kk=kk+1;
    parstr{kk,1} = ['PStor rCutoff=' num2str(xhat.PStor_rCutoff)] ;
end
if isfield(xhat,'PStor_scale')
	kk=kk+1;
    parstr{kk,1} = ['PStor scale=' num2str(xhat.PStor_scale)];
end
if isfield(xhat,'alphaS')
	kk=kk+1;
    parstr{kk,1} = ['alphaS=' num2str(xhat.alphaS)]
end

%% Allocations
% 
Zlevs = [0:0.1:1];


%% plot surface map of L
figure; hold on;
%contourf(lon,lat,C2P(:,:,1)); hold on
imAlpha = ones(size(L(:,:,1)));
imAlpha(isnan(L(:,:,1))) =0;
imagesc(lon,lat,L(:,:,1),'AlphaData',imAlpha)
cb=colorbar;
caxis([Zlevs(1) Zlevs(end)]);
cmocean('matter',2*(length(Zlevs)-1));
%colormap(flipud(summer(2*length(Zlevs)-2)));

%colormap(flipud(summer(15)))
[CC,hh] = contour(lon,lat,L(:,:,1),[0:0.1:0.7],'k');
clabel(CC,hh,'FontName','Times');
title('Allocation to L: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'L [volume fraction of cell]');
axis tight
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');

figTitle = 'Lsurface';
print(gcf,[outPath 'FIG_' figTitle DATEstring '.png'],'-dpng')


%% plot surface map of E
figure; hold on;
imAlpha = ones(size(E(:,:,1)));
imAlpha(isnan(E(:,:,1))) =0;
imagesc(lon,lat,E(:,:,1),'AlphaData',imAlpha)
cb=colorbar;
caxis([Zlevs(1) Zlevs(end)]);
cmocean('matter',2*(length(Zlevs)-1));
[CC,hh] = contour(lon,lat,E(:,:,1),[0:0.1:0.6],'k');
clabel(CC,hh,'FontName','Times');
title('Allocation to E: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'E [volume fraction of cell]');
axis tight
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');

figTitle = 'Esurface';
print(gcf,[outPath 'FIG_' figTitle DATEstring '.png'],'-dpng')

%% plot surface map of S
figure; hold on;
imAlpha = ones(size(S(:,:,1)));
imAlpha(isnan(S(:,:,1))) =0;
imagesc(lon,lat,S(:,:,1),'AlphaData',imAlpha)
cb=colorbar;
caxis([Zlevs(1) Zlevs(end)]);
cmocean('matter',2*(length(Zlevs)-1));
[CC,hh] = contour(lon,lat,S(:,:,1),[0:0.1:1],'k');
clabel(CC,hh,'FontName','Times','Color',[0.8 0.8 0.8]);
title('Allocation to S: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'Structure [volume fraction of cell]');
axis tight
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');

figTitle = 'Ssurface';
print(gcf,[outPath 'FIG_' figTitle DATEstring '.png'],'-dpng')


%% LimType
% Surface Limitation Type: 0=N-Lim; 1=P-Lim; 2&3=Co-Lim
Lim_cmap = [1, 0, 0; ...
    0, 0, 1; ...
    0, 1, 1; ...
    0, 1, 1];
figure; hold on
%plt=pcolor(lon,lat,LimType(:,:,1));
%set(plt,'EdgeColor','none');

imAlpha = ones(size(LimType(:,:,1)));
imAlpha(isnan(LimType(:,:,1))) =0;
colormap(Lim_cmap)
%shifted so smallest value is 1, so can use direct mapping 
image(lon,lat,LimType(:,:,1)+1,'AlphaData',imAlpha)
%colormap(Lim_cmap);
cb=colorbar('Ticks',[1,2,3,4],'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
ylabel(cb,'Limitation Type');
axis tight
%colormap(parula(length(Llevs)))

title('Phytoplankton Nutrient Limitation Type: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
grid off

figTitle = 'LimType_surf';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%----------radius-----------------------------
%% radius
figure; hold on
imAlpha = ones(size(radius(:,:,1)));
imAlpha(isnan(radius(:,:,1))) =0;
imagesc(lon,lat,radius(:,:,1),'AlphaData',imAlpha); hold on
c=colorbar;
cmocean('matter')
%colormap(flipud(summer));
[CC,hh] = contour(lon,lat,radius(:,:,1),[0.2,1,10,100],'k');
clabel(CC,hh,'FontName','Times');
%caxis([Zlevs(1) Zlevs(end)]);
%cmocean('matter',length(Zlevs)-1);

title('Cell radius: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(c,'radius [\mu m]');
axis tight; grid off
annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');

figTitle = 'radius_surface';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')

%% C:P ratio
figure; hold on;
imAlpha = ones(size(C2P(:,:,1)));
imAlpha(isnan(C2P(:,:,1))) =0;
imagesc(lon,lat,C2P(:,:,1),'AlphaData',imAlpha)
cb=colorbar;
cmocean('matter')
%colormap(flipud(summer))
[CC,hh] = contour(lon,lat,C2P(:,:,1),[106 106],'k');
clabel(CC,hh,'FontName','Times');
title('Cell Model C:P Uptake Ratio: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');
ylabel(cb,'C:P [molC/molP]');

annotation('textbox',dim,'String',parstr,'FitBoxToText','on','EdgeColor','none');
axis tight; grid off

figTitle = 'C2Psurface';
print(gcf,[outPath 'FIG_' figTitle '.png'],'-dpng')