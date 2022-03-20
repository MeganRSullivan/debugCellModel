% testing that growth rates are equal in the cell model
% runs CellCNPcopy

%% load data
load('Cellmodel_input_obs.mat')

% need d0 function
% addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/WeiLei_code/my_func')

spd  = 24*60^2 ;  	% seconds per day
spa  = 365*spd ;  	% seconds per year

par.po4obs = obs.po4;
par.no3obs = obs.no3 ;
par.Temp = obs.temp ;
par.PARobs = obs.PAR ;
par.M3d_EZ = obs.wetpoints ;
iprod = obs.iEuphotic ;
grid.xt = obs.gridx ;
grid.yt = obs.gridy ;

lon = grid.xt;
lat = grid.yt;

P0 = inputs.P ;
N0 = inputs.N ;
T0 = inputs.T ;
Irr0 = inputs.PAR ;

M3d_EZ = obs.wetpoints;

%% Set up
% will need the function "PackPar" on your path
% addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/OCIM-BGC-Cell')

on = true; off = false;

par.optim   = on ;
par.Cmodel  = off ;
par.Omodel  = off ;
par.Simodel = off ;
par.Cellmodel = on; % cellular trait model for phyto uptake stoichiometry
par.LoadOpt = off ; % if load optimial par.
par.pscale  = 0.0 ;
par.cscale  = 0.25 ; % factor to weigh DOC in the objective function

% P model parameters
par.opt_sigma = off ;
par.opt_kP_T  = off ;
par.opt_kdP   = off ;
par.opt_bP_T  = off ;
par.opt_bP    = off ;
par.opt_alpha = off ;
par.opt_beta  = off ;
%Cell model parameters
par.opt_Q10Photo = off;
par.opt_fRibE 	 = off;
par.opt_fStorage = off;
par.opt_kST0 	 = off;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale   = off;
par.opt_PLip_PCutoff  = off;
par.opt_PLip_scale    = off;
par.opt_alphaS        = off;
par.opt_gammaDNA      = off;

clear xhat
	%
	% Cell model parameters
	if exist('xhat') & isfield(xhat,'Q10Photo') 
		par.BIO.Q10Photo = real(xhat.Q10Photo);
	else
		par.BIO.Q10Photo = 1.983;		% Q10 of photosynthesis
	end
	if exist('xhat') & isfield(xhat,'fStorage')
		par.BIO.fStorage = real(xhat.fStorage);
	else
		par.BIO.fStorage = exp(-.358);  % strength of luxury P storage [L/molC]
	end
	if exist('xhat') & isfield(xhat,'fRibE')
		par.BIO.fRibE = real(xhat.fRibE);
	else
		par.BIO.fRibE = .618;           % ribosome fraction of biosynthetic apparatus

	end
	if exist('xhat') & isfield(xhat,'kST0')
		par.BIO.kST0 = real(xhat.kST0);
	else
		par.BIO.kST0 =0.185;            % specific synthesis rate of synthetic apparatus at 25degC [1/hr]
	end
	if exist('xhat') & isfield(xhat,'PLip_PCutoff')
		par.BIO.PLip_PCutoff = real(xhat.PLip_PCutoff);
	else
		par.BIO.PLip_PCutoff = exp(-14.408);  % log of [P (mol/L)] below which more PLipids are substituted with Slipids
	end
	if exist('xhat') & isfield(xhat,'PLip_scale')
		par.BIO.PLip_scale = real(xhat.PLip_scale);
	else
		par.BIO.PLip_scale = 3.0e6;  % scale factor for logistic function controlling phospholipid quota (changed from 1.0 b.c. not using log(P). changed from 1e6 to 3e6 to make transition sharper)
	end
	if exist('xhat') & isfield(xhat,'PStor_rCutoff')
		par.BIO.PStor_rCutoff = real(xhat.PStor_rCutoff);
	else
		par.BIO.PStor_rCutoff = 2.25;  % radius [um] above which cell stores luxury phosphorus?
	end
	if exist('xhat') & isfield(xhat,'PStor_scale')
		par.BIO.PStor_scale = real(xhat.PStor_scale);
	else
		par.BIO.PStor_scale = 3.00;  % scale factor for logistic function controlling luxury phosphorus storage (changed default from 1 to 3 to give sharper transition)
	end
	if exist('xhat') & isfield(xhat,'alphaS')
		par.BIO.alphaS = real(xhat.alphaS);
	else
		par.BIO.alphaS = 0.225;          % radius at which cell is all periplasm and membrane [um]
    end

    if exist('xhat') & isfield(xhat,'gammaDNA')
		par.BIO.gammaDNA = real(xhat.gammaDNA);
	else
		par.BIO.gammaDNA = 0.016;          % DNA Fraction of cell
	end

% cell model parameters that don't change
	if (par.Cellmodel==on)
		%par.BIO.gammaDNA = .016;        % DNA fraction of cell
		par.BIO.gammaLipid = .173       % structural Lipid (non-membrane or periplasm) fraction of cell
		%par.BIO.lPCutoff = -7.252;		% log of max [P] for which Plipids will be substituted with Slipids
		%par.BIO.r0Cutoff = 2.25;		% % NEED TO REDEFINE: r0Cutoff =  rFullA; ASK GEORGE % now PStor_rCutoff
		par.BIO.DNT0 = 1e-12*3.6e2*3600;    % Diffusivity of Nitrate at 25degC [m^2/hr]
		par.BIO.DPT0 = 1e-12*3.6e2*3600;    % Diffusivity of Phosphate at 25degC [m^2/hr]
		par.BIO.Q10Diffusivity = 1.5;
		par.BIO.AMin =.05;              % minimal fraction of cell dry mass that is nutrient uptake proteins
		%par.BIO.CStor = 1.00;           %replaced by PStor_scale
		par.BIO.PhiS = .67;             % specific carbon cost of synthesis [gC/gC]
		%%% BIO parameters below should remain fixed
		par.BIO.pDry = .47;             % Dry mass fraction of the cell
		par.BIO.rho = 1e-12;            % cell density [g/um^3]
		par.BIO.fProtM = 0.25;          % protein fraction of cell membranes
		par.BIO.fProtL = .7;            % protein fraction of light harvesting apparatus
		par.BIO.PDNA = .095;            % phosphorus mass fraction in DNA [gP/g]
		par.BIO.PRib = 0.047;           % phosphorus mass fraction in ribosomes [gP/g]
		par.BIO.PPhospholipid = 0.042;  % phosphorus mass fraction in phospholipids [gP/g]
		par.BIO.NProt = .16;            % nitrogen mass fraction in proteins [gN/g]
		par.BIO.NDNA = .16;             % nitrogen mass fraction in DNA [gN/g]
		par.BIO.NRib = .16;             % nitrogen mass fraction in Ribosomes [gN/g]
		par.BIO.CProt = .53;            % carbon mass fraction in proteins [gC/g]
		par.BIO.CDNA = .36;             % carbon mass fraction in DNA [gC/g]
		par.BIO.CPhospholipid = .65;    % carbon mass fraction in phospholipids [gC/g] - why seperate form lipids?
		par.BIO.CLipid = .76;			% carbon mass fraction in other lipids (that are not phospholipids) [gC/g]
		par.BIO.CRib = .419;     		% carbon mass fraction in ribosomes [gC/g] (technically, only correct for eukaryotes)
		par.BIO.alphaPLip = 0.12;       % phospholipid fraction of cell membrane
    end

% PackPar
par = PackPar(par)
x0 = par.p0;

%% remove the density conversion (WOA data multiplied by permil in preprocessing stage)
% rho   = 1024.5       ; % seawater density [kg/m^3];
% permil    = rho*1e-3 ; % density*[1mmol/10^3 umol]; used for conversion from umol/kg to mmol/m3;
% P0 = P0/permil;
% N0 = N0/permil;

%% run cell model
x=x0;
[CellOut, parBIO] = CellCNPcopy(par,x, P0,N0,T0,Irr0);
fprintf('done  \n')
% differences from CellCNP:
% multiplies aP and aN by 1000
% change muP to include fProtAOpt

% mu, radius, C2P seem reasonable. but now the entire ocean is N-Limited.
%% 
muE = M3d_EZ*NaN;
muL = M3d_EZ*NaN;
muN = M3d_EZ*NaN;
muP = M3d_EZ*NaN;

muE(iprod) = CellOut.mu;
muL(iprod) = CellOut.muL;
muN(iprod) = CellOut.muN;
muP(iprod) = CellOut.muP;

mu = muE;

%%
fprintf('The average growth rate from biosynthesis is   %6.4f \n',mean(muE(:),'omitnan'));
fprintf('The average growth rate from photosynthesis is %6.4f \n',mean(muL(:),'omitnan'));
fprintf('The average growth rate from N uptake is       %6.4f \n',mean(muN(:),'omitnan'));
fprintf('The average growth rate from P uptake is       %6.4f \n',mean(muP(:),'omitnan'));

%%
C2P = M3d_EZ*NaN;
C2P(iprod) = CellOut.CP;

N2P = M3d_EZ*NaN;
N2P(iprod) = CellOut.NP;

LimType = M3d_EZ*NaN;
LimType(iprod) = CellOut.LimType;

radius = M3d_EZ*NaN;
radius(iprod) = CellOut.r;

L = M3d_EZ*NaN;
L(iprod) = CellOut.L;

E = M3d_EZ*NaN;
E(iprod) = CellOut.E;

A = M3d_EZ*NaN;
A(iprod) = CellOut.A;

%% quick look: C:P
figure;
imagesc(lon,lat,C2P(:,:,1)); colorbar
axis xy
title('C:P')

%% quick look: cell size
figure;
imagesc(lon,lat,radius(:,:,1)); colorbar
axis xy
title('radius [um]')
%% quick look: limitation type
% 0: Nlimitation
% 1: P limitation
% 2&3: Colimitation

min(CellOut.LimType)
max(CellOut.LimType)

Lim_cmap = [1, 0, 0; ...
    0, 0, 1; ...
    0, 1, 1; ...
    0, 1, 1];
figure; hold on
Llevs= [0,1,2,3];

imAlpha = ones(size(LimType(:,:,1)));
imAlpha(isnan(LimType(:,:,1))) =0;
imagesc(lon,lat,LimType(:,:,1),'AlphaData',imAlpha)
cb=colorbar('Ticks',[0,1,2,3],'TickLabels',{'N-Lim','P-Lim','Co-Lim','Co-Lim-alt'});
colormap(Lim_cmap);
ylabel(cb,'Limitation Type');
axis tight

title('Phytoplankton Nutrient Limitation Type: Surface','Fontsize',18);
xlabel('Longitude');
ylabel('Latitude');


%%
lat = grid.yt;
lon = grid.xt;
figure;
imagesc(lon,lat,muE(:,:,1)-muL(:,:,1)); colorbar
axis xy
title('muE - muL')
% this looks good. machine error difference
%% 
figure;
imagesc(lon,lat,muE(:,:,1)-muP(:,:,1)); colorbar
axis xy
title('muE - muP')
% muP is much larger than muE in subtropics

%% 
figure;
imagesc(lon,lat,muE(:,:,1)-muN(:,:,1)); colorbar
axis xy
title('muE - muN')

%%
figure;
imagesc(lon,lat,muP(:,:,1)-muN(:,:,1)); colorbar
axis xy
title('muP - muN')
% muP is much larger than muN, especially in subtropics