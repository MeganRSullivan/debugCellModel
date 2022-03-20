%% Test optimization
%  optimize cell model to match simple empirical function from Galbraith 
%  & Martiny 2015; using observed fields as funciton inputs
clear all; close all;

%% load data
load('Cellmodel_input_obs.mat')

% need d0 function
%addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/WeiLei_code/my_func')

spd  = 24*60^2 ;  	% seconds per day
spa  = 365*spd ;  	% seconds per year
%rho   = 1024.5       ; % seawater density [kg/m^3];
%permil    = rho*1e-3 ; % density*[1mmol/10^3 umol]; used for conversion from umol/kg to mmol/m3;

par.po4obs = obs.po4;
par.no3obs = obs.no3 ;
par.Temp = obs.temp ;
par.PARobs = obs.PAR ;
par.M3d_EZ = obs.wetpoints ;
iprod = obs.iEuphotic ;
grid.xt = obs.gridx ;
grid.yt = obs.gridy ;

P0 = inputs.P ;
N0 = inputs.N ;
T0 = inputs.T ;
Irr0 = inputs.PAR ;

M3d_EZ = obs.wetpoints;

%% load model output
%  to optimize to the C2P field produced by a previous model run

% OutPath = '/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/v8_Jan2022/';
% 
% % load xhat (load optimized parameter values)
% load([OutPath 'Tv4_PC_DOC0.25_DOP0v8_xhat.mat']);
% xhat_GM15 = xhat;
% clear xhat
% % load model output (needed if sing model PO4 as input to cell model
% load([OutPath 'Tv4_PC_DOC0.25_DOP0v8.mat']);
% model_GM15 = data;
% clear data
% 
% par.C2P_GM15 = M3d_EZ + nan  ;
% par.C2P_GM15(iprod)  = 1./(xhat_GM15.cc* model_GM15.DIP(iprod) + xhat_GM15.dd) ;
% 
% clear model_GM15 xhat_GM15

%% optimize to original GM15 function
%par.C2P_GM15 = 1./(0.0069*par.po4obs+0.0060);
% this is equivalent
par.C2P_GM15 = M3d_EZ + nan  ;
par.C2P_GM15(iprod) = 1000./(6.9.*P0*1e6+6.0);

%% Set up
% will need the function "PackPar" on your path
%addpath('/Users/megansullivan/Documents/UC Irvine/GitHub/OCIM-BGC-Cell')

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
par.opt_Q10Photo = on;
par.opt_fRibE 	 = on;
par.opt_fStorage = off;
par.opt_kST0 	 = off;
par.opt_PStor_rCutoff = off;
par.opt_PStor_scale = off;
par.opt_PLip_PCutoff = off;
par.opt_PLip_scale = off;
par.opt_alphaS = on;
par.opt_gammaDNA = on;

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
        %par.BIO.fStorage = 0;          % set to zero to turn of luxury storage
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
%x = x0;

%% optimization function
x0    = par.p0 ;
myfun = @(x) Cellneglogpost(x, par);

options = optimoptions(@fminunc                  , ...   
                       'Display','iter'          , ...
                       'MaxFunEvals',2000        , ...
                       'MaxIter',2000            , ...
                       'TolX',5e-7               , ...
                       'TolFun',5e-7             , ...
                       'FinDiffType','central'   , ...
                       'PrecondBandWidth',Inf )  ; %    , ...
					   %'SubproblemAlgorithm','factorization') ; %testing this

                       % 'Algorithm','trust-region', ...
                       %'GradObj','off'            , ...
                       %'Hessian','on'            , ...
                        %'DerivativeCheck','off'   , ...

[xsol,fval,exitflag] = fminunc(myfun,x0,options);
	fprintf('----fminunc complete----\n')
%% run Cell model using optimized values
    [f,out] = Cellneglogpost(xsol,par);
	fprintf('----neglogpost solved for final parameter values----\n')

%%
x = xsol
if (par.Cellmodel==on)
        if (par.opt_Q10Photo == on)
            iQ10Photo = par.pindx.lQ10Photo;
            fprintf('current Q10Photo       is  % 3.2e \n', exp(x(iQ10Photo)));
            xhat.Q10Photo = exp(x(iQ10Photo));
        end
		if (par.opt_fStorage == on)
            ifStorage = par.pindx.lfStorage;
            fprintf('current fStorage       is  % 3.2e \n', exp(x(ifStorage)));
            xhat.fStorage = exp(x(ifStorage));
        end
		if (par.opt_fRibE == on)
            ifRibE = par.pindx.tfRibE;
            fprintf('current fRibE          is  % 3.2e \n', 0.5*(1+tanh(x(ifRibE))) ); %fRibE = 0.5*(1+tanh(x(ifRibE)))
            xhat.fRibE = 0.5*(1+tanh(x(ifRibE)));
        end
		if (par.opt_kST0 == on)
            ikST0 = par.pindx.lkST0;
            fprintf('current kST0           is  % 3.2e \n', exp(x(ikST0)));
            xhat.kST0 = exp(x(ikST0));
        end
		if (par.opt_PLip_PCutoff == on)
            iPLip_PCutoff = par.pindx.lPLip_PCutoff;
            fprintf('current PLip_PCutoff   is  % 3.2e \n', exp(x(iPLip_PCutoff)));
            xhat.PLip_PCutoff = exp(x(iPLip_PCutoff));
        end
		if (par.opt_PLip_scale == on)
            iPLip_scale = par.pindx.lPLip_scale;
            fprintf('current PLip_scale     is  % 3.2e \n', exp(x(iPLip_scale)));
            xhat.PLip_scale = exp(x(iPLip_scale));
        end
		if (par.opt_PStor_rCutoff == on)
            iPStor_rCutoff = par.pindx.lPStor_rCutoff;
            fprintf('current PStor_rCutoff  is  % 3.2e \n', exp(x(iPStor_rCutoff)));
            xhat.PStor_rCutoff = exp(x(iPStor_rCutoff));
        end
		if (par.opt_PStor_scale == on)
            iPStor_scale = par.pindx.lPStor_scale;
            fprintf('current PStor_scale    is  % 3.2e \n', exp(x(iPStor_scale)));
            xhat.PStor_scale = exp(x(iPStor_scale));
        end
		if (par.opt_alphaS == on)
            ialphaS = par.pindx.lalphaS;
            fprintf('current alphaS         is  % 3.2e \n', exp(x(ialphaS)));
            xhat.alphaS = exp(x(ialphaS));
        end
        if (par.opt_gammaDNA == on)
            igammaDNA = par.pindx.tgammaDNA;
            fprintf('current gammaDNA       is  % 3.2e \n', 0.5*(1+tanh(x(igammaDNA))) ); %fRibE = 0.5*(1+tanh(x(ifRibE)))
            xhat.gammaDNA = 0.5*(1+tanh(x(igammaDNA)));
        end
end

%%
CellOut = out;
%save('Celloptim2GM15orig_0316.mat','CellOut','xhat','x')

%% make the negative log posterior function
function [f CellOut] = Cellneglogpost(x,par)  
    f=0;
    M3d = par.M3d_EZ;
    iprod = find(M3d(:,:,1:2)); %production in top two layers
    P0 = par.po4obs(iprod)./10^6;   %  convert [ mmol/m^3 --> mol/L]
    N0 = par.no3obs(iprod)./10^6;   %  convert [ mmol/m^3 --> mol/L]
    T0 = par.Temp(iprod);            % [units: degC]
    Irr0 = par.PARobs(iprod);        % [units: umol photon m^-2 s^-1]
    
    %set negative phosphate values to smallest positive concentration.
	    % (negative values break CellCNP code)
	    %fprintf('replacing %d negative Phosphate concentrations with the minimum positive concentration \n',length(P0(P0<0)))
		%negDIPindx = (P0<0);
	    %P0(P0<0)= real(min(P0(P0>=0)));

		%fprintf('replacing %d negative DIN concentrations with the minimum positive concentration \n',length(N0(N0<0)))
		%negDINindx = (N0<0);
	    N0(N0<0)= real(min(N0(N0>=0)));


    [CellOut, parBIO] = CellCNPcopy(par,x, P0,N0,T0,Irr0);
    par.BIO = parBIO;
    
    eCP = CellOut.CP - par.C2P_GM15(iprod);
    WCP = d0(ones([length(CellOut.CP),1])/sum(ones([length(CellOut.CP),1])));

    f  = f + 0.5*(eCP.'*WCP*eCP);

end

