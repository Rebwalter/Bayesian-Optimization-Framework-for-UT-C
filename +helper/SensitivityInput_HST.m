function[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,PropOpticalWall,PropOpticalGround,...
    Person,ParThermalRoof,ParThermalGround,ParThermalWall]=SensitivityInput_HST(MeteoData,Sensitivity,ittm,ParVegGround,ParVegRoof)

%% GEOMETRY OF URBAN AREA  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Height_canyon = Sensitivity(ittm).Height_canyon;
Width_canyon  = Sensitivity(ittm).Width_canyon;

Radius_tree   = Sensitivity(ittm).Radius_tree;
Distance_tree = 0.1+Radius_tree;
Height_tree   = 7.5; 
Width_roof	  =	Sensitivity(ittm).Width_roof;				% Roof width (m), calculated from the land cover fraction and the street width.

fveg_G	        =	Sensitivity(ittm).fveg_G;                   % Vegetated ground raction (-)
fveg_R	        =	0;                   % Vegetated roof raction (-)
albedo_W        =   Sensitivity(ittm).albedo_w; 

%% changed RImp!
% albedo_Rimp     =   Sensitivity(ittm).albedo_r; 



Hcan_max	=	NaN;	% Maximum height of roughness elements (buidlings), (m)
Hcan_std	=	NaN;	% Standard deviation of roughness elements (buildings), (m)

trees   =	1;		% Easy switch to include (=1) and exclude (=0) trees in the urban canyon
ftree	=	1;		% DO NOT CHANGE: Tree fraction along canyon axis

if isnan(Radius_tree)==1	% Tree radius cannot be NaN
	Radius_tree	=	0;
end


% Albedo wall 
emissivity_W	=	0.95;	% Wall emissivity (-)

PropOpticalWall	=	struct('albedo',albedo_W,'emissivity',emissivity_W);


hcanyon			=	Height_canyon/Width_canyon;		% normalized canyon height(-)
wcanyon			=	Width_canyon/Width_canyon;		% normalized canyon width (-)
wroof			=	Width_roof/Width_canyon;		% normalized roof width (-)
htree			=	Height_tree/Width_canyon;		% normalized tree height (-)
radius_tree		=	Radius_tree/Width_canyon;		% normalized tree radius (-)
distance_tree	=	Distance_tree/Width_canyon;		% normalized tree-to-wall distance (-)
ratio			=	hcanyon/wcanyon;				% height to width ratio (-)

wcanyon_norm	=	wcanyon/(wcanyon+wroof);		% normalized canyon width overall (-)
wroof_norm		=	wroof/(wcanyon+wroof);			% normalized roof width overall (-)

Gemeotry_m	=	struct('Height_canyon',Height_canyon,'Width_canyon',Width_canyon,...
				'Width_roof',Width_roof,'Height_tree',Height_tree,...
				'Radius_tree',Radius_tree,'Distance_tree',Distance_tree,...
				'Hcan_max',Hcan_max,'Hcan_std',Hcan_std);

ParTree		=	struct('trees',trees,'ftree',ftree);

geometry	=	struct('hcanyon',hcanyon,'wcanyon',wcanyon,'wroof',wroof,...
				'htree',htree,'radius_tree',radius_tree,'distance_tree',distance_tree,...
				'ratio',ratio,'wcanyon_norm',wcanyon_norm,'wroof_norm',wroof_norm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SURFACE FRACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fimp_R	=	1-fveg_R;	% Impervious roof fraction (-)

Per_runoff_R	=	1;		% Percentage of excess water that leaves the system as runoff, needs to be between 0-1 [-]

FractionsRoof	=	struct('fveg',fveg_R,'fimp',fimp_R,'Per_runoff',Per_runoff_R);

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fbare_G	=	0.0;	% Bare ground raction (-)
fimp_G	=	1-fbare_G-fveg_G;	% Impervious ground raction (-)

Per_runoff_G	=	0.9;	% Percentage of excess water that leaves the system as runoff, needs to be between 0-1 [-]

FractionsGround	=	struct('fveg',fveg_G,'fbare',fbare_G,'fimp',fimp_G,'Per_runoff',Per_runoff_G);

	

% Window parameters -------------------------------------------------------
% Fraction of windows, current scheme does not build for very high window ratios (e.g. do not simulate a glass tower)
% WindowsOn    = 1;    % Include windows in the simulation (1 = yes, 0 = no)
% GlazingRatio = 0.15;  % Glazing to wall ratio (window-to-wall ratio), e.g. 0.2 is 20% windows
% 
% % Heat conduction and capacity of windows
% Uvalue      = 4.95;     % W/m^2K, U-value = Thermal conductivity / thickness (e.g. Bueno et al. 2012)
% lan_windows = NaN;      % Thermal conductivity dry solid [W/m K]
% cv_glass    = 2.1*10^6; % Volumetric heat capacity of solid glass [J/m^3 K], not used at the moment as backcombuted from the u-value
% dztot       = 0.02;     % Total thickness of all the glass layers in a windows [m]
% 
% % Parameters calculating radiation transmission through windows, albedo of
% % windows and absorption of radiatin by the window material
% SHGC                = 0.8;                      % Solar heat gain coefficient, fraction of radiation participating in indoor energy balance, based on literature
% SolarTransmittance  = 0.75.*SHGC;    % Solar radiation transmittance trhough windows, simplified according to Bueno et al. 2012 (0.75 times the solar heat gain coefficient)
% SolarAbsorptivity   = 0;                        % Fraction of solar radiation absorbed by the window material 
% SolarAlbedo         = 1-SolarTransmittance; % Albedo of window is calculated as 1 = Albedo + Transmittance + Absorption


% 
% % PERSON for MRT calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PositionPx		=	Gemeotry_m.Width_canyon/2;	% [m] position within canyon
% PositionPz		=	1.1;						% [m] height of centre of person, usually choose 1.1 m
% 
% % PersonWidth and PersonHeight are not used at the moment
% PersonWidth		=	0.06/2;		% [-] horizontal radius of ellipse describing person (=hip width / 2)
% PersonHeight	=	0.22/2;		% [-] Vertical radius of ellipse describing person (= height / 2)
% 
% % Automatic wind speed calculation at user speficied height
% HeightWind		=	1.1;		% [m] height for wind speed to calculate OTC
% 
% Person		=	struct('PositionPx',PositionPx,'PositionPz',PositionPz,...
% 	'PersonWidth',PersonWidth,'PersonHeight',PersonHeight,'HeightWind',HeightWind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTICAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ROOF OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aveg_R		=	0.2;	% Roof vegetation surface albedo (-)
% aimp_R		=	albedo_Rimp;	% Roof impervious albedo (-)
% albedo_R	=	fveg_R*aveg_R+fimp_R*aimp_R;	% equivalent roof surface albedo (-)
% 
% eveg_R		=	1 - exp(-(ParVegRoof.LAI+ParVegRoof.SAI));	% Roof vegetation emissivity (-) 
% eimp_R		=	0.95;	% Roof impervious emissivity (-)
% emissivity_R=	fveg_R*eveg_R+fimp_R*eimp_R;	% equivalent roof surface emissivity (-)
% 
% PropOpticalRoof	=	struct('aveg',aveg_R,'aimp',aimp_R,'albedo',albedo_R,...
% 					'eveg',eveg_R,'eimp',eimp_R,'emissivity',emissivity_R);

% GROUND OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aveg_G		=	0.2;	% Ground vegetation surface albedo (-)
abare_G		=	0.15;	% Ground vegetation surface albedo (-)
aimp_G		=	0.1;	% Ground impervious albedo (-)
albedo_G	=	fveg_G*aveg_G + fbare_G*abare_G + fimp_G*aimp_G;	% equivalent ground surface albedo (-)

eveg_G		=	1 - exp(-(ParVegGround.LAI+ParVegGround.SAI));	% Ground vegetation emissivity (-) 
ebare_G		=	0.95;	% Ground vegetation surface emissivity (-)
eimp_G		=	0.95;	% Ground impervious emissivity (-)
emissivity_G=	fveg_G*eveg_G + fbare_G*ebare_G + fimp_G*eimp_G;	% equivalent ground surface emissivity (-)

PropOpticalGround	=	struct('aveg',aveg_G,'abare',abare_G,'aimp',aimp_G,'albedo',albedo_G,...
						'eveg',eveg_G,'ebare',ebare_G,'eimp',eimp_G,'emissivity',emissivity_G);

% WALL OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% albedo_W		=	0.4;	% Wall surface albedo (-)

% PERSON for MRT calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PositionPx		=	Gemeotry_m.Width_canyon/2;	% [m] position within canyon
PositionPz		=	1.1;						% [m] height of centre of person, usually choose 1.1 m

% PersonWidth and PersonHeight are not used at the moment
PersonWidth		=	0.06/2;		% [-] horizontal radius of ellipse describing person (=hip width / 2)
PersonHeight	=	0.22/2;		% [-] Vertical radius of ellipse describing person (= height / 2)

% Automatic wind speed calculation at user speficied height
HeightWind		=	1.1;		% [m] height for wind speed to calculate OTC

Person		=	struct('PositionPx',PositionPx,'PositionPz',PositionPz,...
	'PersonWidth',PersonWidth,'PersonHeight',PersonHeight,'HeightWind',HeightWind);


% % Parameters specifying the volumne of the internal mass ------------------
% IntMassOn     = 1;    % 1 = Include building internal mass in the calculation, 0 = no internal mass
% FloorHeight   = 3;    % Average floor height in building (m), which is used to calcualte total floor mass in building for heat storage
% dzFloor       = Sensitivity(ittm).dzFloor;  % Average thickness of floors in building (m), which is used to calculate total floor mass in building for heat storage
% dzWall        = Sensitivity(ittm).dzWall;  % Average thickness of walls in building (m), which is used to calculate total wall mass in building for heat storage
% 
% % Thermal properties of internal building surfaces
% lan_ground_floor    = Sensitivity(ittm).lanGround;		% Building ground: Thermal conductivity dry solid [W/m K]
% cv_ground_floor     = Sensitivity(ittm).cv_Ground;   % Building ground: Volumetric heat capacity solid [J/m^3 K]
% lan_floor_IntMass   = Sensitivity(ittm).lanIntmass;	    % Internal mass within building, Floor: Thermal conductivity dry solid [W/m K]
% cv_floor_IntMass    = Sensitivity(ittm).cv_intMass;	% Internal mass within building, Floor: Volumetric heat capacity solid [J/m^3 K]
% lan_wall_IntMass    = Sensitivity(ittm).lanIntmass;	    % Internal mass within building, Walls:Thermal conductivity dry solid [W/m K]
% cv_wall_IntMass     = Sensitivity(ittm).cv_intMass;	% Internal mass within building, Walls:Volumetric heat capacity solid [J/m^3 K]
% 
% ParThermalBulidingInt		=	struct('IntMassOn',IntMassOn,'FloorHeight',FloorHeight,...
%     'dzFloor',dzFloor,'dzWall',dzWall,'lan_ground_floor',lan_ground_floor,'cv_ground_floor',cv_ground_floor,...
%     'lan_floor_IntMass',lan_floor_IntMass,'cv_floor_IntMass',cv_floor_IntMass,...
%     'lan_wall_IntMass',lan_wall_IntMass,'cv_wall_IntMass',cv_wall_IntMass);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THERMAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_imp_R	=	0.84;		% Thermal conductivity dry solid [W/m K]
cv_s_imp_R		=	1.77*10^6;	% Volumetric heat capacity solid [J/m^3 K]

ParThermalRoof	=	struct('lan_dry_imp',lan_dry_imp_R,'cv_s_imp',cv_s_imp_R);

% GROUND THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_imp_G	=	 0.75;		% Thermal conductivity dry solid [W/m K]
cv_s_imp_G		=	1.94*10^6;	% Volumetric heat capacity solid [J/m^3 K]

ParThermalGround	=	struct('lan_dry_imp',lan_dry_imp_G,'cv_s_imp',cv_s_imp_G);

% WALL THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_W	=	Sensitivity(ittm).lan_dry_W;  % 0.84;			% Thermal conductivity dry solid [W/m K]
cv_s_W		=	Sensitivity(ittm).cv_s_W;  %1.54*10^6;		% Volumetric heat capacity solid [J/m^3 K]

ParThermalWall	=	struct('lan_dry',lan_dry_W,'cv_s',cv_s_W);



% d_leaf_T	=	4;		% Leaf dimension of ground vegetation [cm]
% ParVegTree.d_leaf = d_leaf_T; 
end