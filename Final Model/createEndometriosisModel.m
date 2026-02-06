function m1 = createEndometriosisModel(nva)
arguments




    % Steady state values
    %IGF1...................
    nva.ssIGF1   = 15.6883 ;    % nM - not sure
    nva.ssIGFBP  = 80;          % nM - not sure
    nva.ssIGFBPc = 120;         % nM - not sure
    nva.ssPAPPA  = 40/200;       % PAPPA
    nva.ssPAPPA2 = 0.00132;        % PAPPA2 
    nva.rpIGF1 = 128.7343;%(0.17/24)*15e-02/1;  %colleti etal 2020 (rat values)0.17/24, divide by 7 as at homestasis tumor factor is about 7
    nva.raIGF1 = (0.148/24)*15e-02/1; %0.148/24 % 2x growth over 6month
    % Turnover rates
    nva.trIGF1   = 0.167;       % hours
    nva.trIGFBP  = 16;          % hours
    nva.trIGFBPc = 16;          % hours
    nva.trPAPPA  = 100;         % hours
    nva.trPAPPA2 = 100;         % hours

    nva.Kcat2     = 0.31 * 3600; % kcat = (0.31 +/- 0.06) /s  --> * 3600 1/h
    nva.Km2       = 0.41 * 1e3;  % km   =  0.41 +/- 0.11  uM  --> * 1e3 nM 
    nva.Kcat    = 0.89 * 3600; % kcat = (0.89 +/- 0.08) 1/s --> * 3600 1/h
    nva.Km      = 1.63 * 1e3;  % km   =  1.63 +/- 0.32  uM  --> * 1e3 nM 

    nva.FupmAB = 1; % 5e-4; % where did you get this ???? for mAb this is usually closer to 1.0 
    %SLC71A3
    nva.ssARG   = 10e+4 ;    % nM - not sure
    nva.ssARGi   = 8e+4 ;    % nM - not sure
    nva.ssNO  = 10;          % nM - not sure
    

    % Turnover rates
    nva.trARG   = 2;       % hours
    nva.trARGi   = 1;       % hours
    nva.trNO  = 10;          % hour
    nva.rpARG = 41.6;%(0.17/24)*15e-02/1;  %colleti etal 2020 (rat values)0.17/24, divide by 7 as at homestasis tumor factor is about 7
    nva.raARG = (0.148/24)*15e-02/1; %0.148/24 % 2x growth over 6month
    nva.FupARG = 1; %
    %EP3
    nva.ssPGE2   = 1.1348 ;    % nM - not sure
    nva.ssEP3  = 1.667;      % nM - not sure

    nva.kpainEP3 = 0.2821;  %
    nva.kelpainEP3 = 1; %
    % Turnover rates
    nva.trPGE2   = 0.15;       % hours
    nva.trEP3  = 12;          % hours
    nva.kon = 1;
    nva.koff = 5.8;
    nva.exportFlag = false;
end
%IGF1
nva.kdegIGFBP   = log(2)/nva.trIGFBP ; 
nva.kdegIGFBPc  = log(2)/nva.trIGFBPc ;
nva.kdegPAPPA   = log(2)/nva.trPAPPA;
nva.kdegPAPPA2  = log(2)/nva.trPAPPA2;
nva.kdegIGF1    = log(2)/nva.trIGF1 ; 



r1 = nva.Kcat  * nva.ssIGFBP  * nva.ssPAPPA  / (nva.Km  + nva.ssIGFBP  );
r2 = nva.Kcat2 * nva.ssIGFBPc * nva.ssPAPPA2 / (nva.Km2 + nva.ssIGFBPc );

nva.ksynIGFBP  = nva.kdegIGFBP  * nva.ssIGFBP  + r1;
nva.ksynIGFBPc = nva.kdegIGFBPc * nva.ssIGFBPc + r2; 
nva.ksynIGF1   = nva.kdegIGF1   * nva.ssIGF1 - r1 - r2;
nva.ksynPAPPA  = nva.kdegPAPPA  * nva.ssPAPPA;
nva.ksynPAPPA2 = nva.kdegPAPPA2 * nva.ssPAPPA2;
    
nva
r1
r2
%-------------------------------------------------------
%SLC7A3
nva.kmCAT1 = 110 * 1e3; %110uM TO nM multiply by 1000
nva.kmCAT3 = 200 * 1e3; %110uM TO nM multiply by 1000
nva.VmaxCAT1 = 1.6 ;
nva.VmaxCAT3 = 2 ;
r1 = (nva.VmaxCAT3 * nva.ssARG /(nva.kmCAT3 + nva.ssARG ));
nva.kdegARG   = log(2)/nva.trARG ; 
nva.kdegARGi   = log(2)/nva.trARGi ;
nva.kdegNO  = log(2)/nva.trNO ;
nva.ksynARG  = nva.kdegARG  * nva.ssARG + r1 ;
nva.ksynNO = nva.kdegNO * nva.ssNO ;
%---------------------------------------------------
%EP3 Parameters
nva.kdegPGE2   = log(2)/nva.trPGE2 ; 
nva.kdegEP3  = log(2)/nva.trEP3 ;


nva.ksynPGE2  = nva.kdegPGE2  * nva.ssPGE2;
nva.ksynEP3 = nva.kdegEP3 * nva.ssEP3;
%IGF1 MODEL----------------------------------------------------
if nva.ksynIGF1 < 0
    fprintf('\n\nWARNING !!!!! NEGATIVE VALUE for ksynIGF1: %f\n\n',nva.ksynIGF1);
end

sbiounit('kcell','item',1000);

m1 = sbiomodel('PAPPA2model');
ed = m1.addcompartment('ED', 'Value', 1, 'Units', 'liter', 'constant', 1); % endometrium

% add all species
ed.addspecies('IGFBP',      'Value', nva.ssIGFBP,              'Units', 'nanomolarity');
ed.addspecies('IGFBPc',     'Value', nva.ssIGFBPc,             'Units', 'nanomolarity');
ed.addspecies('IGF1',       'Value', nva.ssIGF1,               'Units', 'nanomolarity');
ed.addspecies('PAPPA',      'Value', nva.ssPAPPA,              'Units', 'nanomolarity');
ed.addspecies('PAPPA2',     'Value', nva.ssPAPPA2,             'Units', 'nanomolarity');
ed.addspecies('totIGFBP',   'Value', nva.ssIGFBP+nva.ssIGFBPc, 'Units', 'nanomolarity');
ed.addspecies('EC',       'Value', 31000,               'Units', 'item');
ed.addspecies('EC0',       'Value', 31000,               'Units', 'item');
ed.addspecies('PainIGF1',      'Value', 0,              'Units', 'item');
ed.addspecies('PainNO',      'Value', 0,              'Units', 'item');
ed.addspecies('PainEP3',      'Value', 0,              'Units', 'item');
% Drug : Mab
ed.addspecies('mAb',        'Value', 0,                        'Units', 'nanomolarity');
ed.addspecies('mAbPAPPA2',  'Value', 0,                        'Units', 'nanomolarity');
%ed.addspecies('Pain',  'Value', 0,                        'Units', 'dimensionless');
% add all parameters
% IGFBP non cleavable by papp-a2 synthesis and degradation
m1.addparameter('ksynIGFBP', 'Value', nva.ksynIGFBP,  'Units','nanomolarity/hour');
m1.addparameter('kdegIGFBP', 'Value', nva.kdegIGFBP,  'Units','1/hour');
% IGFBPc cleavable synthesis and degradation
m1.addparameter('ksynIGFBPc','Value', nva.ksynIGFBPc, 'Units','nanomolarity/hour');
m1.addparameter('kdegIGFBPc','Value', nva.kdegIGFBPc, 'Units','1/hour');
% IGF1 synthesis and degradation
m1.addparameter('ksynIGF1',  'Value', nva.ksynIGF1,   'Units','nanomolarity/hour');
m1.addparameter('kdegIGF1',  'Value', nva.kdegIGF1,   'Units','1/hour');
% PAPPA synthesis and degradation
m1.addparameter('ksynPAPPA', 'Value', nva.ksynPAPPA,  'Units','nanomolarity/hour');
m1.addparameter('kdegPAPPA', 'Value', nva.kdegPAPPA,  'Units','1/hour');
% PAPPA2 synthesis and degradation
m1.addparameter('ksynPAPPA2','Value', nva.ksynPAPPA2, 'Units','nanomolarity/hour');
m1.addparameter('kdegPAPPA2','Value', nva.kdegPAPPA2, 'Units','1/hour');
% Michaelis-Menten Kcat and Km parameters for PAPPAx cleavage of IGFBPx
m1.addparameter('Kcat',      'Value', nva.Kcat,       'Units','1/hour');
m1.addparameter('Kcat2',     'Value', nva.Kcat2,      'Units','1/hour');
m1.addparameter('Km',        'Value', nva.Km,         'Units','nanomolarity');
m1.addparameter('Km2',       'Value', nva.Km2,        'Units','nanomolarity');
% mAb blocking PAPAA2
m1.addparameter('kon_mAb',    'Value', 1,              'Units','1/nanomolarity/hour');
m1.addparameter('koff_mAb',   'Value', 1,              'Units','1/hour');


% add reactions
% IGF1-BP non cleavable only by PAPPA synthesis and degradtion 
m1.addreaction('null -> IGFBP',   'ReactionRate', 'ksynIGFBP',        'Name', 'IGFBP synthesis');
m1.addreaction('IGFBP -> null',   'ReactionRate', 'kdegIGFBP*IGFBP',  'Name', 'IGFBP degradation');
% IGF1-BP cleavable (POI) synthesis and degradtion 
m1.addreaction('null -> IGFBPc',  'ReactionRate', 'ksynIGFBPc',       'Name', 'IGFBPc synthesis');
m1.addreaction('IGFBPc  -> null ','ReactionRate', 'kdegIGFBPc*IGFBPc','Name', 'IGFBPc degradation');
% PAPPA2 synthesis and degradtion 
m1.addreaction('null -> PAPPA2  ', 'ReactionRate','ksynPAPPA2',       'Name', 'PAPPA2 synthesis');
m1.addreaction('PAPPA2  -> null ', 'ReactionRate','kdegPAPPA2*PAPPA2','Name', 'PAPPA2 degradation');
% PAPPA synthesis and degradtion 
m1.addreaction('null -> PAPPA',   'ReactionRate', 'ksynPAPPA',        'Name', 'PAPPA synthesis');
m1.addreaction('PAPPA  -> null',  'ReactionRate', 'kdegPAPPA*PAPPA',  'Name', 'PAPPA Degradation');
% IGF1 cleavage by PAPPA
m1.addreaction('IGFBP -> IGF1',   'ReactionRate', 'Kcat *IGFBP *PAPPA /(Km +IGFBP )', 'Name', 'IGF1 cleavage by PAPPA');
m1.addreaction('IGFBPc -> IGF1',  'ReactionRate', 'Kcat2*IGFBPc*PAPPA2/(Km2+IGFBPc)', 'Name', 'IGF1 cleavage by PAPPA2');
% IGF1 sysnthesis and degradation
m1.addreaction('null -> IGF1',    'ReactionRate', 'ksynIGF1',         'Name', 'IGF1 synthesis');
m1.addreaction('IGF1  -> null',   'ReactionRate', 'kdegIGF1*IGF1',    'Name', 'IGF1 Degradation');
% mAb binding to PAPPA2
m1.addreaction('mAb + PAPPA2 <-> mAbPAPPA2 ', 'ReactionRate', 'kon_mAb*mAb*PAPPA2-koff_mAb*mAbPAPPA2', 'Name', 'mAb binding to PAPPA2');
%IGF1 MODEL--------------------------------------------------------------
%SLC7A3.........................................
% add all species
% ed.addspecies('SLC7A1',      'Value', nva.ssARG,              'Units', 'nanomolarity');
ed.addspecies('Arginine_i',      'Value', 0,              'Units', 'nanomolarity');
ed.addspecies('Arginine_e',      'Value', nva.ssARG,              'Units', 'nanomolarity');
ed.addspecies('NO',            'Value', nva.ssNO,             'Units', 'nanomolarity');
%ed.addspecies('EC',       'Value', 20000,               'Units', 'item');
%ed.addspecies('Pain',      'Value', 0,              'Units', 'item');

% Drug : small molecule
ed.addspecies('Drug',        'Value', 0,                        'Units', 'nanomolarity');

% add all parameters
% Arginine synthesis and degradation
m1.addparameter('ksynARG', 'Value', nva.ksynARG,  'Units','nanomolarity/hour');
m1.addparameter('kdegARG', 'Value', nva.kdegARG,  'Units','1/hour');
m1.addparameter('kdegARGi', 'Value', nva.kdegARGi,  'Units','1/hour')
%transporter KM and Vmax
m1.addparameter('kmCAT1', 'Value', nva.kmCAT1,  'Units','nanomolarity');
m1.addparameter('kmCAT3','Value', nva.kmCAT3, 'Units','nanomolarity');
m1.addparameter('VmaxCAT1', 'Value', nva.VmaxCAT1,  'Units','nanomolarity/hour');
m1.addparameter('VmaxCAT3','Value', nva.VmaxCAT3, 'Units','nanomolarity/hour');
% NO degradation
m1.addparameter('ksynNO', 'Value', nva.ksynNO,  'Units','nanomolarity/hour');
m1.addparameter('kdegNO','Value', nva.kdegNO, 'Units','1/hour');
% NO with Arginine uptake
m1.addparameter('EmaxNO',    'Value', 13,              'Units','dimensionless');
m1.addparameter('EC50NO',   'Value', 0.2,              'Units','nanomolarity');
% Drug
m1.addparameter('Imax',    'Value', 1,              'Units','dimensionless');
m1.addparameter('IC50',   'Value', 1,              'Units','nanomolarity');
% add reactions
% NO synthesis and degradtion 
m1.addreaction('null -> NO',   'ReactionRate', 'ksynNO*((1+ (EmaxNO * (Arginine_i/EC50NO )))/(1+  (Arginine_i/EC50NO )))',        'Name', 'NO synthesis');
m1.addreaction('NO  -> null',  'ReactionRate', 'kdegNO*NO',  'Name', 'NO Degradation');

% ARG sysnthesis and degradation
% m1.addreaction('null -> Arginine',    'ReactionRate', 'ksynARG*(1-(Imax*Drug)/(IC50 + Drug))',         'Name', 'Intracellular Arginine uptake by SLC7A3');
m1.addreaction('null -> Arginine_e',    'ReactionRate', 'ksynARG *(1-(Imax*Drug)/(IC50 + Drug))',         'Name', 'extracellular Arginine ');
m1.addreaction('Arginine_e  -> Arginine_i',   'ReactionRate', '( (VmaxCAT3 * Arginine_e /(kmCAT3 + Arginine_e )))',    'Name', 'ARGININE uptake');
m1.addreaction('Arginine_e  -> null',  'ReactionRate', 'kdegARG*Arginine_e',  'Name', 'ArgNO Degradation');
m1.addreaction('Arginine_i  -> null',  'ReactionRate', 'kdegARGi*Arginine_i',  'Name', 'ArgNOi Degradation');
% Drug binding to SLC7A3
%m1.addreaction('Drug + SLC7A3 <-> DrugSLC7A3 ', 'ReactionRate', 'kon_D*Drug*SLC713-koff_D*DrugSLC7A3', 'Name', 'Drug blocking SLC7A3');
%SLC7A3----------------------------------------------------------
%----------------------------------------------------------
%EP3 Model
%-------------------------------------------------------------
% add all species
ed.addspecies('PGE2',      'Value', nva.ssPGE2,              'Units', 'nanomolarity');
ed.addspecies('EP3',     'Value', nva.ssEP3,             'Units', 'nanomolarity');
ed.addspecies('PGE2EP3',     'Value', 0,             'Units', 'nanomolarity');
%ed.addspecies('Pain',      'Value', 0,              'Units', 'item');
% Drug : Mab
%ed.addspecies('Drug',        'Value', 0,                        'Units', 'nanomolarity');
ed.addspecies('DrugEP3',  'Value', 0,                        'Units', 'nanomolarity');
%ed.addspecies('Pain',  'Value', 0,                        'Units', 'dimensionless');
% add all parameters
% PGE2 synthesis and degradation
m1.addparameter('ksynPGE2', 'Value', nva.ksynPGE2,  'Units','nanomolarity/hour');
m1.addparameter('kdegPGE2', 'Value', nva.kdegPGE2,  'Units','1/hour');
% EP3 synthesis and Degradation
m1.addparameter('ksynEP3','Value', nva.ksynEP3, 'Units','nanomolarity/hour');
m1.addparameter('kdegEP3','Value', nva.kdegEP3, 'Units','1/hour');
% PGE2 and EP3
m1.addparameter('konEP3',    'Value', nva.kon,              'Units','1/nanomolarity/hour');
m1.addparameter('koffEP3',   'Value', nva.koff,              'Units','1/hour');

% Drug blocking EP3
m1.addparameter('kon_drug',    'Value', 1,              'Units','1/nanomolarity/hour');
m1.addparameter('koff_drug',   'Value', 1,              'Units','1/hour');


% add reactions
% PGE2 synthesis and degradtion 
m1.addreaction('null -> PGE2',   'ReactionRate', 'ksynPGE2',        'Name', 'PGE2 synthesis');
m1.addreaction('PGE2 -> null',   'ReactionRate', 'kdegPGE2*PGE2',  'Name', 'PGE2 degradation');
% EP3 synthesis and degradtion 
m1.addreaction('null -> EP3',  'ReactionRate', 'ksynEP3',       'Name', 'EP3 synthesis');
m1.addreaction('EP3  -> null ','ReactionRate', 'kdegEP3*EP3','Name', 'EP3 degradation');
%PGE2 and EP3 interaction
m1.addreaction('PGE2 + EP3 <-> PGE2EP3 ', 'ReactionRate', 'konEP3*PGE2*EP3-koffEP3*PGE2EP3', 'Name', 'PGE2 binding to EP3');
% Drug binding to EP3   
m1.addreaction('Drug + EP3 <-> DrugEP3 ', 'ReactionRate', 'kon_drug*Drug*EP3-koff_drug*DrugEP3', 'Name', 'Drug binding to EP3');
%------------------------------------------------------------------------
%-------------------------------------------------------------
%-------------------------------------------------------------
%-------------------------------------------------------------
%-------------------------------------------------------------
%-------------------------------------------------------------
%TUMOR ------------------------------------------------
%IGF1 and tumor (EC)
m1.addparameter('EmaxIGF1',    'Value', 1.2990,              'Units','dimensionless');
m1.addparameter('EC50IGF1',   'Value', 0.3535,              'Units','nanomolarity');
m1.addparameter('hcIGF1',4.3,'Units','dimensionless');
% FOR lesion growth /apoptosis IGF1 parameters
m1.addparameter('raIGF1',    'Value', nva.raIGF1,              'Units','1/hour');
m1.addparameter('rpIGF1',   'Value', nva.rpIGF1,              'Units','item/hour');

%--------------IGF1 parameters tumor--------------------


%%SLC7A3 and tumor (EC) ---------------------------
m1.addparameter('hcArg',2,'Units','dimensionless');
% tumor proliferation with Arginine uptake
m1.addparameter('EmaxARG',    'Value', 13,              'Units','dimensionless');
m1.addparameter('EC50ARG',   'Value', 0.2,              'Units','nanomolarity');
% FOR lesion growth /apoptosis parameters
m1.addparameter('raARG',    'Value', nva.raARG,              'Units','1/hour');
m1.addparameter('rpARG',   'Value', nva.rpARG,              'Units','item/hour');
%--------------SLC7A3 parameters tumor--------------------


%Final TUMOR growth with ARG and IGF1 as factors
m1.addreaction(' null  -> EC  ', 'ReactionRate', 'rpIGF1*((1+ (EmaxIGF1 * (IGF1/EC50IGF1)^hcIGF1) + (EmaxARG * (Arginine_i/EC50ARG)^hcArg))/(1+  ((IGF1/EC50IGF1)^hcIGF1) +  ((Arginine_i/EC50ARG)^hcArg)))',...
    'Name', 'tumor');
m1.addreaction('EC -> null ', 'ReactionRate', 'raIGF1*EC',...
    'Name', 'basal apoptosisrate Endo');
%-------------------------------------------------------------
%-------------------------------------------------------------
%-------------------------------------------------------------
%flux 
%-------------------------------------------------------------
%-------------------------------------------------------------

%PAIN
%%1. IGF1 and  pain%--------------------------------------------------------------------------------------------------------
m1.addparameter('EmaxIGF1pain',    'Value', 100,              'Units','dimensionless');
m1.addparameter('EC50IGF1pain',   'Value', 59.2907,              'Units','nanomolarity');
m1.addparameter('hcIGF1pain',1.3240,'Units','dimensionless');
m1.addparameter('kpainIGF1',   'Value', 0.1096,              'Units','item/hour');
m1.addparameter('kelpainIGF1', 'Value', 1, 'Units','1/hour');
%IGF1

m1.addreaction('null  -> PainIGF1  ',...
    'ReactionRate', 'kpainIGF1 *((1+ (EmaxIGF1pain * (IGF1/EC50IGF1pain )^hcIGF1pain)))/((1+  (IGF1/EC50IGF1pain )^hcIGF1pain))',...
     'Name', 'Pain pathway');
m1.addreaction('PainIGF1 -> null', 'ReactionRate', 'kelpainIGF1*PainIGF1',...
     'Name', 'pain clearance');
%-------------------------------------------------------------%-------------------------------------------------------------%-------------------------------------------------------------
%pain and IGF1-------------------------------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------
%2. SLC7A3 and pain-----------------------------------------
m1.addparameter('hcpainNO',4.5,'Units','dimensionless');
m1.addparameter('EmaxpainNO',    'Value', 100,              'Units','dimensionless');
m1.addparameter('EC50painNO',   'Value', 74,              'Units','nanomolarity');
m1.addparameter('kpainNO',   'Value', 1,              'Units','item/hour');
m1.addparameter('kelpainNO', 'Value', 1, 'Units','1/hour');

%-------------------------------------------------------------%-------------------------------------------------------------
%NO
 m1.addreaction('null -> PainNO  ', 'ReactionRate', 'kpainNO*((1+ (EmaxpainNO * (NO/EC50painNO)^hcpainNO))/(1+  (NO/EC50painNO)^hcpainNO))',...
     'Name', 'Pain synthesis');
 m1.addreaction('PainNO  -> null ', 'ReactionRate', 'kelpainNO*PainNO',...
     'Name', 'Pain Degradation');
%-------------------------------------------------------------------------------------------------------------------------------
%SLC7A3 and pain------------------------------------------------------------

%---------------------------------------------------------------------------
%3.EP3 and Pain ---------------------------
%EP3 and pain
m1.addparameter('EmaxEP3',    'Value', 200,              'Units','dimensionless');
m1.addparameter('ED50EP3',   'Value', 0.57,              'Units','nanomolarity');%kassuya paper
m1.addparameter('hcEP3',0.88,'Units','dimensionless');
m1.addparameter('kpainEP3',    'Value', nva.kpainEP3,              'Units','item/hour');
m1.addparameter('kelpainEP3',   'Value', nva.kelpainEP3,              'Units','1/hour')
%-------------------------------------------------------------%-------------------------------------------------------------#
%EP3
m1.addreaction(' null  -> PainEP3  ', 'ReactionRate', 'kpainEP3*((1+ (EmaxEP3 * (PGE2EP3/ED50EP3)^hcEP3))/(1+  (PGE2EP3/ED50EP3)^hcEP3))',...
    'Name', 'Allodynia');
m1.addreaction('PainEP3 -> null ', 'ReactionRate', 'kelpainEP3*PainEP3',...
    'Name', 'removal of pain');
%-------------------------------------------------------------%-------------------------------------------------------------
%EP3-------------------------------------------------------------
%Pain as impact of 4 aspects 1.IGF1 2.NO 3.EP3 4.Tumor--------------------

%--PAIN---------------------------------------------------
m1.addparameter('totPain', 'Value', 0, 'Units', 'item', ...
   'Constant', 0);
addrule(m1, 'totPain = (PainIGF1/10 + PainNO/10 + PainEP3/100)', 'RuleType', 'repeatedAssignment');
%-------------------------------------------------------------%-------------------------------------------------------------
m1.addparameter('Pain', 'Value', 0, 'Units', 'dimensionless', ...
   'Constant', 0);
m1.addparameter('onePain', 'Value', 1.5, 'Units', 'item', ...
   'Constant', 0);
addrule(m1, 'Pain = totPain^3/(totPain^3 + onePain^3)', 'RuleType', 'repeatedAssignment');
%--------------------------
%  m1.addreaction('null  -> Pain  ',...
%      'ReactionRate', 'min(0, max(1,kpainIGF1 *((1+ (EmaxIGF1pain * (IGF1/EC50IGF1pain )^hcIGF1pain) + (EmaxpainNO * (NO/EC50painNO)^hcpainNO) + (EmaxEP3 * (PGE2EP3/ED50EP3)^hcEP3)))/((1+  ((IGF1/EC50IGF1pain )^hcIGF1pain)+((NO/EC50painNO)^hcpainNO)+ ((PGE2EP3/ED50EP3)^hcEP3)))))',...
%       'Name', 'Pain pathway');

%  m1.addreaction('Pain -> null', 'ReactionRate', 'kelpainNO*Pain',...
%      'Name', 'pain clearance');
%--PAIN---------------------------------------------------
% Unbound drug
m1.addparameter('FupmAB', 'Value', nva.FupmAB, 'Units', 'dimensionless');
m1.addparameter('dosenM', 'Value', 0, 'Units', 'nanomolarity');
addrule(m1, 'mAb = FupmAB*dosenM', 'RuleType', 'initialAssignment');
% proliferation
% m1.addparameter('proliferation', 'Value', 1, 'Units', 'dimensionless');
% addrule(m1, 'proliferation = (EC-EC0)/EC0 * 100', 'RuleType', 'repeatedAssignment')
addrule(m1, 'totIGFBP = IGFBP + IGFBPc', 'RuleType', 'repeatedAssignment');

% Configure simulation property.
cs = getconfigset(m1, 'default');
cs.CompileOptions.UnitConversion = true;
set(cs,'StopTime',(24*2+100)); %744
set(cs,'TimeUnits','hour');
cs.SolverOptions.AbsoluteTolerance = 1e-8;
cs.SolverOptions.RelativeTolerance = 1e-8;
cs.SolverType = 'ode15s';

verify(m1);

if nva.exportFlag
    mexp = export(m1); accelerate(mexp);
end

end