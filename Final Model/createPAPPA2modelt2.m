function m1 = createPAPPA2modelt2(nva)
arguments

    % IFG1 MW : 7649 Da = g/mol = ng/nanomol
    % (PNID: 9768660) IGF1 Normal ranges for age-matched subjects are 18–274 ng/mL for IGF-I, 0.9–3.3 μg/mL for IGFBP-3, and 10–120 ng/mL for IGFBP-1

    % IFG1: 18–274 ng/mL
    %   MW_IGF1 : 7649 Da = g/mol = ng/nanomol
    %   IGF1    = 120 ng/mL 
    %           = 120 / 7649 * 1e3 nanomole/L = nanomolarity  
    %           = 15.6883 nanomolarity

    % IFGBP1: 10–120 ng/mL 
    %   MW_IGFBP1 : 25 kDa = 25e3 g/mol = ng/nanomol
    %   IGF1BP  = 54 ng/mL  
    %           = 54 / 25e3 * 1e3 nanomole/L = nanomolarity  
    %           = 2.16 nanomolarity

    % IFGBP3: 0.9–3.3 μg/mL
    %   MW_IGFBP3 : 29 kDa = 29e3 g/mol = ng/nanomol
    %   IGF1BP3 = 1856 ng/mL 
    %           = 1856 / 29e3 * 1e3 nanomole/L = nanomolarity  
    %           = 64 nanomolarity

    % IFGBP4: (189 +/- 83 microg/l and 193 +/- 72 microg/l, respectively) (PMID: 11380497) 
    %   MW_IGFBP4 : 28 kDa = 28e3 g/mol = ng/nanomol
    %   IGF1BP4 = 190 ng/mL
    %           = 190 / 28e3 * 1e3 nanomole/L = nanomolarity  
    %           = 6.7857 nanomolarity

    % IFGBP5: (PMID: 24379630) 11.3 (8.0-44.6) ng/mL 
    %   MW_IGFBP5 : 31 kDa = 28e3 g/mol = ng/nanomol
    %   IGF1BP5 = 11.3 ng/mL 
    %          = 11.3 / 31e3 * 1e3 nanomole/L = nanomolarity  
    %          = 0.345 nanomolarity

    % IGFBPc = IGFBP3+IGFBP5 = 64.345 nanomolarity 

    % please check parameters2.xlsx and also 


    % Steady state values
    nva.ssIGF1   = 15.6883 ;    % nM - not sure
    nva.ssIGFBP  = 80;          % nM - not sure
    nva.ssIGFBPc = 120;         % nM - not sure
    nva.ssPAPPA  = 40/200;       % PAPPA
    nva.ssPAPPA2 = 0.00132;        % PAPPA2 
    nva.rp = 128.7343;%(0.17/24)*15e-02/1;  %colleti etal 2020 (rat values)0.17/24, divide by 7 as at homestasis tumor factor is about 7
    nva.ra = (0.148/24)*15e-02/1; %0.148/24 % 2x growth over 6month
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

    nva.Fup = 1; % 5e-4; % where did you get this ???? for mAb this is usually closer to 1.0 

    nva.exportFlag = false;
end

nva.kdegIGFBP   = log(2)/nva.trIGFBP ; 
nva.kdegIGFBPc  = log(2)/nva.trIGFBPc ;
nva.kdegPAPPA   = log(2)/nva.trPAPPA;
nva.kdegPAPPA2  = log(2)/nva.trPAPPA2;
nva.kdegIGF1    = log(2)/nva.trIGF1 ; 

% at steady state - no mAb
% d(IGFBP)/dt     = 0 = ksynIGFBP  - kdegIGFBP*IGFBP   - (Kcat*IGFBP*PAPPA/(Km+PAPPA))
% d(IGFBPc)/dt    = 0 = ksynIGFBPc - kdegIGFBPc*IGFBPc - (Kcat2*IGFBPc*PAPPA2/(Km2+PAPPA2))
% d(IGF1)/dt      = 0 = ksynIGF1   - kdegIGF1*IGF1     + (Kcat*IGFBP*PAPPA/(Km+PAPPA)) + (Kcat2*IGFBPc*PAPPA2/(Km2+PAPPA2)) 
% d(PAPPA)/dt     = 0 = ksynPAPPA  - kdegPAPPA*PAPPA
% d(PAPPA2)/dt    = 0 = ksynPAPPA2 - kdegPAPPA2*PAPPA2


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
ed.addspecies('Pain',      'Value', 1,              'Units', 'item');
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
% FOR lesion growth /apoptosis
m1.addparameter('ra',    'Value', nva.ra,              'Units','1/hour');
m1.addparameter('rp',   'Value', nva.rp,              'Units','item/hour');


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
%IGF1 and tumor (EC) and pain
m1.addparameter('EmaxIGF1',    'Value', 1.2990,              'Units','dimensionless');
m1.addparameter('EC50IGF1',   'Value', 0.3535,              'Units','nanomolarity');
% m1.addparameter('kel_pain', 'Value', 1, 'Units', 'dimensionless');
m1.addparameter('hcIGF1',4.3,'Units','dimensionless');
m1.addreaction(' null  -> EC  ', 'ReactionRate', 'rp*((1+ (EmaxIGF1 * (IGF1/EC50IGF1)^hcIGF1))/(1+  (IGF1/EC50IGF1)^hcIGF1))',...
    'Name', 'tumor');
m1.addreaction('EC -> null ', 'ReactionRate', 'ra*EC',...
    'Name', 'basal apoptosisrate Endo');

%reaction rate for modifiers for tumor
m1.addparameter('EmaxIGF1pain',    'Value', 100,              'Units','dimensionless');
m1.addparameter('EC50IGF1pain',   'Value', 56.1907,              'Units','nanomolarity');
%pain and IGF1
m1.addparameter('hcIGF1pain',1.3995,'Units','dimensionless');
m1.addparameter('kpain',   'Value', 0.5,              'Units','item/hour');
m1.addparameter('kelpain', 'Value', 1, 'Units','1/hour');
m1.addreaction('null  -> Pain  ',...
    'ReactionRate', 'kpain *((1+ (EmaxIGF1pain * (IGF1/EC50IGF1pain )^hcIGF1pain)))/((1+  (IGF1/EC50IGF1pain )^hcIGF1pain))',...
     'Name', 'Pain pathway');
m1.addreaction('Pain -> null', 'ReactionRate', 'kelpain*Pain',...
     'Name', 'pain clearance');
% Unbound drug
m1.addparameter('Fup', 'Value', nva.Fup, 'Units', 'dimensionless');
m1.addparameter('dosenM', 'Value', 0, 'Units', 'nanomolarity');
addrule(m1, 'mAb = Fup*dosenM', 'RuleType', 'initialAssignment');
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