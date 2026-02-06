function m1 = createEP3model(nva)
arguments

    % PGE2 MW : 352.465 Da = g/mol = ng/nanomol
    % (PNID: 9768660) PGE2 Normal ranges for age-matched subjects are
    % 400 pg/mL = 0.4 ng/ml  = 0.4/352.5 *1e3 = 1.1348nM
    % EP3: 18–274 ng/mL
    %   MW_EP3 :60000 Da = g/mol = ng/nanomol 
    %   EP3    = 100 ng/mL 
    %           = 100 / 60000 * 1e3 nanomole/L = nanomolarity  
    %           = 1.667 nanomolarity


    % please check parameters2.xlsx and also 


    % Steady state values
    nva.ssPGE2   = 1.1348 ;    % nM - not sure
    nva.ssEP3  = 1.667;      % nM - not sure

    nva.kpainEP3 = 0.2821;  %
    nva.kelpainEP3 = 1; %
    % Turnover rates
    nva.trPGE2   = 0.15;       % hours
    nva.trEP3  = 12;          % hours
    nva.kon = 1;
    nva.koff = 5.8; %kd = koff/kon = 5.8nM


    nva.Fup = 1e-4; % 5e-4; % where did you get this ???? for mAb this is usually closer to 1.0 

    nva.exportFlag = false;
end



nva.kdegPGE2   = log(2)/nva.trPGE2 ; 
nva.kdegEP3  = log(2)/nva.trEP3 ;


nva.ksynPGE2  = nva.kdegPGE2  * nva.ssPGE2;
nva.ksynEP3 = nva.kdegEP3 * nva.ssEP3;


nva




sbiounit('kcell','item',1000);

m1 = sbiomodel('EP3model');% change to paw size concentration for invitro
ed2 = m1.addcompartment('ED2', 'Value', 1, 'Units', 'liter', 'constant', 1); % endometrium

% add all species
ed2.addspecies('PGE2',      'Value', nva.ssPGE2,              'Units', 'nanomolarity');
ed2.addspecies('EP3',     'Value', nva.ssEP3,             'Units', 'nanomolarity');
ed2.addspecies('PGE2EP3',     'Value', 0,             'Units', 'nanomolarity');
ed2.addspecies('Pain',      'Value', 0,              'Units', 'item');
% Drug : Mab
ed2.addspecies('Drug',        'Value', 0,                        'Units', 'nanomolarity');
ed2.addspecies('DrugEP3',  'Value', 0,                        'Units', 'nanomolarity');
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
% FOR lesion growth /apoptosis
m1.addparameter('kpainEP3',    'Value', nva.kpainEP3,              'Units','item/hour');
m1.addparameter('kelpainEP3',   'Value', nva.kelpainEP3,              'Units','1/hour');


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
%EP3 and pain
m1.addparameter('EmaxEP3',    'Value', 118,              'Units','dimensionless');
m1.addparameter('ED50EP3',   'Value', 0.6117,              'Units','nanomolarity');%kassuya paper
% m1.addparameter('kel_pain', 'Value', 1, 'Units', 'dimensionless');
m1.addparameter('hcEP3',4.4,'Units','dimensionless');
m1.addreaction(' null  -> Pain  ', 'ReactionRate', 'kpainEP3*((1+ (EmaxEP3 * (PGE2EP3/ED50EP3)^hcEP3))/(1+  (PGE2EP3/ED50EP3)^hcEP3))',...
    'Name', 'Allodynia');
m1.addreaction('Pain -> null ', 'ReactionRate', 'kelpainEP3*Pain',...
    'Name', 'removal of pain');

%reaction rate for modifiers for tumor
% m1.addreaction('EC  -> Pain  ',...
%     'ReactionRate', 'kbasal *((1+ (EmaxA * (IGF1/EC50A1 ))))/((1+  (IGF1/EC50A1 )))',...
%     'Name', 'Pain pathway');
% m1.addreaction('Pain -> null', 'ReactionRate', 'kel_pain*Pain',...
%     'Name', 'pain clearance');
% Unbound drug
m1.addparameter('Fup', 'Value', nva.Fup, 'Units', 'dimensionless');
m1.addparameter('dosenM', 'Value', 0, 'Units', 'nanomolarity');
%addrule(m1, 'Drug = Fup*dosenM', 'RuleType', 'initialAssignment');
addrule(m1, 'ED2.PGE2 =  dosenM', 'RuleType', 'initialAssignment');
% proliferation


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