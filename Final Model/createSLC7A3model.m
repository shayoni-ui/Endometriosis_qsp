 function m1 = createSLC7A3model(nva)
arguments

    % ARGININE ss : 80microM = 80000nM : physiological from program team slides
    % turnover for Arginine t1/2 =  1.06 -2  Hours
    % 
    % NO ss: 10nM  https://www.sciencedirect.com/science/article/pii/S0005272899000201
    %turnover for NO= 600 mins = 600/60 = 10 Hours 
   
    % Steady state values
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
    nva.exportFlag = false;
end

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


%nva
%r1




sbiounit('kcell','item',1000);

m1 = sbiomodel('SLC7A3model');
ed = m1.addcompartment('ED', 'Value', 1, 'Units', 'liter', 'constant', 1); % endometrium

% add all species
% ed.addspecies('SLC7A1',      'Value', nva.ssARG,              'Units', 'nanomolarity');
ed.addspecies('Arginine_i',      'Value', 0,              'Units', 'nanomolarity');
ed.addspecies('Arginine_e',      'Value', nva.ssARG,              'Units', 'nanomolarity');
ed.addspecies('NO',            'Value', nva.ssNO,             'Units', 'nanomolarity');
ed.addspecies('EC',       'Value', 20000,               'Units', 'item');
ed.addspecies('Pain',      'Value', 0,              'Units', 'item');

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
% tumor proliferation with Arginine uptake
m1.addparameter('EmaxARG',    'Value', 13,              'Units','dimensionless');
m1.addparameter('EC50ARG',   'Value', 0.2,              'Units','nanomolarity');
% FOR lesion growth /apoptosis
m1.addparameter('raARG',    'Value', nva.raARG,              'Units','1/hour');
m1.addparameter('rpARG',   'Value', nva.rpARG,              'Units','item/hour');
%pain and NO
m1.addparameter('kpainNO',   'Value', 1,              'Units','item/hour');
m1.addparameter('kelpainNO', 'Value', 0.01, 'Units','1/hour');
% Drug
m1.addparameter('Imax',    'Value', 1,              'Units','dimensionless');
m1.addparameter('IC50',   'Value', 1,              'Units','nanomolarity');
%Pain based on NO concentration
m1.addparameter('EmaxpainNO',    'Value', 100,              'Units','dimensionless');
m1.addparameter('EC50painNO',   'Value', 74,              'Units','nanomolarity');


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
%Lesion proliferation due to excell IGF1
%m1.addreaction(' null  -> EC  ', 'ReactionRate', 'rp*((1+ (Emax * (Arginine_i/EC50 )))/(1+  (Arginine_i/EC50 )))',...
%    'Name', 'tumor');
m1.addparameter('hcArg',2,'Units','dimensionless');
m1.addreaction(' null  -> EC  ', 'ReactionRate', 'rpARG*((1+ (EmaxARG * (Arginine_i/EC50ARG)^hcArg))/(1+  (Arginine_i/EC50ARG)^hcArg))',...
    'Name', 'tumor');
m1.addreaction('EC -> null ', 'ReactionRate', 'raARG*EC',...
    'Name', 'basal apoptosisrate Endo');

% Pain due to NO in serum
 m1.addparameter('hcpainNO',4.5,'Units','dimensionless');
 m1.addreaction('null -> Pain  ', 'ReactionRate', 'kpainNO*((1+ (EmaxpainNO * (NO/EC50painNO)^hcpainNO))/(1+  (NO/EC50painNO)^hcpainNO))',...
     'Name', 'Pain synthesis');
 m1.addreaction('Pain  -> null ', 'ReactionRate', 'kelpainNO*Pain',...
     'Name', 'Pain Degradation');

% Unbound drug
m1.addparameter('FupARG', 'Value', nva.FupARG, 'Units', 'dimensionless');
m1.addparameter('dosenM', 'Value', 0, 'Units', 'nanomolarity');
addrule(m1, 'Drug = FupARG*dosenM', 'RuleType', 'initialAssignment');

%addrule(m1, 'Pain = Kpain*NO', 'RuleType', 'repeatedAssignment');

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