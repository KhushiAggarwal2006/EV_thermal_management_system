%% =============================================================
%  EV CABIN HVAC - FORWARD / INVERSE PSYCHROMETRIC MODEL
% FORWARD MODEL: INPUT IS THE DESIRED CABIN CONDITION, OUTPUT IS THE EVAPORATOR POWER AND TEMPERATURE
%INVERSE MODEL: INPUT IS THE GIVEN EVAPORATOR POWER AND TEMPERATURE, OUTPUT IS THE RESULTANT CABIN CONDITION (condition means temp and RH)

%%For getting coolprop to work for MATLAB BROWSER VERSION: (steps 1/2/4 need not be run everytime; only step 3 needs to be run if you are refreshing again)
% 1. Download get-pip.py into temp folder
%websave('/tmp/get-pip.py', 'https://bootstrap.pypa.io/get-pip.py');

% 2. Install pip and CoolProp for user account
%system('python3 /tmp/get-pip.py --user');
%system('python3 -m pip install --user -U CoolProp');

% 3. Dynamically locate user site-packages and add to MATLAB's Python
% path-- this needs to be initialised everytime 
%userSiteDir = char(py.site.getusersitepackages());
%if count(py.sys.path, userSiteDir) == 0
%   insert(py.sys.path, int32(0), userSiteDir);
%end

% 4. Refresh Python module cache
%py.importlib.invalidate_caches();

%  Air path:
%    outdoor air + recirculated cabin air -> MIXING
%      -> EVAPORATOR (cool + dehumidify, BF = 0, saturated outlet)
%      -> PTC HEATER (sensible reheat only) -> CABIN
%
%  PROPERTY SOURCES (all via CoolProp's Python interface)
%    Dry air         : CoolProp 'Air'   (real-air Helmholtz EOS, Lemmon et al.),
%                      evaluated at the dry-air partial pressure P - p_v
%    Water / vapour  : CoolProp 'Water' (IAPWS-95 = the steam tables):
%                      p_sat(T), saturated-vapour h_g(T), saturated-liquid
%                      h_f(T), h_fg(T)
%    Humidity ratio  : Dalton mixture, omega = eps*p_v/(P-p_v),
%                      eps = M_water/M_air taken from CoolProp
%
%  ENTHALPY BASIS (kJ per kg dry air):  h = h_a + omega*h_g
%    h_a = 0 for dry air at 0 deg C (same partial pressure)
%    h_g, h_f = 0 for saturated liquid water at the triple point (0.01 deg C)
%
%  ASSUMPTIONS
%    1. Fresh / recirculated fractions are fractions of the total MOIST-air
%       mass flow. Each stream keeps its own T and RH, so the mixing is done
%       on each stream's dry-air flow: mdot_da,i = mdot_moist,i/(1+omega_i).
%    2. BF = 0: coil outlet is saturated at the ADP, so ADP = T_E.
%    3. Condensate leaves at T_E as saturated liquid; its enthalpy is
%       subtracted from the evaporator load.
%    4. Occupant moisture: mdot_w = Q_latent/h_fg(T_skin), h_fg from steam tables.
%    5. Cabin energy balance: mdot_da*(h_R - h_S) = Q_sens + Q_latent.
%    6. Freezing: T_E must stay above T_freeze, otherwise the cycle is
%       infeasible (frost) and the model stops with an error.
%    7. The PTC can only heat: T_S >= T_E must hold, otherwise the supply
%       air would be supersaturated and the model stops with an error.

%  REQUIREMENT: CoolProp installed in the Python environment used by MATLAB.
%  =============================================================
clear; clc;

%% 1. INPUTS ====================================================

% Atmosphere
par.P = 101325;                   % Pressure [Pa]

% Outdoor air
par.T_O   = 35;                   % [deg C]
par.phi_O = 0.50;                 % [-]

% Cabin target (forward model)
par.T_R_target   = 25;            % [deg C]
par.phi_R_target = 0.50;          % [-]

% Airflow through the evaporator
par.mdot_moist      = 420/3600;   % Total moist-air flow [kg/s] (= 420 kg/h)
par.fresh_fraction  = 0.30;       % Fraction of total moist-air mass flow
par.recirc_fraction = 1 - par.fresh_fraction;

% Occupants and envelope
par.n_people        = 4;
par.Q_sens_person   = 70;         % [W/person]
par.Q_latent_person = 45;         % [W/person]
par.Q_vehicle       = 221.35;     % Envelope heat gain [W]
par.T_skin          = 35;         % Perspiration evaporation temperature [deg C]

% Feasibility limit
par.T_freeze = 0.01;              % Water triple point [deg C]

% Optional cross-check against CoolProp HumidAir
run_validation = true;

%% 2. DERIVED QUANTITIES ========================================

par.Q_sens   = par.n_people*par.Q_sens_person + par.Q_vehicle;   % [W]   %MODIFY THIS
par.Q_latent = par.n_people*par.Q_latent_person;                 % [W]   %MODIFY THIS
%par.Q_sens=400;
%par.Q_latent=200;
par.Q_total  = par.Q_sens + par.Q_latent;                        % [W]


par.h_fg   = steamHfg(par.T_skin);                  % [kJ/kg], steam tables
par.mdot_w = par.Q_latent/(par.h_fg*1000);          % Moisture generation [kg/s]

par.out = psychrometricState(par.T_O, par.phi_O, par.P);   % Outdoor state

%% 3. FORWARD MODEL: cabin state -> required evaporator / PTC ===
F = forwardModel(par, par.T_R_target, par.phi_R_target);
printResults('FORWARD MODEL RESULTS', F);

%% 4. INVERSE MODEL: evaporator (Q, T_E) -> cabin state =========
I = inverseModel(par, F.evap.Q, F.evap.T);
printResults('INVERSE MODEL RESULTS', I);

%% 5. FORWARD-INVERSE CONSISTENCY CHECK =========================
% The inverse recovers the cabin state algebraically; the forward model is
% then re-run on it, so Q_evap and T_E must come back as specified.

fprintf('\n========================================\n');
fprintf('  FORWARD-INVERSE COMPARISON\n');
fprintf('========================================\n');
fprintf('Cabin T   : target %.4f  recovered %.4f  error %.2e deg C\n', ...
        par.T_R_target, I.cabin.T, abs(par.T_R_target - I.cabin.T));
fprintf('Cabin RH  : target %.4f  recovered %.4f  error %.2e pct-pts\n', ...
        par.phi_R_target*100, I.cabin.phi*100, ...
        abs(par.phi_R_target - I.cabin.phi)*100);
fprintf('Q_evap    : given  %.4f  recomputed %.4f  error %.2e W\n', ...
        I.spec.Q_evap, I.evap.Q, abs(I.spec.Q_evap - I.evap.Q));
fprintf('T_E       : given  %.4f  recomputed %.4f  error %.2e deg C\n', ...
        I.spec.T_E, I.evap.T, abs(I.spec.T_E - I.evap.T));


%% 6. OPTIONAL VALIDATION AGAINST COOLPROP HUMIDAIR =============
% Differences of ~0.5 % in omega are expected: HumidAir includes the
% water-vapour enhancement factor, this model does not.

if run_validation
    try
        HA = @(out, T_C, key, val) double(py.CoolProp.CoolProp.HAPropsSI( ...
                 out, 'T', T_C + 273.15, key, val, 'P', par.P));

        w_O_ref  = HA('W', par.T_O, 'R', par.phi_O);
        w_R_ref  = HA('W', F.cabin.T, 'R', F.cabin.phi);
        Tdp_ref  = HA('D', F.mix.T, 'W', F.mix.omega) - 273.15;
        dh_ref   = (HA('Hda', par.T_O, 'R', par.phi_O) ...
                  - HA('Hda', F.cabin.T, 'R', F.cabin.phi))/1000;
        dh_mine  = par.out.h - F.cabin.h;

        fprintf('\n========================================\n');
        fprintf('  VALIDATION VS COOLPROP HUMIDAIR\n');
        fprintf('========================================\n');
        fprintf('omega outdoor : model %.6f  HumidAir %.6f  (%+.2f %%)\n', ...
                par.out.omega, w_O_ref, 100*(par.out.omega/w_O_ref - 1));
        fprintf('omega cabin   : model %.6f  HumidAir %.6f  (%+.2f %%)\n', ...
                F.cabin.omega, w_R_ref, 100*(F.cabin.omega/w_R_ref - 1));
        fprintf('T_dp mixed    : model %.3f  HumidAir %.3f deg C\n', ...
                F.mix.T_dp, Tdp_ref);
        fprintf('h_O - h_R     : model %.3f  HumidAir %.3f kJ/kg\n', ...
                dh_mine, dh_ref);
    catch ME
        fprintf('\nValidation skipped: %s\n', ME.message);
    end
end


%% =============================================================
%  LOCAL FUNCTIONS
%  =============================================================


%% ---------------- A. FORWARD MODEL ----------------------------

function F = forwardModel(par, T_R, phi_R)

    P = par.P;

    if phi_R <= 0 || phi_R > 1
        error('hvac:input', 'Cabin RH must be in (0, 1].');
    end

    % 1. Cabin (return) air state
    cabin = psychrometricState(T_R, phi_R, P);

    % 2. Mixing of cabin and outdoor air (moist-air flow basis)
    mix = mixStreams(par, cabin, par.out);
    mix.T    = solveTemperatureFromEnthalpy(mix.h, mix.omega, P, par.T_freeze);
    mix.phi  = relativeHumidity(mix.T, mix.omega, P);

    % 3. Cabin moisture balance:  mdot_da*(omega_R - omega_S) = mdot_w
    omega_S = cabin.omega - par.mdot_w/mix.mdot_air;

    if omega_S <= 0
        error('hvac:moisture', ...
              'Supply humidity ratio is non-positive. Increase airflow or reduce moisture load.');
    end

    if omega_S >= mix.omega
        error('hvac:noDehumidification', ...
              ['No dehumidification needed (omega_S = %.5f >= omega_M = %.5f). ', ...
               'The BF = 0 saturated-outlet evaporator model does not apply.'], ...
              omega_S, mix.omega);
    end

    % 4. Evaporator outlet (PTC does not change humidity: omega_E = omega_S)
    omega_E = omega_S;

    omega_E_min = saturationHumidityRatio(par.T_freeze, P);
    if omega_E < omega_E_min
        error('hvac:frost', ...
              ['INFEASIBLE CYCLE (frosting): required coil-outlet humidity ratio ', ...
               '%.5f is below saturation at %.2f C (%.5f), i.e. the coil would ', ...
               'have to run below the freezing point of water.'], ...
              omega_E, par.T_freeze, omega_E_min);
    end

    T_E = solveSaturationTemperature(omega_E, P, par.T_freeze);   % BF = 0
    h_E = moistAirEnthalpy(T_E, omega_E, P);

    % 5. Condensate and evaporator load (condensate enthalpy included)
    mdot_cond = mix.mdot_air*(mix.omega - omega_E);          % [kg/s]
    h_w       = steamHliquid(T_E);                           % [kJ/kg]

    Q_air  = mix.mdot_air*(mix.h - h_E)*1000;                % Air-side [W]
    Q_cond = mdot_cond*h_w*1000;                             % Carried off by condensate [W]
    Q_evap = Q_air - Q_cond;                                 % [W]

    % 6. Cabin energy balance -> supply state
    h_S = cabin.h - par.Q_total/(mix.mdot_air*1000);
    T_S = solveTemperatureFromEnthalpy(h_S, omega_S, P, par.T_freeze);

    if T_S < T_E - 1e-9
        error('hvac:supersaturated', ...
              ['Required supply temperature (%.2f C) is below the evaporator outlet ', ...
               '(%.2f C): supply air would be supersaturated and a PTC cannot cool. ', ...
               'Cabin target infeasible for this flow and load.'], T_S, T_E);
    end

    % 7. PTC heater
    Q_PTC = mix.mdot_air*(h_S - h_E)*1000;                   % [W]

    % 8. Results
    F.cabin  = struct('T', T_R, 'phi', phi_R, ...
                      'omega', cabin.omega, 'h', cabin.h);

    F.mix    = struct('T', mix.T, 'phi', mix.phi, 'omega', mix.omega, ...
                      'h', mix.h, 'T_dp', solveSaturationTemperature( ...
                                          mix.omega, P, par.T_freeze), ...
                      'mdot_air_R', mix.mdot_air_R, ...
                      'mdot_air_O', mix.mdot_air_O);

    F.evap   = struct('T', T_E, 'omega', omega_E, 'h', h_E, ...
                      'Q_air', Q_air, 'Q_cond', Q_cond, 'Q', Q_evap, ...
                      'mdot_cond', mdot_cond, 'h_w', h_w, ...
                      'frost_margin', T_E - par.T_freeze);
    F.evap.subcooling = F.mix.T_dp - T_E;

    F.supply = struct('T', T_S, 'phi', relativeHumidity(T_S, omega_S, P), ...
                      'omega', omega_S, 'h', h_S);

    F.Q_PTC      = Q_PTC;
    F.mdot_moist = par.mdot_moist;
    F.mdot_air   = mix.mdot_air;

end


%% ---------------- B. INVERSE MODEL ----------------------------

function I = inverseModel(par, Q_evap, T_E)
% Given evaporator load [W] and outlet temperature [deg C], recover the
% cabin state, then re-run the forward model on it to complete the results.

    P = par.P;
    m = par.mdot_moist;
    r = par.recirc_fraction;
    f = par.fresh_fraction;

    if r <= 0
        error('hvac:input', ...
              'Recirculation fraction must be > 0: with 100%% fresh air the cabin state cannot be recovered.');
    end

    if T_E <= par.T_freeze
        error('hvac:frost', ...
              'INFEASIBLE CYCLE (frosting): evaporator outlet %.2f C is at or below %.2f C.', ...
              T_E, par.T_freeze);
    end

    % 1. Evaporator outlet (saturated, BF = 0)
    omega_E = saturationHumidityRatio(T_E, P);
    h_E     = moistAirEnthalpy(T_E, omega_E, P);

    % 2. Cabin humidity ratio from the moisture balance
    %    mdot_da*(omega_R - omega_E) = mdot_w, with
    %    mdot_da = r*m/(1+omega_R) + A,  A = f*m/(1+omega_O)
    %    -> A*x^2 + a1*x + a0 = 0  (x = omega_R); stable positive root below.
    omega_O = par.out.omega;
    A  = f*m/(1 + omega_O);
    a2 = A;
    a1 = r*m + A - A*omega_E - par.mdot_w;
    a0 = -((r*m + A)*omega_E + par.mdot_w);

    omega_R = -2*a0/(a1 + sqrt(a1^2 - 4*a2*a0));

    % 3. Mixed-air state (moist-air flow basis)
    mdot_air_R = r*m/(1 + omega_R);
    mdot_air_O = A;
    mdot_air   = mdot_air_R + mdot_air_O;

    omega_M = (mdot_air_R*omega_R + mdot_air_O*omega_O)/mdot_air;

    if omega_M <= omega_E
        error('hvac:noDehumidification', ...
              'Recovered omega_M (%.5f) <= omega_E (%.5f): inputs do not describe a dehumidifying coil.', ...
              omega_M, omega_E);
    end

    % 4. Evaporator energy balance (condensate enthalpy included)
    %    Q_evap = mdot_da*(h_M - h_E) - mdot_cond*h_w
    mdot_cond = mdot_air*(omega_M - omega_E);
    h_w       = steamHliquid(T_E);
    h_M       = h_E + (Q_evap/1000 + mdot_cond*h_w)/mdot_air;

    % 5. Mixing energy balance -> cabin enthalpy -> cabin temperature
    h_R = (mdot_air*h_M - mdot_air_O*par.out.h)/mdot_air_R;
    T_R = solveTemperatureFromEnthalpy(h_R, omega_R, P, par.T_freeze);

    phi_R = relativeHumidity(T_R, omega_R, P);

    % 6. Complete the picture with the forward model on the recovered state
    I = forwardModel(par, T_R, phi_R);
    I.spec = struct('Q_evap', Q_evap, 'T_E', T_E);

end


%% ---------------- C. MIXING -----------------------------------

function mix = mixStreams(par, cabin, outdoor)
% Mix recirculated cabin air and fresh outdoor air.
% Fractions apply to the total MOIST-air mass flow; each stream's dry-air
% flow follows from its own humidity ratio.

    mdot_R = par.recirc_fraction*par.mdot_moist;     % moist air [kg/s]
    mdot_O = par.fresh_fraction *par.mdot_moist;

    da_R = mdot_R/(1 + cabin.omega);                 % dry air [kg/s]
    da_O = mdot_O/(1 + outdoor.omega);

    mix.mdot_air_R = da_R;
    mix.mdot_air_O = da_O;
    mix.mdot_air   = da_R + da_O;

    mix.omega = (da_R*cabin.omega + da_O*outdoor.omega)/mix.mdot_air;
    mix.h     = (da_R*cabin.h     + da_O*outdoor.h    )/mix.mdot_air;

end


%% ---------------- D. PSYCHROMETRIC STATES ---------------------

function state = psychrometricState(T_C, RH, P)
% State from dry-bulb temperature and relative humidity.

    p_sat = steamPsat(T_C);
    p_v   = RH*p_sat;

    if p_v >= P
        error('hvac:input', 'Vapour pressure (%.0f Pa) must be below total pressure.', p_v);
    end

    omega = molarMassRatio()*p_v/(P - p_v);

    state.T     = T_C;
    state.RH    = RH;
    state.p_sat = p_sat;
    state.p_v   = p_v;
    state.omega = omega;
    state.h     = moistAirEnthalpy(T_C, omega, P);

end


function RH = relativeHumidity(T_C, omega, P)

    p_v = P*omega/(molarMassRatio() + omega);
    RH  = p_v/steamPsat(T_C);

end


function omega_sat = saturationHumidityRatio(T_C, P)

    p_sat     = steamPsat(T_C);
    omega_sat = molarMassRatio()*p_sat/(P - p_sat);

end


function h = moistAirEnthalpy(T_C, omega, P)
% h = h_a + omega*h_g   [kJ/kg dry air]

    p_v = P*omega/(molarMassRatio() + omega);

    if p_v >= P
        error('hvac:input', 'Invalid humidity ratio %.5f.', omega);
    end

    h = airEnthalpy(T_C, P - p_v) + omega*steamHvapour(T_C);

end


%% ---------------- E. PROPERTY LAYER (CoolProp) ----------------

% ---- Air tables (CoolProp 'Air') ----

function h_a = airEnthalpy(T_C, p_a)
% Dry-air enthalpy [kJ/kg], zero at 0 deg C and the same partial pressure.

    h_T   = double(py.CoolProp.CoolProp.PropsSI( ...
                'H', 'T', T_C + 273.15, 'P', p_a, 'Air'));
    h_ref = double(py.CoolProp.CoolProp.PropsSI( ...
                'H', 'T', 273.15,       'P', p_a, 'Air'));

    h_a = (h_T - h_ref)/1000;

end


% ---- Steam tables (CoolProp 'Water', IAPWS-95) ----

function p_sat = steamPsat(T_C)
% Saturation pressure [Pa].

    p_sat = double(py.CoolProp.CoolProp.PropsSI( ...
                'P', 'T', T_C + 273.15, 'Q', 0, 'Water'));

end


function h_g = steamHvapour(T_C)
% Saturated-vapour enthalpy [kJ/kg], zero at saturated liquid, triple point.

    h_g = (double(py.CoolProp.CoolProp.PropsSI( ...
              'H', 'T', T_C + 273.15, 'Q', 1, 'Water')) ...
           - steamRefEnthalpy())/1000;

end


function h_f = steamHliquid(T_C)
% Saturated-liquid enthalpy [kJ/kg], zero at the triple point.

    h_f = (double(py.CoolProp.CoolProp.PropsSI( ...
              'H', 'T', T_C + 273.15, 'Q', 0, 'Water')) ...
           - steamRefEnthalpy())/1000;

end


function h_fg = steamHfg(T_C)
% Latent heat of vaporisation [kJ/kg].

    h_fg = steamHvapour(T_C) - steamHliquid(T_C);

end


function h0 = steamRefEnthalpy()
% Reference: saturated liquid at the triple point (273.16 K). Cached.

    persistent h_ref
    if isempty(h_ref)
        h_ref = double(py.CoolProp.CoolProp.PropsSI( ...
                    'H', 'T', 273.16, 'Q', 0, 'Water'));
    end
    h0 = h_ref;

end


function eps_w = molarMassRatio()
% M_water/M_air from CoolProp (~0.62196). Cached.

    persistent v
    if isempty(v)
        v = double(py.CoolProp.CoolProp.PropsSI('molar_mass', 'Water')) ...
          / double(py.CoolProp.CoolProp.PropsSI('molar_mass', 'Air'));
    end
    eps_w = v;

end


%% ---------------- F. ROOT FINDERS -----------------------------

function T_C = solveTemperatureFromEnthalpy(h_target, omega, P, T_lo)
% Solve h(T, omega) = h_target on [T_lo, 100] deg C.

    fun = @(T) moistAirEnthalpy(T, omega, P) - h_target;

    if fun(T_lo)*fun(100) > 0
        error('hvac:range', ...
              ['State (h = %.3f kJ/kg, omega = %.5f) has T outside [%.2f, 100] C. ', ...
               'Below %.2f C the cycle is infeasible (freezing).'], ...
              h_target, omega, T_lo, T_lo);
    end

    T_C = fzero(fun, [T_lo, 100]);

end


function T_C = solveSaturationTemperature(omega_target, P, T_lo)
% Solve omega_sat(T) = omega_target on [T_lo, 60] deg C.

    fun = @(T) saturationHumidityRatio(T, P) - omega_target;

    if fun(T_lo)*fun(60) > 0
        error('hvac:range', ...
              ['omega = %.5f has no saturation temperature in [%.2f, 60] C ', ...
               '(below %.2f C the cycle is infeasible: freezing).'], ...
              omega_target, T_lo, T_lo);
    end

    T_C = fzero(fun, [T_lo, 60]);

end


%% ---------------- G. REPORTING --------------------------------

function printResults(titleStr, S)

    bar = repmat('=', 1, 40);
    fprintf('\n%s\n  %s\n%s\n', bar, titleStr, bar);

    fprintf('\nCABIN (RETURN) AIR\n');
    fprintf('  Temperature          = %9.3f deg C\n',     S.cabin.T);
    fprintf('  Relative humidity    = %9.3f %%\n',        S.cabin.phi*100);
    fprintf('  Humidity ratio       = %9.6f kg/kg\n',     S.cabin.omega);
    fprintf('  Enthalpy             = %9.3f kJ/kg\n',     S.cabin.h);

    fprintf('\nAIRFLOW\n');
    fprintf('  Moist-air flow       = %9.6f kg/s\n',      S.mdot_moist);
    fprintf('  Dry-air flow         = %9.6f kg/s\n',      S.mdot_air);
    fprintf('    from cabin         = %9.6f kg/s\n',      S.mix.mdot_air_R);
    fprintf('    from outdoor       = %9.6f kg/s\n',      S.mix.mdot_air_O);

    fprintf('\nMIXED AIR (EVAPORATOR INLET)\n');
    fprintf('  Temperature          = %9.3f deg C\n',     S.mix.T);
    fprintf('  Relative humidity    = %9.3f %%\n',        S.mix.phi*100);
    fprintf('  Humidity ratio       = %9.6f kg/kg\n',     S.mix.omega);
    fprintf('  Enthalpy             = %9.3f kJ/kg\n',     S.mix.h);
    fprintf('  Dew point            = %9.3f deg C\n',     S.mix.T_dp);

    fprintf('\nEVAPORATOR OUTLET (SATURATED, ADP = T_E)\n');
    fprintf('  Temperature          = %9.3f deg C\n',     S.evap.T);
    fprintf('  Margin above freezing= %9.3f K\n',         S.evap.frost_margin);
    fprintf('  Humidity ratio       = %9.6f kg/kg\n',     S.evap.omega);
    fprintf('  Enthalpy             = %9.3f kJ/kg\n',     S.evap.h);
    fprintf('  Subcooling below T_dp= %9.3f K\n',         S.evap.subcooling);

    fprintf('\nSUPPLY AIR (AFTER PTC)\n');
    fprintf('  Temperature          = %9.3f deg C\n',     S.supply.T);
    fprintf('  Relative humidity    = %9.3f %%\n',        S.supply.phi*100);
    fprintf('  Humidity ratio       = %9.6f kg/kg\n',     S.supply.omega);
    fprintf('  Enthalpy             = %9.3f kJ/kg\n',     S.supply.h);

    fprintf('\nCONDENSATE\n');
    fprintf('  Rate                 = %9.6f kg/s  (%.4f kg/h)\n', ...
            S.evap.mdot_cond, S.evap.mdot_cond*3600);
    fprintf('  Liquid enthalpy h_w  = %9.3f kJ/kg\n',     S.evap.h_w);

    fprintf('\nCAPACITIES\n');
    fprintf('  Air-side enthalpy removal = %9.3f W\n',    S.evap.Q_air);
    fprintf('  Condensate enthalpy       = %9.3f W\n',    S.evap.Q_cond);
    fprintf('  Evaporator cooling load   = %9.3f W\n',    S.evap.Q);
    fprintf('  PTC heating               = %9.3f W\n',    S.Q_PTC);

end
