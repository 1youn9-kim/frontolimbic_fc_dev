clear; clc;

%% 1. Load Data
raw = readtable('beh/mab01.txt');

vars = {
    'src_subject_id', 'interview_age', 'sex', ...
    'presmed_alcohol', 'devhis_8_alchohol_avg', 'devhis_9_alchohol_avg', ...
    'presmed_smoke', 'devhis_8_cigs_per_day', 'devhis_9_cigs_per_day', ...
    'cocaine_gestation', 'marijuan', 'heroin_gestation', ...
    'preg4_2', 'devhis_10b_p', 'devhis_10i_p', 'quesmoth_probanemia', 'devhis_10k_p', ...
    'ph_11', 'birth_weight_lbs', 'birth_weight_oz', 'neo_cyanosis', ...
    'devhis_14b_p', 'devhis_14c_p', 'devhis_14d_p', 'devhis_14f_p'
};

data = standardizeMissing(raw(:, vars), {999, -888, '999', "Don't Know"});

w_grams = (fillmissing(data.birth_weight_lbs, 'constant', 0) * 453.592) + ...
          (fillmissing(data.birth_weight_oz, 'constant', 0) * 28.35);
w_grams(w_grams == 0) = NaN;

R = table();

% Alcohol: 'Yes' OR avg drinks > 0
R.Alcohol = ismember(data.presmed_alcohol, 'Yes') | ...
            data.devhis_8_alchohol_avg > 0 | data.devhis_9_alchohol_avg > 0;

% Tobacco: 'Yes' OR cigs/day > 0
R.Tobacco = ismember(data.presmed_smoke, 'Yes') | ...
            data.devhis_8_cigs_per_day > 0 | data.devhis_9_cigs_per_day > 0;

% Drugs: 'Yes' check
R.Cocaine   = ismember(data.cocaine_gestation, 'Yes');
R.Marijuana = ismember(data.marijuan, 'Yes');
R.Heroin    = ismember(data.heroin_gestation, 'Yes');

% Medical
R.LowBirthWeight = w_grams < 2500;
R.Premature      = data.ph_11 == 1;
R.Preeclampsia   = data.preg4_2 == 1;
R.Bleeding       = data.devhis_10b_p == 1;
R.GestDiabetes   = data.devhis_10i_p == 1;
R.Anemia         = ismember(data.quesmoth_probanemia, 'Yes');
R.Placenta       = data.devhis_10k_p == 1;
R.Cyanosis       = ismember(data.neo_cyanosis, 'Yes');
R.SlowHeart      = data.devhis_14b_p == 1;
R.NoBreath       = data.devhis_14c_p == 1;
R.Convulsions    = data.devhis_14d_p == 1;
R.OxygenReq      = data.devhis_14f_p == 1;

risk_mat = table2array(R);

data.prenatal_risk_score = mean(risk_mat, 2, 'omitnan');

data.src_subject_id = cellstr(string(data.src_subject_id));

writetable(data, fullfile(Top, 'derivatives', 'perinatal_risk_score.csv'));