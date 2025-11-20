clear; clc;


T = readtable('beh/socdem01.txt');

keys = {
    '4TH GRADE;', '5TH GRADE;', '7TH GRADE;', '8TH GRADE;', '9TH GRADE;', ...
    '10TH GRADE;', '11TH GRADE;', 'HIGH SCHOOL GRADUATE;', 'SOME COLLEGE, NO DEGREE;', ...
    'ASSOCIATE DEGREE: OCCUPATIONAL, TECHNICAL, OR VOCATIONAL PROGRAM;', ...
    'BACHELOR''S DEGREE (EXAMPLE: BA, AB, BS, BBA);', ...
    'MASTER''S DEGREE (EXAMPLE: MA, MS, MEng, MEd, MBA);', ...
    'DOCTORAL DEGREE (EXAMPLE:PhD, EdD);'
    };
vals = [4, 5, 7, 8, 9, 10, 11, 12, 14, 14, 16, 20, 22];

[tf1, idx1] = ismember(T.cg1_bkgrnd_education, keys);
ed_cg1 = nan(height(T), 1);
ed_cg1(tf1) = vals(idx1(tf1));

[tf2, idx2] = ismember(T.ptner_grade, keys);
ed_ptner = nan(height(T), 1);
ed_ptner(tf2) = vals(idx2(tf2));

caregiver_ed = max([ed_cg1, ed_ptner], [], 2, 'omitnan');

mom_ed = double(T.mother_edu_cat);
dad_ed = double(T.father_edu_cat);

adult_ed = max([mom_ed, dad_ed], [], 2, 'omitnan');

T.parental_education_years = caregiver_ed;
missing_rows = isnan(T.parental_education_years);
T.parental_education_years(missing_rows) = adult_ed(missing_rows);

inc_str = string(T.annual_fam_inc);
inc_str(inc_str == "-999999") = "NaN";
inc_str(inc_str == "999") = "NaN";
inc_str(inc_str == "98") = "NaN";
inc_str(inc_str == "") = "NaN";
T.family_income = str2double(inc_str);

valid_rows = ~isnan(T.parental_education_years) & ~isnan(T.family_income);

T.parental_ses_score = nan(height(T), 1);

z_ed = zscore(T.parental_education_years(valid_rows));
z_inc = zscore(T.family_income(valid_rows));

T.parental_ses_score(valid_rows) = z_ed + z_inc;

writetable(T, 'SES.csv');