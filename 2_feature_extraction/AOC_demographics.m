% AOC Demographics
demos = readtable('/Volumes/methlab_vp/OCC/AOC/AOC_VPs.xlsx');
demos = demos(:, {'ID', 'Gender', 'Alter', 'H_ndigkeit', 'Ausbildung_1_6_Andere_', 'OcularDominance'});
demos = table2struct(demos(1:120, :));

histogram([demos.Alter])