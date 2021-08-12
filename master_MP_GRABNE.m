% analysis the GRAB_NE signal 

clearvars;
close all;
setup_figprop;

root_path = 'F:\GRAB_analysis';
%% matching pennies behavior 

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('-----------------------------------------------------------');

% Look for data files and create a database index
logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex = MP_createBehMatFiles(dataIndex);

%% fluorescent signal preprocessing
MP_GRAB_preprocess(dataIndex);

%% simple plots
MP_GRAB_simpleplots(dataIndex);
