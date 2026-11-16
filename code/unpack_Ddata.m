function out = unpack_Ddata(fname)

% Unpacks a text file with step heating experiment data. 
% 
% Argument fname is just a string with the filename.
%
% This is the general-use version that is adaptable with regard to header
% lines, columns, etc. and so on. 
%
% The file to be unpacked should have the following format: Data lines with
% columns. The first three cols are step number, Temp (Deg C), and time
% (min). Subsequent pairs of columns are amount of diffusant and
% uncertainty therein, ideally the unit should be in moles. The last column
% (col #6) is the data use and is optional. If exist, this column should
% have zeroes and ones, indicate which data step to be used in the
% analysis. If absent the script will automatically generate a generic data
% use column which exclude the first and last step only from the analysis.
%
% Returns structure data:

% data.samplename   % Diffusion sample name
% data.step         % heating step number
% data.TC           % Temp (Deg C)
% data.tmin          % Heating time (min)
% data.M            % Diffusant amount (mol) 
% data.delM         % Diffusant uncertainty (+/-)
% data.n            % Total number of heating steps
% data.use          % Data to be used for diffusion exploration

% This code is provided for the sole purpose of evaluation of an
% accompanying scientific paper preprint entitled ’Diffusion kinetics of
% 3He in pyroxene and plagioclase and applications to cosmogenic exposure
% dating and paleothermometry in mafic rocks’ by Marie Bergelin and three
% other authors. It is not licensed for any other reuse or redistribution.
%
% Written by Marie Bergelin 
% Contact: mbergelin@bgc.org
% Last modified: 2025.02.27

impdata = importdata(fname);

% get sample name
namesplit = strsplit(fname,'.'); 
data.samplename = namesplit{1};

% Make space for the data block
data.step  = [];
data.TC = [];
data.tmin = [];
data.M = [];
data.delM = [];

% Assign the data block
if isnumeric(impdata) % if .txt file does not contain headers
    data.step  = impdata(:,1);
    data.TC = impdata(:,2);
    data.tmin = impdata(:,3);
    data.M = impdata(:,4);
    data.delM = impdata(:,5);
    data.n = length(data.step);

    % generate data to use mask if missing from imported data
    if width(impdata) > 5
        data.use = impdata(:,6);
    else data.use = ones(data.n,1); end

elseif isstruct(impdata) % if .txt contains headers
    data.step  = impdata.data(:,1);
    data.TC = impdata.data(:,2);
    data.tmin = impdata.data(:,3);
    data.M = impdata.data(:,4);
    data.delM = impdata.data(:,5);
    data.n = length(data.step);
    
    % generate data to use mask if missing from imported data
    if width(impdata.data) > 5
        data.use = impdata.data(:,6);
    else data.use = ones(data.n,1); end

end

out = data;