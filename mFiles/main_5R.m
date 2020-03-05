function main_5R(illumina_files_dir,db_dir, results_filename,kmer_len)

if isdeployed && exist('kmer_len')
  kmer_len = str2num(kmer_len);
end

if nargin<4 || isempty(kmer_len)
    kmer_len = 100;
   
    
end

disp('Input params are: ')
disp(['Working on files in directory: ',illumina_files_dir]);
if kmer_len==100
  disp('Using the default k-mer len of 100bp')
else
  disp(['Reconstruction using kmers of length: ',num2str(kmer_len)])
end


% ************** Split to directories (if required) ****************
batch_dir = split_files2directories(illumina_files_dir);


% *************************** Get CONFIGS **************************
[PrepConfig,AlgoConfig,DbFiles,primers_seq] = get_configs(db_dir, kmer_len);


% **************************** LOAD DBs  ***************************
% Load the 16S sequences
uniS16_file = [DbFiles.uniS16_dir '/' DbFiles.file_prefix '_headers.mat'];
load(uniS16_file,'Header_uni')

% Load the taxonomy
taxaS16_file = [DbFiles.uniS16_dir '/' DbFiles.rdp_name_calls_file];
if ~exist('taxa_name_calls','var') && exist(taxaS16_file,'file')
    load(taxaS16_file,'taxa_name_calls','ranks_to_extract')
else
    taxa_name_calls = {};
    ranks_to_extract = {};
end

% ***************************  RECONSTRUCT SAMPLES  ***********************
directories = dir(batch_dir);
batch_samples_list = {};
for dd = 3:length(directories)
    if directories(dd).isdir == 0
        continue
    end
    
    disp(['WORKING ON SAMPLE ' num2str(ss) ': ' directories(dd).name])
    SampleConfig = struct;
    SampleConfig.dname = [batch_dir '/' directories(dd).name];
    SampleConfig.primers_seq = primers_seq;
   
    main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig,Header_uni,[],taxa_name_calls,ranks_to_extract)
    batch_samples_list{end+1} = SampleConfig.dname;
end

% Save the reconstructions file
scott_format_newer_func(batch_samples_list,results_filename)

