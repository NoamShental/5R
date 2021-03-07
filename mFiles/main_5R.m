function main_5R(illumina_files_dir, db_dir, db_name, results_filename, kmer_len, generate_group_results)

filepath = fileparts(results_filename);
if ~exist(filepath, 'dir')
    disp(['Output folder: ' filepath ' does not exist!'])
    return
end
    

if isdeployed && exist('kmer_len', 'var')
  kmer_len = str2num(kmer_len);
end

if isdeployed && exist('generate_group_results', 'var')
  generate_group_results = str2num(generate_group_results);
end

if nargin<5 || isempty(kmer_len)
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
[PrepConfig,AlgoConfig,DbFiles,primers_seq] = get_configs(db_dir, db_name, kmer_len);
if exist('generate_group_results', 'var')
    AlgoConfig.write_groups_fasta = generate_group_results;
end

% **************************** LOAD DBs  ***************************
% Load the 16S sequences and/or headers
if exist('generate_group_results', 'var')
    uniS16_file = [DbFiles.uniS16_dir '/' DbFiles.file_prefix '.mat'];
    load(uniS16_file,'Header_uni','Sequence_uni')
else
    try
        uniS16_file = [DbFiles.uniS16_dir '/' DbFiles.file_prefix '_headers.mat'];
        load(uniS16_file,'Header_uni')
    catch
        uniS16_file = [DbFiles.uniS16_dir '/' DbFiles.file_prefix '.mat'];
        load(uniS16_file,'Header_uni')        
    end
    Sequence_uni = [];
end

% Load the taxonomy
taxaS16_file = [DbFiles.uniS16_dir '/' DbFiles.rdp_name_calls_file];
if ~exist('taxa_name_calls','var') && exist(taxaS16_file,'file')
    load(taxaS16_file,'taxa_name_calls','ranks_to_extract')
    if AlgoConfig.write_groups_fasta == 0
        scott_levels_2save = {'species'};
    else
        scott_levels_2save = {'species','groups'};
    end
    
else
    ranks_to_extract = {'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'};
    taxa_name_calls = cell(length(Header_uni), length(ranks_to_extract));
    [taxa_name_calls{:}] = deal('');
    scott_levels_2save = {'groups'};
end

% ***************************  RECONSTRUCT SAMPLES  ***********************
directories = dir(batch_dir);
batch_samples_list = {};
for dd = 3:length(directories)
    if directories(dd).isdir == 0
        continue
    end
    
    disp(['WORKING ON SAMPLE ' num2str(dd) ': ' directories(dd).name])
    SampleConfig = struct;
    SampleConfig.dname = [batch_dir '/' directories(dd).name];
    SampleConfig.primers_seq = primers_seq;
   
    err_code = main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig,Header_uni,Sequence_uni,taxa_name_calls,ranks_to_extract);
    if err_code == 0
        batch_samples_list{end+1} = SampleConfig.dname;
    else
        disp(['Directory ' SampleConfig.dname ' is not a sample directory with paired end fastq files.'])
    end
end

% Save the reconstructions file
if ~isempty(taxa_name_calls)
    scott_format_newer_func(batch_samples_list,results_filename, scott_levels_2save)
end
