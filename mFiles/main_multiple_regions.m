function [err_code] = main_multiple_regions(PrepConfig,AlgoConfig,SampleConfig,Header_uni,Sequence_uni,taxa_name_calls,ranks_to_extract)

err_code = 0;

% Configs processing
sample_dir = SampleConfig.dname;
primers_seq = SampleConfig.primers_seq;


% Get the sample name
files = dir([sample_dir '/*fastq*']);
if isempty(files) || 2*fix(length(files)/2) ~= length(files)
    err_code = 1;
    return
end
sample_name = extract_sample_name(files(1).name);

% Clean the results directory
resDir = [sample_dir '/resDir'];
if exist(resDir,'dir')
    rmdir(resDir, 's')
end
mkdir(resDir)


% Add primers handling params to AlgoConfig
AlgoConfig.primers_len = PrepConfig.with_primer_flag*zeros(size(primers_seq)) + (1-PrepConfig.with_primer_flag)*cellfun(@length,primers_seq);
AlgoConfig.with_primer_flag = PrepConfig.with_primer_flag;
if ~isfield(AlgoConfig,'use_regions')
    AlgoConfig.use_regions = (1:size(primers_seq,1));
end


% Reads stats struct
readsStatsObj = ReadsStats(size(primers_seq,1));


% Algo PE flag
algo_pe_flag = PrepConfig.algo_pe_flag;
if PrepConfig.pe_flag == 0
    algo_pe_flag = 0;
end


% Quality filter
read_fastq_save_unireads(sample_dir, sample_name, PrepConfig, algo_pe_flag, readsStatsObj)


% Split reads to regions
read_unireads_save_split_to_regions(resDir, sample_name, primers_seq, PrepConfig.max_err_inprimer, PrepConfig.with_primer_flag, algo_pe_flag, readsStatsObj)


% Save the read count statistics
matlab_filename = [resDir '/sample_' sample_name '_sampPrepReadStats.mat'];
save(matlab_filename, 'readsStatsObj')


% Community reconstruction
reconstruction_func(resDir, sample_name, AlgoConfig, readsStatsObj);



% Save the read count statistics
matlab_filename = [resDir '/sample_' sample_name '_readStats.mat'];
save(matlab_filename, 'readsStatsObj')

% Save the reconstruction parameters
matlab_filename = [resDir '/sample_' sample_name '_paramsStructs.mat'];
save(matlab_filename, 'PrepConfig','AlgoConfig','SampleConfig')


% Save the reconstruction
save_reconstruction_new_nogroups(sample_dir, sample_name, taxa_name_calls, ranks_to_extract, Header_uni,Sequence_uni,AlgoConfig.write_groups_fasta)

