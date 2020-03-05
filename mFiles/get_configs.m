function [PrepConfig,AlgoConfig,DbFiles,primers_seq] = get_configs(db_dir, kmer_len)


% Set the GG DB files
DbFiles.uniS16_dir = db_dir; %'/data1/fefuks/Green_Genes_201305/unique_up_to_3_ambiguous_16S';
DbFiles.file_prefix = 'GreenGenes_201305_unique_up_to_3_ambiguous_16S';
DbFiles.rdp_name_calls_file = 'taxonomy_db.mat';


% 5 Regions FFPE multiplex primers
primers_seq{1,1} = 'TGGCGAACGGGTGAGTAA';     primers_name{1,1} = '143-16S_DS-F1b';
primers_seq{1,2} = 'CCGTGTCTCAGTCCCARTG';    primers_name{1,2} = '145-16S_DS-R1a';
primers_seq{2,1} = 'ACTCCTACGGGAGGCAGC';     primers_name{2,1} = '130-16S_DS-F2 ';
primers_seq{2,2} = 'GTATTACCGCGGCTGCTG';     primers_name{2,2} = '131-16S_DS-R2 ';
primers_seq{3,1} = 'GTGTAGCGGTGRAATGCG';     primers_name{3,1} = '132-16S_DS-F3 ';
primers_seq{3,2} = 'CCCGTCAATTCMTTTGAGTT';   primers_name{3,2} = '133-16S_DS-R3 ';
primers_seq{4,1} = 'GGAGCATGTGGWTTAATTCGA';  primers_name{4,1} = '134-16S_DS-F4 ';
primers_seq{4,2} = 'CGTTGCGGGACTTAACCC';     primers_name{4,2} = '135-16S_DS-R4 ';
primers_seq{5,1} = 'GGAGGAAGGTGGGGATGAC';    primers_name{5,1} = '149-16S_DS-F5b';
primers_seq{5,2} = 'AAGGCCCGGGAACGTATT';     primers_name{5,2} = '153-16S_DS-R5c';

DB_kmer_len = 160;
if kmer_len > DB_kmer_len
    kmer_len = DB_kmer_len;
end


% Set the kmers DB parameters
allowed_mm = 2;

% Generate kmers DB path and filename
suffix = ['_ffpe5regions_' num2str(allowed_mm) 'mm_RL' num2str(DB_kmer_len)];
dbPath = DbFiles.uniS16_dir;
dbFileName = [DbFiles.file_prefix suffix];

% suffix = ['_ffpe5regions_' num2str(allowed_mm) 'mm_RL' num2str(DB_kmer_len)];
% dbPath = [DbFiles.uniS16_dir suffix];
% 
% dot_ind = find(DbFiles.headers_filename == '.',1,'last');
% dbFileName = [DbFiles.headers_filename(1:dot_ind-1) suffix];
% 

% ********************** SAMPLE PREP PARAMETERS ********************
PrepConfig.read_len = kmer_len;

% Quality filter
PrepConfig.qual_th = 30;
PrepConfig.prc_high_qual = 0.75; 
PrepConfig.low10_th = 3;
PrepConfig.pe_flag = 1;
PrepConfig.max_num_Ns = 0; 

% Algo related
PrepConfig.algo_pe_flag = 1;
PrepConfig.max_err_inprimer = 2;
PrepConfig.with_primer_flag = 0;


% ********************** ALGORITHM  PARAMETERS ********************
% DB path
AlgoConfig.dbPath = dbPath;
AlgoConfig.dbFileName = dbFileName;

% General parameters
AlgoConfig.verbose = 1;
AlgoConfig.write_groups_fasta = 0;


% Reads filter params
AlgoConfig.filter_reads = 1; % NOTICE: We cancel this filter for LONG reads
AlgoConfig.min_read_freq = 1e-4;
AlgoConfig.min_read_count = 2;

% Reads to kmers alignment params
AlgoConfig.pe = 0.005;
AlgoConfig.nMM_cut = 2; 

% Reconstruction params
AlgoConfig.do_filter = 1;
AlgoConfig.filter_included_bacteria = 1;

AlgoConfig.read_type = char('PE'*PrepConfig.algo_pe_flag + 'SE'*(1-PrepConfig.algo_pe_flag));
AlgoConfig.readLen = kmer_len;

AlgoConfig.tol = 5e-7;
AlgoConfig.numIter = 10000;
