function scott_format_newer_func(batch_directories_list, results_filename, scott_levels_2save)

sl_ind = find(results_filename=='/',1,'last');

qiime_taxa_format = 0;
write_counts_flag = 0;
print_resolution = 0;
reduce_sample_name = 0;
nR = 5;
nStats = 3;
nRegStats = 5;

min_freq = 0;
cut_freq_th = 1.1; 


nS = length(batch_directories_list);

% STATS ?????
stats_mat = zeros(nStats,nS);
regions_stats_mat = zeros(nRegStats,nR,nS);
samples_names = cell(1,nS);


% Build Scott list
levels_list = {'domain','phylum','class','order','family','genus','species','groups'};
% scott_levels_2save = {'domain','phylum','class','order','family','genus','species','groups'};

for sct = 1:length(scott_levels_2save)
    disp(['Building the Scott files for level: ' scott_levels_2save{sct}])
    [TF,tl] = ismember(scott_levels_2save{sct},levels_list);
    
    % Load the Groups and select the correct data
    all_answer_cell = cell(1,nS);
    all_gr_headers = cell(1,nS);
    all_passed_filt = cell(1,nS);
    samples_names = cell(1,nS);
    for nn = 1:length(batch_directories_list)
        
        fastq_files = dir([batch_directories_list{nn} '/*fastq*']);
        nn_sample_name = extract_sample_name(fastq_files(1).name);
        
        % Load Groups of one sample
        clear Groups
        load([batch_directories_list{nn} '/resDir/sample_' nn_sample_name '_reconstruction_new_nogroups.mat'],'Groups')
        
        
        % Generate the answer cell
        if isempty(Groups)
            all_answer_cell{nn} = {};
            all_gr_headers{nn} = {};
            group_freq = [];
        else
            % Hard Decision
            hd = min(tl,7);
            fractions = arrayfun(@(x) x.fractions{hd},Groups,'UniformOutput',false);
            max_fraction = arrayfun(@(x) max(x.fractions{hd}),Groups,'UniformOutput',false);
            max_ind = cellfun(@(x,y) find(x==y,1),fractions,max_fraction,'UniformOutput',false);
            
            answer_cells = arrayfun(@(x) x.answer_cell{hd},Groups,'UniformOutput',false);
            tmp_answer_cell = cellfun(@(x,y) x(y,:),answer_cells,max_ind,'UniformOutput',false);
            if tl == 8
                hash_cell = arrayfun(@(x) datahash(x.hash{hd}),Groups,'UniformOutput',false);
                all_answer_cell{nn} = [cat(1,tmp_answer_cell{:}) hash_cell'];
            else
                all_answer_cell{nn} = cat(1,tmp_answer_cell{:});
            end
            
            headers = arrayfun(@(x) x.headers{hd},Groups,'UniformOutput',false);
            %                         all_gr_headers{nn} = cellfun(@(x,y) x{y},headers,max_ind,'UniformOutput',false);
            all_gr_headers{nn} = cellfun(@(x) [x{:}],headers,'UniformOutput',false);
            
            group_freq = [Groups.freq];
            all_answer_cell{nn}(:,1) = cellfun(@num2str,num2cell(group_freq),'UniformOutput',false);
        end
        
        
        % Keep only significant frequencies (those that pass the min_freq threshold (instead of what was done after the scott reduction below)
        remove_ind = find(group_freq < min_freq);
        all_answer_cell{nn}(remove_ind,:) = [];
        all_gr_headers{nn}(remove_ind) = [];
        group_freq(remove_ind) = [];
        
        
        % Keep only significant frequencies (those that sum to "cut_freq_th" of the frequency) (instead of what was done at reconstruction)
        [sorted_freq,If] = sort(group_freq,'descend');
        cut_ind = find(cumsum(sorted_freq)>cut_freq_th,1);
        if isempty(cut_ind)
            cut_ind = length(sorted_freq);
        end
        passed_freq_thresh = sort(If(1:cut_ind));
        all_passed_filt{nn} = false(size(all_answer_cell{nn},1),1);
        all_passed_filt{nn}(passed_freq_thresh) = true;
        
        
        % Samples names
        if ~reduce_sample_name
            samples_names{nn} = nn_sample_name;
        else
            s_ind = find(nn_sample_name=='_',1,'last');
            samples_names{nn} = nn_sample_name(1:s_ind-1);
        end
        
        
        % Load read count stats
        load([batch_directories_list{nn} '/resDir/sample_' nn_sample_name '_readStats.mat'])
        if ~isempty(readsStatsObj.stats_counts)
            nr = length(readsStatsObj.stats_counts);
            stats_mat(1:nr,nn) = readsStatsObj.stats_counts;
        end
        if ~isempty(readsStatsObj.regions_stats_counts)
            [nr,nc] = size(readsStatsObj.regions_stats_counts);
            regions_stats_mat(1:nr,1:nc,nn) = readsStatsObj.regions_stats_counts;
        end
        
        disp(['Loaded ' samples_names{nn}])
    end
    row_headers = readsStatsObj.row_headers;
    regions_row_headers = readsStatsObj.regions_row_headers;
    
    % Number of reads per sample
    num_reads_per_sample = zeros(1,nS);
    str_reads_per_sample = cell(1,nS);
    for ss = 1:nS
        if ~isempty(all_answer_cell{ss})
            num_reads_per_sample(ss) = sum(cellfun(@str2num,all_answer_cell{ss}(:,2)));
            str_reads_per_sample{ss} = num2str(num_reads_per_sample(ss));
        else
            str_reads_per_sample{ss} = '0';
        end
    end
    
    
    [new_compact_freq_mat,avg_num_in_scott_cell,gr_in_scott_cell,rows_description] = ...
        build_scott_list_new(tl,all_answer_cell,all_gr_headers,all_passed_filt);
    
    [N,M] = size(new_compact_freq_mat);
    
    % Save the reconstruction file to matlab format
    %     save_filename = [save_dir '/' save_name '_GroupsHeaders_' upper(scott_levels_2save{sct}) '_cutFreq' num2str(min_freq) '.mat'];
    %     save(save_filename,'-v7.3','gr_in_scott_cell','new_compact_freq_mat','samples_names','rows_description','num_reads_per_sample','avg_num_in_scott_cell')
    
    
    % Save to file   
    if write_counts_flag == 0
        % Save to file frequencies
        new_compact_save_cell = [[{'Total # of reads'} cellfun(@(x) '',cell(1,tl-1),'UniformOutput',false) str_reads_per_sample];[levels_list(1:tl) samples_names];[rows_description num2cell(new_compact_freq_mat)]];
    else
        % Save to file reads
        new_compact_count_mat = round(bsxfun(@times,new_compact_freq_mat,num_reads_per_sample));
        if qiime_taxa_format == 1
            q_pfx = 'kpcofgs';
            qiime_names = cellfun(@(x) ['k__' x],rows_description(:,1),'UniformOutput',false);
            for tt = 2:tl
                qiime_names = cellfun(@(x,y) [x ';' q_pfx(tt) '__' y],qiime_names,rows_description(:,tt),'UniformOutput',false);
            end
            qiime_names = cellfun(@(x) strrep(x,' ','_'),qiime_names,'UniformOutput',false);
            new_compact_save_cell = [[{'Total # of reads',''} str_reads_per_sample];[{'Taxonomy','Count'} samples_names];[qiime_names num2cell([sum(new_compact_count_mat,2) new_compact_count_mat])]];
        else
            new_compact_save_cell = [[{'Total # of reads'} cellfun(@(x) '',cell(1,tl),'UniformOutput',false) str_reads_per_sample];[levels_list(1:tl) {'Count'} samples_names];[rows_description num2cell([sum(new_compact_count_mat,2) new_compact_count_mat])]];
        end
    end
    curr_level_results_filename = [results_filename(1:sl_ind) upper(scott_levels_2save{sct}) '_' results_filename(sl_ind+1:end)];
    saveCellFile(new_compact_save_cell, curr_level_results_filename)
    
    % Save to file group sizes
    if print_resolution == 1
        new_compact_save_cell = [[{'Total # of reads'} cellfun(@(x) '',cell(1,tl-1),'UniformOutput',false) str_reads_per_sample];[levels_list(1:tl) samples_names];[rows_description num2cell(round(avg_num_in_scott_cell))]];
        saveCellFile(new_compact_save_cell, [results_filename(1:sl_ind) 'GroupSizes_' upper(scott_levels_2save{sct}) '_' results_filename(sl_ind+1:end)])
    end
    
end
% %------------------------------------------------------------
% 
% Save to disk regions read count stats
stats_header_row = row_headers';
for rr = 1:nR
    stats_header_row = [stats_header_row strcat(regions_row_headers',cellfun(@num2str,num2cell(rr*ones(1,nRegStats)),'UniformOutput',false))];
end


% Stats
compact_count_mat = nan(nS, nStats + nR*nRegStats);
[TF,loc] = ismember(stats_header_row,row_headers);
compact_count_mat(:, TF) = stats_mat(loc(TF),:)';

% Regions stats
for rr = 1:nc
    tmp_mat = squeeze(regions_stats_mat(:,rr,:))';
    compact_count_mat(:,nStats+(rr-1)*nRegStats+(1:nRegStats)) = tmp_mat;
end
    

% compact_count_mat = [stats_mat' reshape(shiftdim(regions_stats_mat,2),[nS nRegStats*nR])];
compact_save_cell = [[{'Sample name'}; samples_names'] [stats_header_row; cellfun(@num2str,num2cell(compact_count_mat),'UniformOutput',false)]];
saveCellFile(compact_save_cell, [results_filename(1:sl_ind) 'ReadCountStats_' results_filename(sl_ind+1:end)])
% %------------------------------------------------------------
