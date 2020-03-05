function  [new_compact_freq_mat,avg_num_in_scott_cell,gr_in_scott_cell,TAXA_list] = ...
    build_scott_list_new(tl, all_answer_cell, all_gr_headers, all_passed_filt)

nS = length(all_answer_cell);

% Concatenate all groups into one list
one_compact_cell = cat(1,all_answer_cell{:});
nN = cellfun(@length,all_gr_headers);
if size(one_compact_cell,1) ~= sum(nN)
    error('Number of answers and headers aree not the same')
end
sI = cumsum([0 nN(1:end-1)])+1;
sE = cumsum(nN(1:end));

% Find the unique list to report
full_names_table = cell2table(one_compact_cell(:,3:tl+2));
[Us,~,Js] = unique(full_names_table,'rows');
TAXA_list = cellfun(@(x) strrep(x,'Assigned','Unknown'),table2cell(Us),'UniformOutput',false);
nT = size(TAXA_list,1);

% Map frequencies and groups of each sample on the combined list of bacteria
passed_filt_mat = zeros(nT, nS);
new_compact_freq_mat = zeros(nT, nS);
gr_in_scott_cell = cell(nT, nS);
weighted_num_in_scott_cell = zeros(nT, nS);
for ss = 1:nS
    tic
    if ~isempty(all_answer_cell{ss})
        this_sample_ind = Js(sI(ss):sE(ss));
        this_sample_freq = cellfun(@str2num,all_answer_cell{ss}(:,1));
        for ii = 1:length(this_sample_ind)
            
            new_compact_freq_mat(this_sample_ind(ii),ss) = ...
                new_compact_freq_mat(this_sample_ind(ii),ss) + this_sample_freq(ii);
            
            gr_in_scott_cell{this_sample_ind(ii),ss} = ...
                [gr_in_scott_cell{this_sample_ind(ii),ss} all_gr_headers{ss}{ii}];
            
            weighted_num_in_scott_cell(this_sample_ind(ii),ss) = ...
                weighted_num_in_scott_cell(this_sample_ind(ii),ss) + this_sample_freq(ii)*length(unique(fix(all_gr_headers{ss}{ii})));
            
            passed_filt_mat(this_sample_ind(ii),ss) = ...
                passed_filt_mat(this_sample_ind(ii),ss) + all_passed_filt{ss}(ii);
            
        end
    end
    %     disp(['Mapped sample ' num2str(ss) ' in ' num2str(toc) ' seconds'])
end
% gr_in_scott_cell = cellfun(@unique,gr_in_scott_cell,'UniformOutput',false);
avg_num_in_scott_cell = weighted_num_in_scott_cell./(new_compact_freq_mat+eps);

% Filter empty row (no sample passed the filter)
nS_passed_filt = sum(passed_filt_mat>0,2);
remove_rows = find(nS_passed_filt==0);
new_compact_freq_mat(remove_rows,:) = [];
avg_num_in_scott_cell(remove_rows,:) = [];
gr_in_scott_cell(remove_rows,:) = [];
TAXA_list(remove_rows,:) = [];





