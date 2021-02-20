function [bact_freq, bactMetaGroups, keep_col] = solve_iterative_noisy(dat0,Config,readsStatsObj)

if strcmp(Config.read_type,'PE')
    warning('Make sure PE is supported properly')
end

if strcmp(Config.read_type,'PE')
    SE_PE_flag = 1;
elseif strcmp(Config.read_type,'SE')
    SE_PE_flag = 2;
end

nR = length(dat0.A);
nB_all = size(dat0.A{1},2);


% Some preperations 
sumAPerBactPerRegion = zeros(nR,nB_all);
num_perfect_matches =  zeros(nR,nB_all);
for i = 1:nR

    rL = size(dat0.reads{i},2);
    
    % This is the minimal probabilty mass that should be present in reads
    % (it could be larger than this for a cirtain bacteria if some other kmers are mapping to kmer of given bacteria)
    total_pr_thresh = ((1-Config.pe).^(rL-Config.nMM_cut)).*((Config.pe/3).^Config.nMM_cut);
    perfect_match_pr = (1-Config.pe)^rL/SE_PE_flag;
    max_error_pr = total_pr_thresh/SE_PE_flag;
    
    % Remove zero rows of A (even if there is a read there)
    % If we have a row of zeros in A and a non zero read - we assume that this
    % is due to the fact that pair of this read was not observed due to low
    % coverage and we erased the single one from A matrix because if the region
    % was amplified we expect to get 2 reads (for SE) - so now we just erase
    % this row - assuming that the double read (PE) was not observed
    keep_ind = find(sum(dat0.A{i},2)>0);
    
    F_vec{i} = dat0.F{i}(keep_ind);
    A_mat{i} = dat0.A{i}(keep_ind,:);
    
    addStats(readsStatsObj,'Number of reads mathced to DB', length(F_vec{i}), sum(F_vec{i}), i);
    addStats(readsStatsObj,'Count of maximal unmatched read', 1, max([0;dat0.F{i}(sum(dat0.A{i},2)==0)]), i);
    addStats(readsStatsObj,'Count of all unmatched read', sum(sum(dat0.A{i},2)==0), sum(dat0.F{i}(sum(dat0.A{i},2)==0)), i);
    disp(['Region ' num2str(i) ' out of ' num2str(nR)])
    disp(['Keeping reads matched to DB: ' num2str(round(100*length(F_vec{i})/length(dat0.F{i}))) '% of reads'])
    disp(['Keeping reads matched to DB: ' num2str(round(100*sum(F_vec{i})/sum(dat0.F{i}))) '% of counts'])
    disp('--------------------------------------------')
    
    
    
    % Calculate the total observed probabilty per bacteria in region
    sumAPerBactPerRegion(i,:) = sum(A_mat{i},1);
    
    % Calculate the number of perfect matches per region per bacteria
    num_perfect_matches(i,:) = sum(A_mat{i} > perfect_match_pr-0.1*max_error_pr-eps,1);
end


% Set the regions normalization factors
comb_mat = dat0.indInSeqs'>0;
regions_norm_factor = sum(comb_mat,1);


% Filter columns (bacterias)
if Config.verbose
    disp('Filter out columns (bacteria)')
end

missing_0MM_region = zeros(nR,nB_all);
if Config.do_filter == 1
    % We keep bacteria that are present in all the 0MM regions, any of 2MM regions
    missing_0MM_region = (dat0.is_perfect_match' & num_perfect_matches<SE_PE_flag) | ...
        (~dat0.is_perfect_match' & num_perfect_matches>0 & num_perfect_matches<SE_PE_flag);
    
    
    % THIS ROW WAS COMENTED ON 12/02/2016. 
    % I DONT UNDERSTAND WHY IF SE WE DONT FILTER OUT REGIONS WITH ONE READ
    % IT CAUSES FALSE POSITIVES
    %    missing_0MM_region = dat0.is_perfect_match' & num_perfect_matches<SE_PE_flag;
elseif Config.do_filter == 2
    % We keep bacteria that are present in all the 0MM regions, any of 2MM regions
    missing_0MM_region = dat0.is_perfect_match' & num_perfect_matches<Config.num_reads_matched;
end
keep_col = find(sum(missing_0MM_region,1)==0 & sum(sumAPerBactPerRegion,1) > 0 & regions_norm_factor > 0);
nB = length(keep_col);



% Normalize frequency counts
% warning('TAKE PROPER CARE OF NOT AMPLIFIED REGIONS')
if Config.verbose
    disp('Normalize frequency counts')
end
numReadsPerRegion = cellfun(@sum,F_vec);
freq_vec = {};
for i = 1:nR
    
    % Normalize the read counts
    freq_vec{i} = F_vec{i}/sum(numReadsPerRegion);
    
end



% Write one matrix for all the measurements
if Config.verbose
    disp('Build matrix A_L2')
end
A_L2 = zeros(sum(cellfun(@(x) size(x,1),A_mat)),nB);
y_vec = zeros(size(A_L2,1),1);
k = 0;
for i = 1:nR
    A_L2(k+1:k+size(A_mat{i},1),1:nB) = A_mat{i}(:,keep_col);
    y_vec(k+1:k+size(A_mat{i},1)) = freq_vec{i};
    k = k+ size(A_mat{i},1);
end
A_L2 = sparse(A_L2);

% Normalize to regions probability
A_L2 = bsxfun(@rdivide,A_L2,regions_norm_factor(keep_col));

% Find unique columns
if Config.verbose
    disp('Making columns of A unique...')
end
[~,I,J] = unique(A_L2','rows');

% Record the groups
bactGroupsInUnique = [];
for ii = 1:length(I)
    bactGroupsInUnique(ii).group = keep_col(J==ii);
    bactGroupsInUnique(ii).unique = keep_col(I(ii));
end

% It is important to mark all the bacteria as present at the end since due
% to filter could be mismatch
A_L2 = A_L2(:,I);
keep_col = keep_col(I);
nB = length(keep_col);

    
% Filter out bacterias "included" in other bacteria
if Config.filter_included_bacteria == 1
    if Config.verbose
        disp('Removing included bacterias...')
    end
    A_L2_TF = double((A_L2 > 0));
    num_non_zeros = full(sum(A_L2_TF,1));
    n_ii_jj_mat = full(A_L2_TF'*A_L2_TF);
    remove_included = false(1,nB);
    for ii = 1:nB
        jj_ii_diff = num_non_zeros - num_non_zeros(ii);
        this_row_min_n = num_non_zeros - (jj_ii_diff>0).*jj_ii_diff;
        if sum((n_ii_jj_mat(ii,:)-this_row_min_n) == 0 & jj_ii_diff > 0)
            remove_included(ii) = true;
        end
        %         for jj = ii+1:nB
        %             if n_ii_jj_mat(ii,jj) == this_row_min_n(jj)
        %                 if num_non_zeros(ii) < num_non_zeros(jj)
        %                     remove_included(ii) = true;
        %                 elseif num_non_zeros(ii) > num_non_zeros(jj)
        %                     remove_included(jj) = true;
        %                 end
        %             end
        %         end
    end
    A_L2(:,remove_included) = [];
    keep_col(remove_included) = [];
    bactGroupsInUnique(remove_included) = [];
    if Config.verbose
        disp(['Removed ' num2str(sum(remove_included)) ' out of ' num2str(length(remove_included))])
    end
end

num_reads_per_bact = sum(A_L2>0,1);
non_even_count = sum(num_reads_per_bact~=2*fix(num_reads_per_bact/2));
if Config.verbose
    warning(['Found ' num2str(non_even_count) ' bacterias with non even number of reads mapped'])
    disp('Starting iterations...')
end


% *- *- *- *- Reconstruct the bacterial comunity -* -* -* -* -*
[bact_freq, use_col_ind] = ml_em_iterative(A_L2, y_vec, Config);
use_ind = find(bact_freq > 1e-10);
bact_freq = bact_freq(use_ind);
keep_col = keep_col(use_col_ind(use_ind));
bactGroupsInUnique = bactGroupsInUnique(use_col_ind(use_ind));
A_L2 = A_L2(:,use_col_ind(use_ind));


% Record the groups
bactMetaGroups = [];
for ii = 1:length(bactGroupsInUnique)
    %     bactMetaGroups(ii).orig_db_ind = bactGroupsInUnique(ii).group;
    bactMetaGroups(ii).db_ind = bactGroupsInUnique(ii).group;
    bactMetaGroups(ii).comb_vec = comb_mat(:,keep_col(ii))';
end


bact_freq = bact_freq./regions_norm_factor(keep_col)';
bact_freq = bact_freq/sum(bact_freq);
% --------------------- DONE --------------------

% Calculate the log - likelihood removing one col at a time
if 0
    theta_mat = zeros(size(A_L2));
    for jj = 1:size(A_L2,2)
        tmp_Q_mat = [A_L2(:,1:jj-1) A_L2(:,jj+1:end)];
        [tmp_bact_freq, use_col] = ml_em_iterative(tmp_Q_mat , y_vec, Config);
        theta(:,jj) = tmp_Q_mat*tmp_bact_freq;
        LL(jj) = y_vec'*log(theta(:,jj)+eps);
    end
end


