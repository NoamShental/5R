function [dat0] = build_A_matrices(bactData, reads, AlgoConfig, readsStatsObj)

% if strcmp(AlgoConfig.read_type,'SE') && sum(AlgoConfig.primers_len(:)>0)
%     error('SE without primers is not supported')
% end

if 0 %~AlgoConfig.barcoded_regions
    dat0 = [];
    return
end

pe = AlgoConfig.pe;
nMM_cut = AlgoConfig.nMM_cut;

nR = length(bactData.kmers);
nB = size(bactData.indInSeqs,1);


for rr = 1:nR
    disp(['Region ' num2str(rr) ' out of ' num2str(nR)])
    
    % Filter low abundance reads
    if AlgoConfig.filter_reads == 1
        tot_number = length(reads(rr).uniqueReads_count);
        tot_count = sum(reads(rr).uniqueReads_count);
        filter_thresh = max(AlgoConfig.min_read_freq*tot_count,AlgoConfig.min_read_count);
        remove_ind = reads(rr).uniqueReads_count < filter_thresh;
        reads(rr).uniqueReads(remove_ind,:) = [];
        reads(rr).uniqueReads_count(remove_ind) = [];
        
        addStats(readsStatsObj,'After low abundance filter', length(reads(rr).uniqueReads_count), sum(reads(rr).uniqueReads_count), rr);
        disp(['Keep high freq: ' num2str(round(100*length(reads(rr).uniqueReads_count)/tot_number)) '% of reads'])
        disp(['Keep high freq: ' num2str(round(100*sum(reads(rr).uniqueReads_count)/tot_count)) '% of counts'])
    end
    
    
    % Build M matrix
    if AlgoConfig.verbose
        disp('Building matrix M')
    end
    bamp_in_reg = find(bactData.indInSeqs(:,rr)>0);
    
    if strcmp(AlgoConfig.read_type,'PE')
        Kd = size(bactData.kmers{rr},1);
        rr_M = sparse(bactData.indInSeqs(bamp_in_reg,rr),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd,nB);
        rr_kmers = bactData.kmers{rr};
    elseif strcmp(AlgoConfig.read_type,'SE')
        % Fwd
        fwd_readLen = AlgoConfig.readLen - (1-AlgoConfig.const_len_flag)*AlgoConfig.primers_len(rr,1);
        [values_fwd, ~, indInUni_fwd] = unique(bactData.kmers{rr}(:,1:fwd_readLen),'rows');
        Kd_fwd = length(values_fwd);
        indInValue_fwd = zeros(nB,1);
        indInValue_fwd(bamp_in_reg) = indInUni_fwd(bactData.indInSeqs(bamp_in_reg,rr));
        Ad_fwd = sparse(indInValue_fwd(bamp_in_reg),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd_fwd,nB);
        
        % Rvs
        [values_rvs, ~, indInUni_rvs] = unique(bactData.kmers{rr}(:,fwd_readLen+1:end),'rows');
        Kd = length(values_rvs);
        indInValue_rvs = zeros(nB,1);
        indInValue_rvs(bamp_in_reg) = indInUni_rvs(bactData.indInSeqs(bamp_in_reg,rr));
        Ad_rvs = sparse(indInValue_rvs(bamp_in_reg),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd,nB);
        
        % Write the matrix
        rr_M = [Ad_fwd;Ad_rvs];
        rr_kmers = [values_fwd; values_rvs];
    end
    sumMPerBactPerRegion = sum(rr_M,1);
    rr_M = bsxfun(@rdivide,rr_M,(sumMPerBactPerRegion+eps));
    
    
    % Count errors between reads and kmers
    if AlgoConfig.verbose
        disp('Building matrix A')
    end
    nY = length(reads(rr).uniqueReads_count);
    cell_A_i = cell(1,nY);
    cell_A_j = cell(1,nY);
    cell_A_s = cell(1,nY);
    non_zero_counts = zeros(1,nY);

    nK = size(rr_M,1);
    debug_dist_mat = zeros(nY,size(rr_kmers,1));
    rrL = size(rr_kmers,2);
%     dbfix_thresh = min_dbfix_freq*tot_count;
    for yy = 1:nY
        if AlgoConfig.verbose && yy == 1000*fix(yy/1000)
            disp(['Mapped ' num2str(yy) ' reads out of ' num2str(nY)])
        end
        
        distvec = sum(bsxfun(@ne,rr_kmers,reads(rr).uniqueReads(yy,:)),2)-sum(reads(rr).uniqueReads(yy,:)=='N');
        debug_dist_mat(yy,:) = distvec;
        min_dist_2db = min(distvec);
        
            
        % Augment the database
        if 0 %reads(rr).uniqueReads_count(yy) > dbfix_thresh && min_dist_2db>0 && min_dist_2db<=max_nMM_fix
            disp('OK-fix DB')
            augment_ind = find(distvec == min_dist_2db);
            %             rr_M(end+1,:) = ???
            rr_kmers(end+1,:) = reads(rr).uniqueReads(yy,:);
%             min_dist_2db = 0;
%             distvec(end+1) = 0;
        end
        
        if min_dist_2db <= nMM_cut
            %             Pr_read_given_kmer = ((1-pe).^(rL-distvec)).*((pe/3).^distvec);
            %             Pr_read_given_kmer(distvec>nMM_cut) = 0;
            %             %         Pr_read_given_kmer(Pr_read_given_kmer<p_cut) = 0;
            %             Pr_read_given_j = (Pr_read_given_kmer'*rr_M)./(sumMPerBactPerRegion+eps);
            
            mapped_kmer = distvec<=nMM_cut;
            tmp_Pr_read_given_kmer = ((1-pe).^(rrL-distvec(mapped_kmer))).*((pe/3).^distvec(mapped_kmer));
            Pr_read_given_kmer = sparse(find(mapped_kmer),ones(sum(mapped_kmer),1),tmp_Pr_read_given_kmer,nK,1);
%             Pr_read_given_j = (Pr_read_given_kmer'*rr_M)./(sumMPerBactPerRegion+eps);
            Pr_read_given_j = Pr_read_given_kmer'*rr_M;
            
            mapped_bact = Pr_read_given_j>0;
            cell_A_i{yy} = yy*ones(1,sum(mapped_bact));
            cell_A_j{yy} = find(mapped_bact);
            cell_A_s{yy} = full(Pr_read_given_j(mapped_bact));
            non_zero_counts(yy) = sum(mapped_bact);
        end
        
    end
    
    
    % Build one sparse matrix
    i_all = [cell_A_i{:}];
    j_all = [cell_A_j{:}];
    s_all = [cell_A_s{:}];
    dat0.A{rr} = sparse(i_all,j_all,s_all,nY,nB);
    dat0.reads{rr} = reads(rr).uniqueReads;
    dat0.F{rr} = reads(rr).uniqueReads_count;

    
    disp('--------------------------------------------')
    
    if 0 % AlgoConfig.reduce_memory == 0
        dat0.kmers{rr} = rr_kmers;
        dat0.M{rr} = rr_M;       
    end
end
dat0.is_perfect_match = bactData.is_perfect_match;
% dat0.is_amplified = bactData.indInSeqs>0;
dat0.indInSeqs = bactData.indInSeqs;

