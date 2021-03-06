function save_reconstruction_new_nogroups(readsPath, sample_num, taxa_name_calls, levels_list, Header_uni,Sequence_uni,write_groups_fasta)

% Load the algorithm results
matlab_filename = [readsPath '/resDir/sample_' sample_num '_results.mat'];
load(matlab_filename, 'found_bacteria','bactMetaGroups')


% **************** Build the RECONSTRUCTIONS file ************
if 0 %isempty(taxa_name_calls)
    warning(['TAXA file was NOT FOUND. Not saving reconstruction file to disk for sample: ' sample_num])
else

    L = length(levels_list);
    
    Groups = struct('freq', [], 'reads', [], 'fractions', [], 'hash', [], 'headers', [], 'answer_cell', []);
    Groups(1) = [];
    n_bact = length(found_bacteria.frequency);
    for uu = 1:n_bact
        
        
        bact_ind = sort(bactMetaGroups(uu).db_ind);
        this_group_taxa_table = cell2table(taxa_name_calls(bact_ind,:));
        
        Groups(end+1).freq = found_bacteria.frequency(uu);
        Groups(end).reads = found_bacteria.assigned_reads(uu);
        for ll = 1:L
            [U,~,J] = unique(this_group_taxa_table(:,1:ll),'rows');
            full_names = table2cell(U);
            N = hist(J,1:size(U,1));
            fractions = N'/sum(N);
            Groups(end).fractions{ll} = fractions;
            for rr = 1:length(fractions)
                Groups(end).headers{ll}{rr} = [Header_uni{bact_ind(J==rr)}];
                Groups(end).hash{ll}{rr} = datahash(num2str(bact_ind(J==rr)));
            end
            Groups(end).answer_cell{ll} = [cellstr(num2str(Groups(end).freq*fractions)),cellstr(num2str(round(Groups(end).reads*fractions))), full_names];
        end
        
    end
    
    % Save the MAT file with reconstructed groups
    groups_filename = [readsPath '/resDir/sample_' sample_num '_reconstruction_new_nogroups.mat'];
    save(groups_filename,'Groups','levels_list')
end
% ----------------------------------------------------------------------------------

% **************** Write GROUPS sequences to disk ************
if ~exist('write_groups_fasta','var') || isempty(Sequence_uni)
    write_groups_fasta = 0;
end

if write_groups_fasta
    if ~exist([readsPath '/resDir/groups'],'dir')
        mkdir([readsPath '/resDir/groups'])
    else
        delete([readsPath '/resDir/groups/*.fasta'])
    end
    
    % Write the sequences fasta (group by group)
    if ~isempty(bactMetaGroups)
        for ii = 1:length(bactMetaGroups)
            fasta_filename = [readsPath '/resDir/groups/sample_' sample_num '_group' num2str(ii) '.fasta'];
            fastawrite(fasta_filename, cellfun(@num2str,Header_uni(bactMetaGroups(ii).db_ind),'UniformOutput',false),Sequence_uni(bactMetaGroups(ii).db_ind))
        end
    end
    % ----------------------------------------------------------------------------------
end



