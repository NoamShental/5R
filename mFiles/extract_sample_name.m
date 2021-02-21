function [sample_name] = extract_sample_name(file_name)


% Get the sample name
try
    % L_ind = find(files(ff).name(1:end-1) == '_' & files(ff).name(2:end) == 'L',1);
    L_ind = strfind(file_name,'_L0');
    R_ind = find(file_name(1:end-1) == '_' & file_name(2:end) == 'R',1);
    if isempty(L_ind) || isempty(R_ind)
        sample_name = '';
    else
        str_ind = find(file_name == '_');
        jj = find(str_ind(1:end-1)==L_ind & str_ind(2:end)==R_ind);
        sample_name = file_name(1:str_ind(jj)-1);
    end
catch
    error([file_name ' is not a valid Illumina file name(e.g. SAMPLE_NAME_L001_R1_001.fastq).'])
end

