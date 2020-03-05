function [dest_dir] = split_files2directories(data_dir)

files = dir(data_dir);
is_dir_vec = [files.isdir];
if sum(is_dir_vec) > 2
    dest_dir = data_dir;
    return
end

dest_dir = [data_dir '/Samples'];


% ALL FILES IN ONE DIRECTORY
files = dir([data_dir '/*.gz']);
if isempty(files)
    files = dir([data_dir '/*.fastq']);
end

for ff = 1:length(files)
    
    file_name = strrep(files(ff).name,'_R3','_R2');
    
    dname = extract_sample_name(file_name);
    if isempty(dname)
        continue;
    end
    
    if ~exist([dest_dir '/' dname],'dir')
        mkdir([dest_dir '/' dname])
    end
    copyfile([data_dir '/' files(ff).name], [dest_dir '/' dname '/' file_name])    
end

