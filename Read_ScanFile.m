function [scans] = Read_ScanFile(data_path, data_prefix)
% Read_ScanFile: Import every Scanfile in the directory
    
    % list files in folder
    file_list = dir(strcat(data_path, filesep, data_prefix, '*.dat'));
    
    % extract file names
    file_names = extractfield(file_list, 'name');
    
    scan_count = length(file_names);
    
    % allocate scans
    scans = cell(1,scan_count);
    
    % read each file
    for i=1:scan_count
        cur_file = strcat(data_path, filesep, file_names{i}); 
        % import the data contained in the file
        fid=fopen(cur_file);
        cur_data = textscan(fid,'%d%f%f%f'); 
        scans{i}=[cur_data{2}, cur_data{3}, cur_data{4}];
        fclose(fid); 
    end

end