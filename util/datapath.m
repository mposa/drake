function outpath = datapath(file)
wd = pwd; 
outpath = strcat(wd(1:strfind(pwd,'/drake')), 'research_data/', wd(strfind(pwd,'/drake')+1:end),'/',file);
end

