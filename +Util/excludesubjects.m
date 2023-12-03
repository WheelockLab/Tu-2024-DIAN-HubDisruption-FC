function data = excludesubjects(data,exclusionid)
fs = fields(data);
for i = 1:length(fs)
    if size(data.(fs{i}),1)>100 % arbitrary number to determine if the data is formatted into Ncols of subjects
       data.(fs{i})(exclusionid,:) = [];
    end   
end

end