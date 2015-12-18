function write_gridfile(grid_struct,filename,precision_x,precision_y)
% write_gridfile(struct,filename,precision)
%
% Output a structure to a comma delimited file with column headers
%
%           struct : any structure composed of one or more matrices and cell arrays
%         filename : file name
%        precision : Number of decimal places to output with numbers
%
%      Given s:
%
%          s.Alpha = { 'First', 'Second';
%                      'Third', 'Fourth'};
%
%          s.Beta  = [[      1,       2;
%                            3,       4]];
%          
%          s.Gamma = {       1,       2;
%                            3,       4};
%
%          s.Epsln = [     abc;
%                          def;
%                          ghi];
% 
%      STRUCT2CSV(s,'any.csv') will produce a file 'any.csv' containing:
%
%         "Alpha",        , "Beta",   ,"Gamma",   , "Epsln",
%         "First","Second",      1,  2,      1,  2,   "abc",
%         "Third","Fourth",      3,  4,      3,  4,   "def",
%                ,        ,       ,   ,       ,   ,   "ghi",
%    
%      v.0.9 - Rewrote most of the code, now accommodates a wider variety
%              of structure children
%
% Written by James Slegers, james.slegers_at_gmail.com
% Covered by the BSD License
%

FID = fopen(filename,'w');
headers = fieldnames(grid_struct);
num_fieldnames = length(headers);
field_sizes = zeros(num_fieldnames,2);

num_struct_elements = length(grid_struct);

for current_element = 1:num_struct_elements
  %  fprintf('%f percent done writing grid...', rr/t*100);
    l = '';
    for rowindex = 1:num_fieldnames
        field_sizes(rowindex,:) = size(grid_struct(current_element).(headers{rowindex}));   
        if ischar(grid_struct(current_element).(headers{rowindex}))
            field_sizes(rowindex,2) = 1;
        end
        l = [l,'"',headers{rowindex},'",'];%,repmat(',',1,sz(ii,2)-1)]; %%HERE
    end

    l = [l,'\n'];

    fprintf(FID,l);

    max_num_field_rows = max(field_sizes(:,1));

    for rowindex = 1:max_num_field_rows
        l = '';
        for current_field_index = 1:num_fieldnames
            c = grid_struct(current_element).(headers{current_field_index});
            str = '';
            
            %If current field has at less than the current rowindex number of rows, don't write anything for this
            if field_sizes(current_field_index,1) < rowindex 
                str = repmat(',',1,field_sizes(current_field_index,2));
            else
                if isnumeric(c)
                    if findstr('cornerlocs_x',headers{current_field_index}) %If this is a grid x node number
                        for kk = 1:field_sizes(current_field_index,2)
                            str = [str,num2str(c(rowindex,kk),precision_x),','];
                        end
                    elseif findstr('cornerlocs_y',headers{current_field_index}) %If this is a grid y node number
                        for kk = 1:field_sizes(current_field_index,2)
                            str = [str,num2str(c(rowindex,kk),precision_y),','];
                        end
                    else
                        for kk = 1:field_sizes(current_field_index,2)
                            str = [str,num2str(c(rowindex,kk)),','];
                        end
                    end
                elseif islogical(c)
                    for kk = 1:field_sizes(current_field_index,2)
                        str = [str,num2str(double(c(rowindex,kk))),','];
                    end
                elseif ischar(c)
                    str = ['"',c(rowindex,:),'",'];
                elseif iscell(c)
                    if isnumeric(c{1,1})
                        for kk = 1:field_sizes(current_field_index,2)
                            str = [str,num2str(c{rowindex,kk},precision),','];
                        end
                    elseif islogical(c{1,1})
                        for kk = 1:field_sizes(current_field_index,2)
                            str = [str,num2str(double(c{rowindex,kk})),','];
                        end
                    elseif ischar(c{1,1})
                        for kk = 1:field_sizes(current_field_index,2)
                            str = [str,'"',c{rowindex,kk},'",'];
                        end
                    end
                end
            end
            l = [l,str];
        end
        l = [l,'\n'];
        fprintf(FID,l);
    end
    fprintf(FID,'\n');
end

fclose(FID);
