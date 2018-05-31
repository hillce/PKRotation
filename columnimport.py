def columnimport(searchfolder, fileending, column_index, delimiter):

#=========================================================================================================================
# Imports values of a specified column from all text files in the searchfolder and subfolders with a specified fileending
# into two lists.
# One output list is called "all_values", which contains all values from all files.
# The other list is called "values_by_file" and contains sublists of values from each separate files 
#=========================================================================================================================
                             
    import os
    all_paths =[]
    for root, dirs, files in os.walk(searchfolder):
        for file in files:
            if file.endswith(fileending):
                all_paths.append(os.path.join(root,file))
    print('Files imported:')
    for path in all_paths:
        print(path)

    values = []
    values_by_file = []
    all_values = []
    for path in all_paths:
        with open(path,"r") as openfile:
            lines = openfile.readlines()
            for x in lines:
                all_values.append(x.split(delimiter)[column_index])
                values.append(x.split(delimiter)[column_index])
        values_by_file.append(values)
        values = []
 
    for i in range(len(all_values)):
        all_values[i] = float(all_values[i])
        
    j = 0
    for value_list in values_by_file:
        for i in range(len(value_list)):
            value_list[i] = float(value_list[i])
        
        values_by_file[j] = value_list
        j +=1

    print('Values imported:')
    for sublist in values_by_file:
        print(len(sublist))
    return values_by_file, all_values