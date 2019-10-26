function path_py=formatPathToPython(path)
    path_py= regexprep(path, '\','\\\');
end