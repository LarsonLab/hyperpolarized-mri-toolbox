% This is a function for Andy to get filenames in a folder
function [filenames,len] = A_get_filenames(in_dir)
	dirs=dir([in_dir,'/*.xlsx']);
	len=length(dirs);
	dircell=struct2cell(dirs)';
	filenames = dircell(:,1);