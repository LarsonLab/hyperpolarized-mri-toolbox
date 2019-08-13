% This is a function for Andy to get filenames in a folder
function [filenames,len] = A_get_filenames2(in_dir,ext)
	dirs=dir([in_dir,'/*',ext]);
	len=length(dirs);
	dircell=struct2cell(dirs)';
	filenames = dircell(:,1);