function ss_xlatebin(root_fname, N)

currentdir = pwd;

for k = 1:N;
  system(sprintf('ssh ese2 ''cd %s; mv %s%03d.rho %s%03d-nohd.rho; xlatebin -o %s%03d.rho %s%03d-nohd.rho'' ', currentdir, root_fname,k, root_fname,k, root_fname,k, root_fname,k));
  system(sprintf('ssh ese2 ''cd %s; mv %s%03d.pha %s%03d-nohd.pha; xlatebin -o %s%03d.pha %s%03d-nohd.pha'' ', currentdir, root_fname,k, root_fname,k, root_fname,k, root_fname,k));
  if k==N
    system(sprintf('ssh ese2 ''cd %s; mv %s%03d.grd %s%03d-nohd.grd; xlatebin -o %s%03d.grd %s%03d-nohd.grd'' ', currentdir, root_fname,k, root_fname,k, root_fname,k, root_fname,k));
  else
      system(sprintf('rm %s%03d.grd', root_fname,k));
  end

  disp(int2str(k))
end

