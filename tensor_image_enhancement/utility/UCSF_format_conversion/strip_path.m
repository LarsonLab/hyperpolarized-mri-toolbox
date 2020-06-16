function filename = stip_path(fullname)

filename = fullname;

filename = fullname;
for (i = length(fullname):-1:1)
  if (fullname(i) == '/')      
    filename = fullname(i+1:end);
    break;
  end
end