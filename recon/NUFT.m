function output = NUFT(input,tr,st)

scale = sqrt(prod(prod(st.Kd))/numel(input(:)));
if (strcmp(tr,'transp')),
    %input=input/length(input(:));
    output = nufft_adj(input,st)*scale;
    output = output(:);
elseif (strcmp(tr,'notransp')),
    %input=input/length(input(:));
    output = reshape(input,st.Nd);
    output = nufft(output,st)*scale;
    output=output(:);   
else
    error('Transpose flag not appropriately defined');
end

return


