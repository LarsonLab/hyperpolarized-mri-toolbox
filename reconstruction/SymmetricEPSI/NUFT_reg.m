function output = NUFT_reg(input,tr,st,lambda)
% l2 norm regularized 

scale = sqrt(prod(prod(st.Kd))/numel(input(:)));

if (strcmp(tr,'transp')),
    %input=input/length(input(:));
    output = nufft_adj(input,st)/(sqrt(prod(nufft_st.Kd))))*scale;
    output = output(:);
elseif (strcmp(tr,'notransp')),
    %input=input/length(input(:));
    output = reshape(input,st.Nd);
    output = nufft(output,st)/(sqrt(prod(nufft_st.Kd))))*scale;
    output=output(:);   
else
    error('Transpose flag not appropriately defined');
end

return