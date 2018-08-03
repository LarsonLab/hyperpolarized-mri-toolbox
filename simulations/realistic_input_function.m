function input_function = realistic_input_function(N, TR, Tarrival, Tbolus);

t = [1:N]*TR;
t_input = t-Tarrival;
input_function = gampdf(t_input,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
