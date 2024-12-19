function num=randn_rng(range, int_flag)
    % samples random number from normal dist where range gives -3 and +3
    % std dev marks of the distribution
    %
    % inputs:
    %         range - 1x2 vector where mean of 2 numbers is the mean of the
    %                 dist, and the two numbers are +/-3 std devs
    %         int_flag - true/false if true output is int (default=false)
    % outputs:
    %         num  - random number sampled from normal dist described by
    %                range
    %
    arguments
        range (1,2) {mustBeNumeric}
        int_flag = false
    end
    
    mu = mean(range);
    sig = (mu - range(1)) ./ 3;

    num = sig.*randn + mu;

    if int_flag
        num = round(num);
    end
end
