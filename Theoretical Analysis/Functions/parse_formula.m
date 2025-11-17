function [x, y, z] = parse_formula(formula)
    % Parse CxHyOz formula to extract x, y, z values
    
    x = 0; y = 0; z = 0;
    
    % Extract number after C
    C_match = regexp(formula, 'C(\d*)', 'tokens');
    if ~isempty(C_match)
        if isempty(C_match{1}{1})
            x = 1;
        else
            x = str2double(C_match{1}{1});
        end
    end
    
    % Extract number after H
    H_match = regexp(formula, 'H(\d*)', 'tokens');
    if ~isempty(H_match)
        if isempty(H_match{1}{1})
            y = 1;
        else
            y = str2double(H_match{1}{1});
        end
    end
    
    % Extract number after O
    O_match = regexp(formula, 'O(\d*)', 'tokens');
    if ~isempty(O_match)
        if isempty(O_match{1}{1})
            z = 1;
        else
            z = str2double(O_match{1}{1});
        end
    end
end