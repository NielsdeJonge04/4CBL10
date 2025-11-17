function Cr = combustion_coefficients(fuel_formula)
    % Parse combustion reaction and return combustion coefficients
    % fuel_formula: string of fuel formula, e.g., 'CH4'
    % Cr: combustion coefficients for {'Fuel','O2','N2','CO2','H2O'}
    %     (negative for reactants, positive for products)
    
    % Parse fuel formula to get x, y, z (CxHyOz)
    [x, y, z] = parse_formula(fuel_formula);
    
    % Calculate combustion coefficients for complete combustion
    % Reaction: CxHyOz + a*O2 -> x*CO2 + (y/2)*H2O
    a_O2 = x + y/4 - z/2;  % oxygen coefficient
    
    % Initial coefficients: [Fuel, O2, N2, CO2, H2O]
    Cr = [-1, -a_O2, 0, x, y/2];
end