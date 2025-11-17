% Example 1: Methane
Cr = combustion_coefficients('CH4')
% Output: Cr = [-1, -2, 0, 1, 2]

% Example 2: Ethanol
Cr = combustion_coefficients('C2H6O')
% Output: Cr = [-1, -3, 0, 2, 3]

% Example 3: Propane
Cr = combustion_coefficients('C5H9O2')
% Output: Cr = [-1, -5, 0, 3, 4]

% Example 4: Octane (gasoline)
Cr = combustion_coefficients('C8H18')
% Output: Cr = [-1, -12.5, 0, 8, 9]