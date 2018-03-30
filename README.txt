There are SIX pieces of code in this project currently:

1. "tov_rad_eulerian" solves the TOV equations parameterised by 'r' for a polytropic equation of state using Euler's method.

2. "tov_rad_eulerian2" is a modification of (1) and has the added utility of reading a tabulated EOS.

3. "tov_rad" solves the TOV equations parameterised by 'r' using RK4 algorithm.

4. "enth_calc" computes the central specific enthalpy for a given value of central density and EOS. Currently, this only works for a polytropic EOS.

5. "tov_enth" solves the TOV equations parameterised by specific enthalpy using RK4 algorithm. Works for polytropic EOS only.

6. "tov_rad_interpol" is an extension of (3) to use tabulated EOS.

Apart from the above codes and the README, the repository contains the following two files.

7. "EOS.txt" contains the table for SLy4 EOS.

8. "Theory.pdf" explains the theoretical background and the equations used in this code.