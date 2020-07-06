# Greek-Default
"Inflation-Default Trade-off Without a Nominal Anchor: The Case of Greece", joint work with Enes Sunel, https://www.sciencedirect.com/science/article/pii/S1094202520300338
1) Description of files for table and figure replications:

    a. Solving model for CCD economy (Column 4 of Table 2 in the manuscript):
    
        i. Compile and run file CCD.f90
        
    b. Solving model for MCD economy (Column 2 of Table 2 in the manuscript)
    
        i. Compile and run file MCD.f90. After the deviations between two iterations are small and stationary, change lines 1988 and 1540 by setting indicator_external =1, and indicator_global_search = 0. 
        
        ii. Computation time takes around one week in Windows environment.
        
    c. hpfilter.m: Matlab file that filters data using the Hodrick-Prescott filter. 
    
    d. simulate.m: Matlab file that produces the moments that are reported in the tables of the paper. To obtain the CCD economy moments, change indicator_CCD = 1 on line 8. For MCD economy moments, change indicator_CCD = 0.
    
    e. Solving model for alternative parameterizations and to replicate the moments in Table 3.
    
        i. Replace the corresponding parameters in the FORTRAN file. Note that Χ, θ1 and θ2 are denoted by alpha0, psi and psi2, respectively, in the code.

2) How to run the codes

FORTRAN files require having an access to the IMSL library from which the files invoke some subroutines from. Save your files in a subdirectory named "graphs" for each FORTRAN run. Save matlab files in the same subdirectory. Necessary seeds for random number generation are all specified in the corresponding FORTRAN files. We use a random number generator to draw sequences of income shocks and probability of exclusion. We keep the draws so that the same values can be used for each sample. Routine RNSET of IMSL is used to initialize the seed, which is set to be 139719. Data.xlsx includes data on inflation, level and currency composition of sovereign debt and bond yields.
