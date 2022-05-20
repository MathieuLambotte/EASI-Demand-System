# EASI-Demand-System
R Codes for an EASI Demand System in R, using the linear iterated 3SLS method with IV from Pendakur 2008 (http://www.sfu.ca/~pendakur/EASI%20made%20Easier.pdf).
This code reused some of the code of the discontinued EASI package on the CRAN http://www2.uaem.mx/r-mirror/web/packages/easi/index.html.
I updated the matrix related computations to make it faster and included a way to deal with censored budget share à la Shonkwiler & Yen (https://doi.org/10.2307/1244339).
I adapted the computation of price and income elasticities. 
