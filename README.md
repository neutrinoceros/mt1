# mt1
n bodies scholar collaboration project in fortran90

authors : Clément Robert and Riwan Kherouf, supervised by Valery Lainey

##To do

* finish implentation of walk subroutine
* make Masses a parameter array
* verify data units consistency and that of the subroutines (careful with masses in GM and "forces" computations) 
* use the RADAU integrator (or use MERCURY) 
* make basic tests on the current code

##Optional 

* write a decent Makefile (compilation is currently crudely made with a hard-to-edit shell script compil.sh)

##Done recently 

* got initial-conditions data from pdf (in doc, table 5 and table 8) to .dat
* implemented loadIC subroutine to read from the .dat file
