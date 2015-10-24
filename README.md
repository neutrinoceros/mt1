# mt1
N bodies scholar collaboration project in fortran90. Work under progress.

**Authors** : Clémnent Robert (clement.robert@protonmail.com) & Riwan Kherouf, supervised by Valery Lainey

![illustration](img/6m7bodies_t.png?raw=true)


##Goals

Learn numerical methods for ephemerids computations.

##Git Memento
**get the code**
  
    $ git clone https://github.com/neutrinoceros/mt1.git

**...or *refresh* your copy**

    $ cd mt1/
    $ git pull

**run the program**

    $ cd mt1/
    $ ./compil.sh
    $ ./program

**submit your changes**
  
    $ git add [all modified files you want to submit]
    $ git commit -m "[short descritpion of the changes]"
    $ git push origin master

##To do

**debbuging**

* ~~none of the planets ever move except for the sun, which trajectory is not even a straight line and depends on integration parameters~~
* ~~check use of RADAU and computations made in subroutine Forces~~
* angular momentum never varies even the slightest, that is to good to be true.
* plottraj.py scripts fails if we use more than 7 bodies (probably an issue in fortran when writing column-heavy files)

**mandatory**

* [ ] verify data units consistency and that of the subroutines (careful with masses in GM and "forces" computations) (*done, but needs a second opinion to ensure all calculations were done right*)
* [ ] make basic tests on the current code :
  - [ ] check deviation for a rapidly evolving planet (Mercury or the Moon) with a 100yr back-and-forth integration to adjust time step (criterion is 100m deviation max)
  - [ ] check conservation of total energy and total angular momentum with 11 bodies (Sun + all planets + moon) (*this one is almost done, we need to redo it after above steps have been checked*)
* [x] get initial-conditions data from pdf (in doc, table 5 and table 8) to .dat
* [x] implement loadIC subroutine to read from the .dat file.
* [x] make Masses and initial condition parameter arrays in module data_planets.f90 generated by data/convtool.py.
* [x] finish implentation of walk subroutine.
* [x] debug writing results.dat

**optional**

* [ ] write a decent Makefile (compilation is currently crudely made with a hard-to-edit shell script compil.sh)
* [ ] check conservation of (a,e,i,omega,Omega,M(t)) for a 2 bodies pb (Sun + Mercury)
* [ ] use double precision

