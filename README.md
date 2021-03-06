# mt1
N bodies scholar collaboration project in fortran90. Work under progress.

**Authors** : Clément Robert (clement.robert@protonmail.com) & Riwan Kherouf (kherouf@crans.org), supervised by Valery Lainey

![illustration](img/6m7bodies_t.png?raw=true)


##Goals

Learn numerical methods for ephemerids computations.

##Documentation used
available on [Dropbox](https://www.dropbox.com/sh/48ggibduzgidf6v/AAB1_qRgjvUp0z6cd0Wd-8Wna?dl=0)

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

##install SPICE 

 1. download the [toolkit](http://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html)
 2. move it to the directory that contains mt1/
 3. download the [data file](http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp) and move it to toolkit/data

**NAIF integer code ID :**

[naif.jpl.nasa.gov](http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html)

##To do

**debbuging**

* branch 'ajustement' : le calcul des corrections plante après avoir atteint la ligne de alld.dat dont le numéro égale "N\_EVAL * 3 * N\_BOD", ce qui est la dimension du vecteur "OminusC".
 
* ~~none of the planets ever move except for the sun, which trajectory is not even a straight line and depends on integration parameters~~
* ~~check use of RADAU and computations made in subroutine Forces~~
* check angular momentum and energy computations
* ~~plottraj.py scripts fails if we use more than 7 bodies (probably an issue in fortran when writing column-heavy files)~~
* ~~the moon is not initialized correctly~~

**mandatory**

* [x] verify data units consistency and that of the subroutines (careful with masses in GM and "forces" computations) (*done, but needs a second opinion to ensure all calculations were done right*)
* [ ] make basic tests on the current code :
  - [x] check deviation for a rapidly evolving planet (Mercury or the Moon) with a 100yr back-and-forth integration to adjust time step (criterion is 100m deviation max)
  - [x] check conservation of total energy and total angular momentum with 11 bodies (Sun + all planets + moon) (*this one is almost done, we need to redo it after above steps have been checked*)
  - [x] check conservation of (a,e,i,omega,Omega,M(t)) for a 2 bodies pb (Sun + Mercury)
* [x] get initial-conditions data from pdf (in doc, table 5 and table 8) to .dat
* [x] implement loadIC subroutine to read from the .dat file.
* [x] make Masses and initial condition parameter arrays in module data_planets.f90 generated by data/convtool.py.
* [x] finish implentation of walk subroutine.
* [x] debug writing results.dat
* [ ] optimisation :
  - [ ] : multiplication of double abd simple (2*truc with truc double : baaaaaad !) (Division worse !)
  - [ ] : no test in boucle !!
  - [ ] : always boucle in left column (do i =... // do j =...     A(i,j)          =  BAD idea) 

**optional**

* [ ] write a decent Makefile (compilation is currently crudely made with a hard-to-edit shell script compil.sh)
* [x] use double precision
* [ ] use kepler routine for 11-body problem (and see perihelie evolution)
* [ ] optionnal computation (secular, post-fit, corrections)


TIPS :

* You need compile spicelib (chmod +x makeall.csh) in toolkit