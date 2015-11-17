module data_planets 

use parameters
	character(len=10),parameter,dimension(N_BOD) :: NAMES=(/ &
		"Sun       ",&
		"Mercury   "/)

	real(8),parameter,dimension(N_BOD) :: MASSES=(/ &
		1.0,&
		1.66011415305e-07/)

	real(8),parameter,dimension(3*N_BOD) :: IPOSITIONS=(/ &
		0.00450250878464055477,0.00076707642709100705,0.00026605791776697764,&
		0.36176271656028195477,-0.09078197215676599295,-0.08571497256275117236/)

	real(8),dimension(3*N_BOD),parameter :: IVELOCITIES=(/ &
		-0.00000035174953607552,0.00000517762640983341,0.00000222910217891203,&
		0.00336749397200575848,0.02489452055768343341,0.01294630040970409203/)

end module data_planets
