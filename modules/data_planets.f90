module data_planets 

use parameters
	character(len=10),parameter,dimension(N_BOD) :: NAMES=(/ &
		"Sun       ",&
		"Mercury   ",&
		"Venus     ",&
		"Earth     ",&
		"Mars      ",&
		"Jupiter   ",&
		"Saturn    ",&
		"Uranus    ",&
		"Neptune   ",&
		"Pluto     ",&
		"Moon      "/)

	real(8),parameter,dimension(N_BOD) :: MASSES=(/ &
		1.0,&
		1.66011415305e-07,&
		2.44783828779e-06,&
		3.00348961492e-06,&
		3.22715603755e-07,&
		0.000954791915634,&
		0.000285885672705,&
		4.36624373585e-05,&
		5.15138377262e-05,&
		7.36178160609e-09,&
		3.69430331069e-08/)


end module data_planets
