module formats

! line formats for out files
!----------------------------
character(len=30) :: OFMT1   = "(3E30.16E3)",&   ! fmt of ipms.dat and imps_back.dat
     OFMT2   = "(34E30.16E3)",&  ! fmt of traj.dat, traj_back.dat, 
     !        vel.dat, vel_back.dat
     OFMT3   = "(7E30.16E3)",&   ! fmt of traj_spice.dat
     OFMT4   = "(7E30.16E3)",&   ! fmt of 2bodipms_back.dat and _back
     OFTM44  = "(66E30.16E3)"  ! fmt of alld.dat

end module formats
