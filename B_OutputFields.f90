Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,V,GradP,divV,lapP,divVP,curlV)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P,divV,lapP,divVP,curlV
  Real,Dimension(0:NI,0:NJ,2)::V,GradP
  character(:),allocatable::varss

  varss = 'VARIABLES = "X", "Y", "P", "U", "V", "GradPX", "GradPY", '
  varss = varss // '"DivV", "LapP", "DivVP", "curlV"'
  Write(IO,*) varss
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-40]=CELLCENTERED)'
  Write(IO,'(100F15.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F15.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F15.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F15.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F15.7)') V(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F15.7)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F15.7)') GradP(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F15.7)') divV(1:NI-1,1:NJ-1)
  Write(IO,'(100F15.7)') LapP(1:NI-1,1:NJ-1)
  Write(IO,'(100F15.7)') divVP(1:NI-1,1:NJ-1)
  Write(IO,'(100F15.7)') curlV(1:NI-1,1:NJ-1)
  
End Subroutine 
