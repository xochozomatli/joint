      program test
      common vect(5,5)
      open(unit=1,file="testdat")
      do 2 i=1,4
      do 3 j=1,4
      read(1,*)vect(i,j)
      write(6,*)vect(i,j)
    3 continue
    2 continue
      end
