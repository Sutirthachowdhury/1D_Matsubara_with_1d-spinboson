Program cxxaverage
      implicit none

      integer l,d,lf,nd
      !parameter(lf=2000)
      parameter(nd=99)
      real*8 time,sum11,ave11
      real*8 state1
      character(len=80) :: fn


      !do l=1,lf
         sum11=0.d0


         do d=0,nd
            write(fn,10)d
            open(d,file=fn)
            read(d,*) state1
            
            sum11=sum11+state1
                       
         enddo
         ave11=sum11/(nd+1)
         
         write(121,'(5g15.7)') ave11
      !enddo

10    format('mapp-',I2.2,'/fort.3001')
    end Program cxxaverage
