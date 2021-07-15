Program average
      implicit none

      integer l,d,lf,nd
      parameter(lf=3000)
      parameter(nd=99)
      real*8 time,sum11,ave11,ave22,sum22
      real*8 state1,state2
      character(len=80) :: fn


      do l=1,lf
         sum11=0.d0
         sum22=0.0d0


         do d=0,nd
            write(fn,10)d
            open(d,file=fn)
            read(d,*) time,state1,state2
            
            sum11=sum11+state1
            sum22=sum22+state2
                       
         enddo
         ave11=sum11/(nd+1)
         ave22=sum22/(nd+1)
         
         write(120,'(5g15.7)') time,ave11,ave22
      enddo

10    format('mapp-',I2.2,'/fort.200')
    end Program average
