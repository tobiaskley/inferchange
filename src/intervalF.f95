subroutine intervalF(n,dec,dep,ilen,nint,nsum,bou) bind(c,name="intervalF_")
         use, intrinsic :: ISO_C_BINDING
         implicit none
         integer(kind= c_int), intent(in)          :: n,dep,nsum
         integer(kind= c_int), intent(in)          :: nint(dep-1)
         integer(kind= c_int)                      :: i,j,k
         real(kind=c_double),  intent(in)          :: dec,ilen(dep-1)
         real(kind=c_double)                       :: ran,temp
         integer(kind= c_int), intent(out)         :: bou(nsum,2)
   10000 j=2
   10002 continue
         bou(1,1)=1
         bou(1,2)=n
   10003 do 10006 i=1,dep-1
         ran = (n-ilen(i))/(nint(i)-1)
         !ran2 = (n-ilen(i))/(nint(i)-1)
         !temp=0
   10004 do 10005 k=1,nint(i)
         bou(j,1)=floor(1+(k-1)*ran)
         bou(j,2)=min(ceiling(ilen(i)+(k-1)*ran-1e-12),n)
         !bou(j,1)=floor(1+temp)
         !bou(j,2)=ceiling(ilen(i)+temp)
         !temp=temp+ran
         j=j+1
   10005 continue
   10006 continue
end subroutine intervalF

