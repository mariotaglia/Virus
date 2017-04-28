
subroutine readinput
use system
use molecules
use const
use bulk
use MPI
use ellipsoid
use chainsdat
use inputtemp
use mprotein
use saves
use ematrix
use kai
use probe
implicit none

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer ios
integer line, linemax
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
character basura
integer ndi
real*8 ndr

! not defined variables, change if any variable can take the value

ndi = -1e5
ndr = -1.0d10


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
!

vtkflag = ndi
dimx = ndi
dimy = ndi
dimz = ndi
ncha = ndi
long = ndi
cuantas = ndi
infile = ndi
nst = ndi
randominput = ndi
kapd = ndi
Xulimit = ndi
kaptype = ndi

xkapbulk = ndr
dielS = ndr
pHbulk = ndr
dielP = ndr
delta = ndr
csalt = ndr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)print*, 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


 select case (label)

 case ('vtkflag')
   read(buffer, *, iostat=ios) vtkflag
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('probe')
   read(buffer, *, iostat=ios) probenum
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('dimx')
   read(buffer, *, iostat=ios) dimx
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('delta')
   read(buffer, *, iostat=ios) delta
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)


 case ('dimy')
   read(buffer, *, iostat=ios) dimy
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('dimz')
   read(buffer, *, iostat=ios) dimz
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('ncha')
   read(buffer, *, iostat=ios) ncha
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
   cpp = ncha/size
 

  if(mod(ncha, size).ne.0) then
  print*, 'Cannot divide', size, 'processors among ',ncha, 'chains'
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
  endif


 case ('long')
   read(buffer, *, iostat=ios) long
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('cuantas')
   read(buffer, *, iostat=ios) cuantas
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('dielP')
   read(buffer, *, iostat=ios) dielP
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('dielS')
   read(buffer, *, iostat=ios) dielS
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('csalt')
   read(buffer, *, iostat=ios) csalt
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('pHbulk')
   read(buffer, *, iostat=ios) pHbulk
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('infile')
   read(buffer, *, iostat=ios) infile
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('nst')
   read(buffer, *, iostat=ios) nst
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
   read(fh, *)(st0(i),i=1,nst)

 case ('randominput')
   read(buffer, *, iostat=ios) randominput
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('Kapd')
   read(buffer, *, iostat=ios) Kapd
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

   if(mod(Kapd,2).ne.1) then
   print*, 'Kapd should be odd... stoping'
   stop
   endif

 case ('Monomerskap')
 read(fh, *)(ntypes(i),i=1,N_monomer)

 case ('xkapbulk')
   read(buffer, *, iostat=ios) xkapbulk
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('Xulimit')
   read(buffer, *, iostat=ios) Xulimit
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

 case ('kaptype')
   read(buffer, *, iostat=ios) kaptype
   if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

   select case (kaptype)
    case(1) 
     read(fh, *) basura
     read(fh, *)NNN

     if(NNN.ne.0) then

     call allocateell
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), Rell(1,j), Rell(2,j), Rell(3,j)
     if(rank.eq.0)print*,'Set particle',j,'pos to',  Rell(1,j), Rell(2,j), Rell(3,j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *), rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *), rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         print*,'Set particle',j,'rotation to:'
         print*, rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         print*, rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         print*, rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), echarge(j)
     if(rank.eq.0)print*,'Set particle',j,'charge to', echarge(j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), eeps(j)
     if(rank.eq.0)print*,'Set particle',j,'hydrophobicity to', echarge(j)
     enddo

     endif ! NNN

    case(2)
      NNN = 1
      call allocateell
      read(fh, *), basura
      read(fh, *), Rell(1,1), Rell(2,1), Rell(3,1)
      if(rank.eq.0)print*,'Set particle',1,'pos to',  Rell(1,1), Rell(2,1), Rell(3,1)
      read(fh, *), basura
      read(fh, *), rotmatrix(1,1,1), rotmatrix(1,2,1), rotmatrix(1,3,1)
      read(fh, *), rotmatrix(2,1,1), rotmatrix(2,2,1), rotmatrix(2,3,1)
      read(fh, *), rotmatrix(3,1,1), rotmatrix(3,2,1), rotmatrix(3,3,1)
      if(rank.eq.0) then
         print*,'Set particle',1,'rotation to:'
         print*, rotmatrix(1,1,1), rotmatrix(1,2,1), rotmatrix(1,3,1)
         print*, rotmatrix(2,1,1), rotmatrix(2,2,1), rotmatrix(2,3,1)
         print*, rotmatrix(3,1,1), rotmatrix(3,2,1), rotmatrix(3,3,1)
      endif
   endselect
endselect
endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
! 


if(vtkflag.eq.ndi)call stopundef('vtkflag')
if(dimx.eq.ndi)call stopundef('dimx')
if(dimy.eq.ndi)call stopundef('dimy')
if(dimz.eq.ndi)call stopundef('dimz')
if(ncha.eq.ndi)call stopundef('ncha')
if(long.eq.ndi)call stopundef('long')
if(cuantas.eq.ndi)call stopundef('cuantas')
if(infile.eq.ndi)call stopundef('infile')
if(nst.eq.ndi)call stopundef('nst')
if(randominput.eq.ndi)call stopundef('randominput')
if(kapd.eq.ndi)call stopundef('kapd')
if(Xulimit.eq.ndi)call stopundef('Xulimit')
if(kaptype.eq.ndi)call stopundef('kaptype')

if(delta.eq.ndr)call stopundef('delta')
if(xkapbulk.eq.ndr)call stopundef('xkapbulk')
if(dielS.eq.ndr)call stopundef('dielS')
if(dielP.eq.ndr)call stopundef('dielP')
if(csalt.eq.ndr)call stopundef('csalt')
if(pHbulk.eq.ndr)call stopundef('pHbulk')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

subroutine stopundef(namevar)
character(len=*) :: namevar
print*, 'Variable ', namevar, ' is undefined '
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
end

