subroutine gravity(a,  b, c)
    integer, intent(in) :: a
    integer, intent(out) ::  b
    integer, intent(out) :: c
    b = a**2
    c = b**2 
    
end subroutine gravity
!hello
program example
real::N
N = 1e6
print*, N
end program example
