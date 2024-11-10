subroutine gravity(a,  b)
    integer, intent(in) :: a
    integer, intent(out) ::  b
    b = a**2
    
end subroutine gravity
!hello
program example
real::N
N = 1e6
print*, N
end program example
