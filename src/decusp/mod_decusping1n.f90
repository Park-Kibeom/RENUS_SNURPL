module decusping1n
    use const
    implicit none
    
    integer                   :: ntfine=10
    
    type xsset
        real                  :: xstr
        real                  :: xsd
        real                  :: xsa
        real                  :: xst
        real                  :: xsf
        real                  :: xsnf
        real                  :: xskp
    end type
end module