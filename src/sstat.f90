!************************************
!* subroutine: sstat
!* version: 2.0
!* date: 06/22/2014
!* last modified: 
!* author: Angelo MeleF
!* description: computes the sufficient stats
!*              given a vector of change statistics
!*              
!************************************
subroutine sstat(ps,ns,qs,g,x,dts,ts)

! SSTAT: functions that computes ALL sufficient statistis
!        given a vector of change statistics
!
! INPUT
! ps   = vector with number of suff stats (ppu,ppm,ppv) (INTEGER(3))
! ns   = number of players (INTEGER)
! qs   = number of variables in xxx (columns of xxx), INTEGER
! g    = matrix with the network (integer, dimension(ns,ns))
! x    = matrix with exogenous variables (real, dimesion(ns,qs))
! dts  = vector containing the specification (INTEGER, DIMENSION(ppp,2)
!
! OUTPUT
! ts   = vector of sufficient statistics


implicit none

! inpur vbls
integer, dimension(3), intent(in):: ps
integer, intent(in):: ns
integer, intent(in):: qs
integer, intent(in):: g(ns,ns)
real(8), dimension(ns,qs), intent(in):: x
integer, dimension(ps(3),2), intent(in):: dts

! output vbls
real(8), dimension(ps(3)), intent(out):: ts
! working vbls
integer i,j,j1,k
real(8) dtfu,dtfm,dtfvq


ts(1:ps(3))=0.0


! statistics for u
do i=1,ns
	do j = 1,ns
        if (g(i,j)==1) then
        do k=1,ps(1)
    		ts(k) = ts(k) + dtfu(ns,qs,x,dts(k,1),i,j,dts(k,2))
        enddo !k
        endif
	enddo !j
enddo !i


! statistics for m
do i=1,ns
	do j=i,ns
		if (g(i,j)==1 .and. g(j,i)==1) then	
        do k=ps(1)+1,ps(2)
			ts(k) = ts(k) + dtfm(ns,qs,x,dts(k,1),i,j,dts(k,2)) 
		enddo ! k
		endif
	enddo !j
enddo !i


! statistics for v
do i=1,ns
	do j = 1,ns
        if (g(i,j)==1) then
        do k=ps(2)+1,ps(3)
			ts(k) = ts(k) + dtfvq(ns,qs,g,x,dts(k,1),i,j,dts(k,2))
        enddo !k
        endif !(g(i,j)==1)
	enddo !j
enddo !i

!ts(:) = 10.0
end subroutine sstat



!************************************
!* function: dtfu
!* version: 2.0
!* date: 05/22/2009
!* last modified: 05/22/2009
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics of u_ij
!*              implemented with sparse network matrix
!************************************

function dtfu(nu,qu,x,xcolu,iu,ju,tu)

! dtfu: contains the formula for all the change statistics of u_ij
!
! INPUT
! nu  = number of agents in the network (INTEGER) 
! iu  = player i (INTEGER)
! ju  = player j (INTEGER)
! xcolu = column of x matrix to use
! tu = type of statistics (INTEGER)
!       1 = single link
!       2 = difference: xu(i)-xu(j)
!       3 = absolute difference: |xu(i)-xu(j)|
!       4 = sum: xu(i)+xu(j)
!       5 = indicator: I(x(i)=0,x(j)=0)
!       6 = indicator: I(x(i)=0,x(j)=1)
!       7 = indicator: I(x(i)=1,x(j)=0)
!       8 = indicator: I(x(i)=1,x(j)=1)
!       9 = indicator: I(x(i)=x(j))
!      10 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions

implicit none

!**** input vbls ****
integer, intent(in)::nu,qu
real(8), intent(in), dimension(nu,qu)::x
integer, intent(in):: iu,ju,tu, xcolu
!**** output vbls ****
real(8) dtfu
!**** working vbls ****
!real(8) indic


dtfu=0.0
if (tu==1) then
	dtfu=1.0         !**** edges ****
else if (tu==2) then
	dtfu = x(iu, xcolu) - x(ju, xcolu)       !**** difference ****
else if (tu==3) then
	dtfu = abs(x(iu, xcolu) - x(ju, xcolu))         !**** absolute difference ****
else if (tu==4) then
	dtfu = x(iu, xcolu) + x(ju, xcolu)       !**** sum ****
else if (tu==5) then
	if (x(iu, xcolu)==0 .and. x(ju, xcolu)==0) dtfu = 1.0    !**** indicator both zeros ****
else if (tu==6) then
	if (x(iu, xcolu)==0 .and. x(ju, xcolu)==1) dtfu = 1.0      !**** indicator zero/one ****
else if (tu==7) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==0) dtfu = 1.0      !**** indicator one/zero ****
else if (tu==8) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1) dtfu = 1.0      !**** indicator both ones ****
else if (tu==9) then
	if (x(iu, xcolu)==x(ju, xcolu))  dtfu = 1.0    !**** same ****
else if (tu==10) then 
	if (x(iu, xcolu) .ne. x(ju, xcolu))  dtfu = 1.0     !**** different ****
else if (tu==11) then
	dtfu = x(iu, xcolu)*x(ju, xcolu)        !**** product ****
else if (tu==12) then
	dtfu = x(iu, xcolu)             !**** x_i ****
else if (tu==13) then
	dtfu = x(ju, xcolu)            !**** x_j ****

! **** DIFFERENTIAL HOMOPHILY ***
else if (tu==21) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1)  dtfu = 1.0
else if (tu==22) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==2)  dtfu = 1.0
else if (tu==23) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==24) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==25) then
	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==31) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==32) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==33) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==34) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==35) then
	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==41) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==42) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==43) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==44) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==45) then
	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==51) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==52) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==53) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==54) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==55) then
	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==5) dtfu = 1.0
else if (tu==61) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==1) dtfu = 1.0
else if (tu==62) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==2) dtfu = 1.0
else if (tu==63) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==3) dtfu = 1.0
else if (tu==64) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==4) dtfu = 1.0
else if (tu==65) then
	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==5) dtfu = 1.0

!! **** DIFF HOMOPHILY interacted with SHARES ****
!else if (tu==101) then
!	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1)  dtfu = 1.0*share(1)
!else if (tu==102) then
!	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==2) dtfu = 1.0*share(2)
!else if (tu==103) then
!	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==3) dtfu = 1.0*share(3)
!else if (tu==104) then
!	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==4) dtfu = 1.0*share(4)
!else if (tu==105) then
!	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==5) dtfu = 1.0*share(5)
!else if (tu==111) then
!	if (x(iu, xcolu)==1 .and. x(ju, xcolu)==1) dtfu = 1.0*(share(1)**2)
!else if (tu==112) then
!	if (x(iu, xcolu)==2 .and. x(ju, xcolu)==2) dtfu = 1.0*(share(2)**2)
!else if (tu==113) then
!	if (x(iu, xcolu)==3 .and. x(ju, xcolu)==3) dtfu = 1.0*(share(3)**2)
!else if (tu==114) then
!	if (x(iu, xcolu)==4 .and. x(ju, xcolu)==4) dtfu = 1.0*(share(4)**2)
!else if (tu==115) then
!	if (x(iu, xcolu)==5 .and. x(ju, xcolu)==5) dtfu = 1.0*(share(5)**2)
!
!! **** SHARES 
!else if (tu==201) then
!	dtfu = share(1)
!else if (tu==202) then
!	dtfu = share(2)
!else if (tu==203) then
!	dtfu = share(3)
!else if (tu==204) then
!	dtfu = share(4)
!else if (tu==205) then
!	dtfu = share(5)
!
!! **** SHARES SQUARED
!else if (tu==211) then
!	dtfu = share(1)**2
!else if (tu==212) then
!	dtfu = share(2)**2
!else if (tu==213) then
!	dtfu = share(3)**2
!else if (tu==214) then
!	dtfu = share(4)**2
!else if (tu==215) then
!	dtfu = share(5)**2

endif

!dtfu = dtfu/(nu*(nu-1))
!dtfu = dtfu/nu
!dtfu = dtfu/10000.0
return
end  function dtfu





!************************************
!* function: dtfm
!* version: .0
!* date: 05/18/2009
!* last modified: 05/18/2009
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics of m_ij
!************************************

function dtfm(nm,qm,x,xcolm,im,jm,tm)

! dtfm: contains the formula for all the change statistics of m_ij
!
! INPUT
! nm  = number of agents in the network (INTEGER) 
! qm  = number of columns in x (INTEGER)
! x   = matrix of exogenous variables (REAL, dimension(nm,qm))
! xcolm =  column of matrix x to use
! im  = player i (INTEGER)
! jm  = player j (INTEGER)
! tm  = type of statistics (INTEGER)
!       1 = single link
!       2 = absolute difference: |x(i)-x(j)|
!       3 = sum: xu(i)+xu(j)
!       4 = indicator: I(x(i)=0,x(j)=0)
!       5 = indicator: I(x(i)=1,x(j)=1)
!       6 = indicator: I(x(i)=x(j))
!       7 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions

implicit none

! input vbls
integer, intent(in)::nm,qm
real(8), intent(in), dimension(nm,qm)::x
integer, intent(in):: im,jm,tm, xcolm
! output vbls
real(8) dtfm
! working vbls
!real(8) indic

dtfm=0.0
if (tm==1) then
	dtfm=1.0         !**** edges ****
else if (tm==2) then
	dtfm = abs(x(im, xcolm) - x(jm, xcolm))         !**** absolute difference ****
else if (tm==3) then
	dtfm = x(im, xcolm) + x(jm, xcolm)       !**** sum ****
else if (tm==4) then
	if (x(im, xcolm)==0 .and. x(jm, xcolm)==0) dtfm = 1.0     !**** indicator both zeros ****
else if (tm==5) then
	if (x(im, xcolm)==1 .and. x(jm, xcolm)==1) dtfm = 1.0     !**** indicator both ones ****
else if (tm==6) then
	if (x(im, xcolm)==x(jm, xcolm)) dtfm = 1.0
	!dtfm = indic(x(im, xcolm),x(jm, xcolm))     !**** same ****
else if (tm==7) then 
	if (x(im, xcolm) .ne. x(jm, xcolm)) dtfm = 1.0
	!dtfm = 1.0 - indic(x(im, xcolm),x(jm, xcolm))      !**** different ****
else if (tm==8) then
	dtfm = x(im, xcolm)*x(jm, xcolm)      !**** different ****
else if (tm==9) then
	dtfm = x(im, xcolm)      !**** xi ****
else if (tm==10) then
	dtfm = x(jm, xcolm)      !**** xj ****

! **** DIFFERENTIAL HOMOPHILY ***
else if (tm==21) then
	if (x(im, xcolm)==1 .and. x(jm, xcolm)==1) dtfm = 1.0
else if (tm==22) then
	if (x(im, xcolm)==2 .and. x(jm, xcolm)==2) dtfm = 1.0
else if (tm==23) then
	if (x(im, xcolm)==3 .and. x(jm, xcolm)==3) dtfm = 1.0
else if (tm==24) then
	if (x(im, xcolm)==4 .and. x(jm, xcolm)==4) dtfm = 1.0
else if (tm==25) then
	if (x(im, xcolm)==5 .and. x(jm, xcolm)==5) dtfm = 1.0

endif
!dtfm = dtfm/((nm*(nm-1))/2)
!dtfm = dtfm/(nm*(nm-1))
!dtfm = dtfm/(nm/2)
!dtfm = dtfm/10000.0

return
end  function dtfm


!************************************
!* function: dtfv
!* version: 2.0
!* date: 05/22/2009
!* last modified: 6/25/2014
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics for v_ij  
!*              
!************************************

function dtfv(nv,qv,g,x,xcolv,iv,jv,tvq)

! dtfvq: contains the formula for all the change statistics
!
! INPUT
! nv   = number of agents in the network (INTEGER) 
! qv   = number of columns in x (integer)
! g    = network matrix (integer, dimension(nv,nv))
! x    = matrix with exogenous characteristics (real, dimension(nv,qv))
! xcolv= column of x matrix to use
! iv   = player i (INTEGER)
! jv   = player j (INTEGER)
! tv   = type of statistics (INTEGER)
!       1 = single link
!       2 = difference: x(i)-x(j)
!       3 = absolute difference: |x(i)-x(j)|
!       4 = sum: x(i)+x(j)
!       5 = indicator: I(x(i)=0,x(j)=0)
!       6 = indicator: I(x(i)=0,x(j)=1)
!       7 = indicator: I(x(i)=1,x(j)=0)
!       8 = indicator: I(x(i)=1,x(j)=1)
!       9 = indicator: I(x(i)=x(j))
!      10 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions


implicit none

! input vbls
integer, intent(in)::nv,qv
integer, intent(in), dimension(nv,nv)::g
real(8), intent(in), dimension(nv,qv)::x
integer, intent(in)::xcolv
integer, intent(in):: iv,jv,tvq
! output vbls
real(8) dtfv
! working vbls
integer i0,k0,ch
!real(8) indic

dtfv=0.0


!*********************
!*    single link    *
!*********************
if (tvq==1) then
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfv = dtfv + 1.0
            endif !(g(jv,k0)==1)           
        endif !(k0 .ne. iv)
        
        if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                dtfv = dtfv + 1.0
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)
 
    enddo !k0
	
!********************
!*    difference    *
!********************
else if (tvq==2) then  
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfv = dtfv + ( x(iv, xcolv) - x(k0, xcolv) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)

        if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                dtfv = dtfv + ( x(k0, xcolv) - x(jv, xcolv) )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)

    enddo !k0


!***********************
!* absolute difference *
!***********************
else if (tvq==3) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfv = dtfv + abs( x(iv, xcolv) - x(k0, xcolv) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
        if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                dtfv = dtfv + abs( x(k0, xcolv) - x(jv, xcolv) )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)

 
    enddo !k0

!*************
!*    sum    *
!*************
else if (tvq==4) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfv = dtfv + ( x(iv, xcolv) + x(k0, xcolv) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
       if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                dtfv = dtfv + ( x(k0, xcolv) + x(jv, xcolv) )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)

         
    enddo !k0
	
!*****************************
!*    Indicator zero/zero    *
!*****************************
else if (tvq==5) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolv)==0 .and. x(k0, xcolv) ==0 ) then
                    dtfv = dtfv + 1.0
                endif !( x(iv, xcolv)==0 .and. x(k0, xcolv) ==0 )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
       if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                if ( x(k0, xcolv)==0 .and. x(jv, xcolv) ==0 ) then
                    dtfv = dtfv + 1.0
                endif !( x(k0, xcolv)==0 .and. x(jv, xcolv) ==0 )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)

         
    enddo !k0
	
!!*****************************
!!*    Indicator zero/one     *
!!*****************************
!else if (tvq==6) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==0 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==0)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!	
!!*****************************
!!*    Indicator one/zero     *
!!*****************************
!else if (tvq==7) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==0 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==0 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!	
!*****************************
!*    Indicator one/one      *
!*****************************
else if (tvq==8) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolv)==1 .and. x(k0, xcolv) ==1 ) then
                    dtfv = dtfv + 1.0
                endif !( x(iv, xcolv)==1 .and. x(k0, xcolv) ==1 )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
       if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                if ( x(k0, xcolv)==1 .and. x(jv, xcolv) ==1 ) then
                    dtfv = dtfv + 1.0
                endif !( x(k0, xcolv)==1)*(x(jv, xcolv) ==1 )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)
    enddo !k0

!*****************************
!*    Indicator x(i)=x(j)    *
!*****************************
else if (tvq==9) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolv) == x(k0, xcolv) ) then
                    dtfv = dtfv + 1.0
                endif !( x(iv, xcolv) == x(k0, xcolv) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
       if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                if ( x(k0, xcolv)==x(jv, xcolv) ) then
                    dtfv = dtfv + 1.0
                endif !( x(k0, xcolv)==x(jv, xcolv) )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)
    enddo !k0

!*****************************
!*    Indicator x(i)/=x(j)   *
!*****************************
else if (tvq==10) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolv) .ne. x(k0, xcolv) ) then
                    dtfv = dtfv + 1.0
                endif ! ( x(iv, xcolv) .ne. x(k0, xcolv) )               
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
       if (k0 .ne. jv) then
            if (g(k0,iv)==1) then
                if ( x(k0, xcolv) .ne. x(jv, xcolv) ) then
                    dtfv = dtfv + 1.0
                endif !( x(k0, xcolv) .ne. x(jv, xcolv) )
            endif !(g(k0,iv)==1)
        endif !(k0 .ne. jv)
    enddo !k0

!*****************************
!!*    Product x(i)*x(k)      *
!!*****************************
!else if (tvq==11) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)*x(k0, xcolv) )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)*x(jv, xcolv) )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!*****************************
!!*     x(j)      *
!!*****************************
!else if (tvq==12) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + x(jv, xcolv)
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!        if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + x(k0, xcolv)
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!   enddo !k0
!
!!******************************
!!* Differential Homophily 1/1 *
!!******************************
!else if (tvq==21) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==1)
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/2 *
!!******************************
!else if (tvq==22) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/3 *
!!******************************
!else if (tvq==23) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/4 *
!!******************************
!else if (tvq==24) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 1/5 *
!!******************************
!else if (tvq==25) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==1 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==1)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/1 *
!!******************************
!else if (tvq==31) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/2 *
!!******************************
!else if (tvq==32) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/3 *
!!******************************
!else if (tvq==33) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/4 *
!!******************************
!else if (tvq==34) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/5 *
!!******************************
!else if (tvq==35) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==2 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==2)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/1 *
!!******************************
!else if (tvq==41) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv= dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/2 *
!!******************************
!else if (tvq==42) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/3 *
!!******************************
!else if (tvq==43) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/4 *
!!******************************
!else if (tvq==44) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/5 *
!!******************************
!else if (tvq==45) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==3 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==3)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/1 *
!!******************************
!else if (tvq==51) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/2 *
!!******************************
!else if (tvq==52) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/3 *
!!******************************
!else if (tvq==53) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/4 *
!!******************************
!else if (tvq==54) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/5 *
!!******************************
!else if (tvq==55) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==4 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==4)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/1 *
!!******************************
!else if (tvq==61) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==1 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/2 *
!!******************************
!else if (tvq==62) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==2 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/3 *
!!******************************
!else if (tvq==63) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==3 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/4 *
!!******************************
!else if (tvq==64) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==4 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/5 *
!!******************************
!else if (tvq==65) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfv = dtfv + ( x(iv, xcolv)==5 )*( x(k0, xcolv) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!       if (k0 .ne. jv) then
!            if (g(k0,iv)==1) then
!                dtfv = dtfv + ( x(k0, xcolv)==5)*(x(jv, xcolv) ==5 )
!            endif !(g(k0,iv)==1)
!        endif !(k0 .ne. jv)
!    enddo !k0
!
!else
!	dtfv=0.0
!endif
!

!*********************
!*    cyclic triangles   *
!*********************
elseif (tvq==100) then
    do k0=1,nv
        if (k0 .ne. iv) then
            dtfv = dtfv + 1.0*g(jv,k0)*g(k0,iv)
        endif !(k0 .ne. iv)
 
    enddo !k0
	
else
	dtfv=0.0
endif

!dtfv = dtfv/(nv*(nv-1)*(nv-2))
!dtfv = dtfv/(nv*(nv-1))
!dtfv = dtfv/10000.0

return
end function dtfv




!************************************
!* function: dtfvq
!* version: 2.0
!* date: 05/22/2009
!* last modified: 05/22/2009
!* author: Angelo Mele
!* description: contains the formula for the 
!*              change statistics for v_ij 
!*              ONLY TO COMPUTE POTENTIAL FUNCTION
!*              implemented with sparse network matrix
!************************************

function dtfvq(nv,qv,g,x,xcolvq,iv,jv,tvq)

! dtfvq: contains the formula for all the change statistics
!
! INPUT
! nv   = number of agents in the network (INTEGER) 
! qv   = number of columns in x (integer)
! g    = network matrix (integer, dimension(nv,nv))
! x    = matrix with exogenous characteristics (real, dimension(nv,qv))
! xcolvq= column of matrix x to use
! iv   = player i (INTEGER)
! jv   = player j (INTEGER)
! tv   = type of statistics (INTEGER)
!       1 = single link
!       2 = difference: x(i)-x(j)
!       3 = absolute difference: |x(i)-x(j)|
!       4 = sum: x(i)+x(j)
!       5 = indicator: I(x(i)=0,x(j)=0)
!       6 = indicator: I(x(i)=0,x(j)=1)
!       7 = indicator: I(x(i)=1,x(j)=0)
!       8 = indicator: I(x(i)=1,x(j)=1)
!       9 = indicator: I(x(i)=x(j))
!      10 = indicator: I(x(i) .ne. x(j))
!
! NOTE: if you want to add some statistics you can do it here
!       by adding additional terms to the IF conditions

implicit none

! input vbls
integer, intent(in)::nv,qv
integer, intent(in), dimension(nv,nv)::g
real(8), intent(in), dimension(nv,qv)::x
integer, intent(in):: iv,jv,tvq, xcolvq
! output vbls
real(8) dtfvq
! working vbls
integer i0,k0,ch
!real(8) indic

dtfvq=0.0


!*********************
!*    single link    *
!*********************
if (tvq==1) then
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfvq = dtfvq + 1.0
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0
	
!********************
!*    difference    *
!********************
else if (tvq==2) then  
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfvq = dtfvq + ( x(iv, xcolvq) - x(k0, xcolvq) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0


!***********************
!* absolute difference *
!***********************
else if (tvq==3) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfvq = dtfvq + abs( x(iv, xcolvq) - x(k0, xcolvq) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!*************
!*    sum    *
!*************
else if (tvq==4) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfvq = dtfvq + ( x(iv, xcolvq) + x(k0, xcolvq) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0
	
!*****************************
!*    Indicator zero/zero    *
!*****************************
else if (tvq==5) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if (x(iv,xcolvq)==0 .and. x(k0,xcolvq)==0) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0
	
!*****************************
!*    Indicator zero/one     *
!*****************************
else if (tvq==6) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==0 .and. x(k0, xcolvq) ==1 ) then
					dtfvq = dtfvq + 1.0 
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0
	
!*****************************
!*    Indicator one/zero     *
!*****************************
else if (tvq==7) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==1 .and. x(k0, xcolvq) ==0 ) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0
	
!*****************************
!*    Indicator one/one      *
!*****************************
else if (tvq==8) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==1 .and. x(k0, xcolvq) ==1 ) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!*****************************
!*    Indicator x(i)=x(j)    *
!*****************************
else if (tvq==9) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)== x(k0, xcolvq) ) then
					dtfvq = dtfvq + 1.0
				endif	
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!*****************************
!*    Indicator x(i)/=x(j)   *
!*****************************
else if (tvq==10) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq) .ne. x(k0, xcolvq) ) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!*****************************
!*    Product x(i)*x(k)      *
!*****************************
else if (tvq==11) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfvq = dtfvq + ( x(iv, xcolvq)*x(k0, xcolvq) )
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!*****************************
!*     x(j)      *
!*****************************
else if (tvq==12) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                dtfvq = dtfvq + x(jv, xcolvq)
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!******************************
!* Differential Homophily 1/1 *
!******************************
else if (tvq==21) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==1 .and.  x(k0, xcolvq) ==1) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!******************************
!* Differential Homophily 1/2 *
!******************************
else if (tvq==22) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==1 .and. x(k0, xcolvq) ==2 ) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!******************************
!* Differential Homophily 1/3 *
!******************************
else if (tvq==23) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==1 .and. x(k0, xcolvq) ==3 ) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!******************************
!* Differential Homophily 1/4 *
!******************************
else if (tvq==24) then          
    do k0=1,nv
        if (k0 .ne. iv) then
            if (g(jv,k0)==1) then
                if ( x(iv, xcolvq)==1 .and. x(k0, xcolvq) ==4 ) then
					dtfvq = dtfvq + 1.0
				endif
            endif !(g(jv,k0)==1)
        endif !(k0 .ne. iv)
    enddo !k0

!!******************************
!!* Differential Homophily 1/5 *
!!******************************
!else if (tvq==25) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==1 )*( x(k0, xcolvq) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/1 *
!!******************************
!else if (tvq==31) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==2 )*( x(k0, xcolvq) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/2 *
!!******************************
!else if (tvq==32) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==2 )*( x(k0, xcolvq) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/3 *
!!******************************
!else if (tvq==33) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==2 )*( x(k0, xcolvq) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/4 *
!!******************************
!else if (tvq==34) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==2 )*( x(k0, xcolvq) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 2/5 *
!!******************************
!else if (tvq==35) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==2 )*( x(k0, xcolvq) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/1 *
!!******************************
!else if (tvq==41) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==3 )*( x(k0, xcolvq) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/2 *
!!******************************
!else if (tvq==42) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==3 )*( x(k0, xcolvq) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/3 *
!!******************************
!else if (tvq==43) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==3 )*( x(k0, xcolvq) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/4 *
!!******************************
!else if (tvq==44) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==3 )*( x(k0, xcolvq) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 3/5 *
!!******************************
!else if (tvq==45) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==3 )*( x(k0, xcolvq) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/1 *
!!******************************
!else if (tvq==51) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==4 )*( x(k0, xcolvq) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/2 *
!!******************************
!else if (tvq==52) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==4 )*( x(k0, xcolvq) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/3 *
!!******************************
!else if (tvq==53) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==4 )*( x(k0, xcolvq) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/4 *
!!******************************
!else if (tvq==54) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==4 )*( x(k0, xcolvq) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 4/5 *
!!******************************
!else if (tvq==55) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==4 )*( x(k0, xcolvq) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/1 *
!!******************************
!else if (tvq==61) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==5 )*( x(k0, xcolvq) ==1 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/2 *
!!******************************
!else if (tvq==62) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==5 )*( x(k0, xcolvq) ==2 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/3 *
!!******************************
!else if (tvq==63) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==5 )*( x(k0, xcolvq) ==3 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/4 *
!!******************************
!else if (tvq==64) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==5 )*( x(k0, xcolvq) ==4 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!!******************************
!!* Differential Homophily 5/5 *
!!******************************
!else if (tvq==65) then          
!    do k0=1,nv
!        if (k0 .ne. iv) then
!            if (g(jv,k0)==1) then
!                dtfvq = dtfvq + ( x(iv, xcolvq)==5 )*( x(k0, xcolvq) ==5 )
!            endif !(g(jv,k0)==1)
!        endif !(k0 .ne. iv)
!    enddo !k0
!
!*********************
!*    cyclic triangles   *
!*********************
elseif (tvq==100) then
    do k0=1,nv
        if (k0 .ne. iv) then
            dtfvq = dtfvq + 1.0*g(jv,k0)*g(k0,iv)/3.0
        endif !(k0 .ne. iv)
 
    enddo !k0
	
else
	dtfvq=0.0
endif




!dtfvq = dtfvq/(nv*(nv-1)*(nv-2))
!dtfvq = dtfvq/(nv*(nv-1))
!dtfvq = dtfvq/10000.0

return
end function dtfvq




!************************************
!* subroutine: sample_ast_rc
!* version: 1.0
!* date: 11/09/2012
!* last modified:  06/26/2014
!* author: Angelo Mele
!* description: Samples using Accelerated Simulated
!*              Tempering as in Li et al (2004).
!               introduces steps where row or column
!               are modified
!*******************************************

!subroutine sample_ast_rc(nn,qq,pp,g,x,aa,dtt,skip,sample,&
!    &mtemp,qup,tout,kout,seed1,pr,pc)
subroutine sample_ast_rc(nn,qq,pp,g,x,aa,dtt,skip,sample,&
    &mtemp,qup,tout,seed1,pr,pc)

! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! gsim   = network matrix for simulations
! x      = matrix of exogenous characteristics
! aa     = vector of parameters for simulations
! dtt    = vector of statistics for simulations
! skip   = number of samples to skip
! sample = number of networks to sample
! mtemp  = number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g

! output variables
! tout   = vector of sampled network stats
! kout   = vector of temperatures sampled
!***********************************************************

implicit none

! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::skip
integer, intent(in)::sample
real(8), intent(in)::qup  
integer, intent(in)::mtemp 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc


! output vbls and pars
!real(8), dimension(pp(3),sample), intent(inout)::tout
real(8), dimension(sample, pp(3)), intent(inout)::tout
!real(8), dimension(sample), intent(inout)::kout

! working variables
integer, dimension(12):: seedn
integer, dimension(nn,nn)::gsim
integer i,j,k,rr,t, kk  ! indicators for loops
real(8), dimension(pp(3))::tsim, tsim2 ! current simulated network stats
integer ktemp !indicator of current temperature
real(8) invp(mtemp) ! probability of inversion for each temperature
real(8) temper(mtemp) ! vector with temperature levels
real(8) rhotemp(mtemp) ! vector of pseudopriors
real(8) rho0,nn0 ! constants for stochastic approximation of pseudopriors (see Geyer-Thompson 1995 for details)
integer totsim  ! maximum nummber of simulations
real(8) ux,uux ! random number
real(8) odds ! odds of two networks
real(8) dtfu, dtfm, dtfv ! network statistics functions
real(8) dttemp(pp(3))


! set prob of inversion
invp(1) = 0.0 ! .01 !.01 
do k=2,mtemp
invp(k)=invp(k-1)+.05
enddo

! set temperatures
temper(1)=1
do k=2,mtemp
	temper(k)=1 + .1*(k-1)
enddo

! set pseudopriors
rhotemp(:)=100000.0
rho0 = 1.0
nn0 = 10000000.0



! set initial simulated network the one provided inline
gsim = g

! compute observed network stats
tsim(:) = 0.0
call sstat(pp,nn,qq,gsim,x,dtt,tsim)

!tsim = tobs


!** random numbers seed
CALL RANDOM_SEED(size =  k)
seedn = seed1
CALL RANDOM_SEED(put = seedn(1:k) )



!**************************************************!
!******* MAIN SIMULATION LOOP STARTS HERE *********!
!**************************************************!

tout(:,:)=0.0  ! initialize 
rr=1   ! sample loop indicator
ktemp =1  ! initialize temperature
totsim = skip*sample ! 10000000 ! max number of simulations

! main loop
do t=1, totsim
    
    call random_number(uux)
    !print*, uux
    ! row step
    if (uux<pr) then
        !print*, 'row'
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose row/player i
        gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
		gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii

        endif ! (log(ux)<=min(0.0,odds))
    
    ! column step    
    elseif (pr <= uux .and. uux< pr+pc) then
        !print*, 'column' 
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose column i
        gsim(:,i) = 1 - gsim(:,i) !invert all links of column i
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
			   !print*, 'accepted'
		else
			gsim(:,i) = 1 - gsim(:,i) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
		
        endif ! (log(ux)<=min(0.0,odds))
        
   
    ! one link step
    !elseif (uux <= pr+pc) then
     else   
        !print*, 'link'
        
        
        ! one-link per iteration 
        call random_number(ux)
        i = ceiling(nn*ux) ! choose player i
        call random_number(ux)
        j = ceiling(nn*ux) ! choose agent j
        
        do while (i==j)
            call random_number(ux)
            j = ceiling(nn*ux) ! choose agent j
        enddo ! while loop
        
        call random_number(ux)
        odds = 0.0
        dttemp(:) = 0.0 
        
        ! direct utility update
        do kk=1,pp(1)
            dttemp(kk) = dtfu(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        ! mutual utility update        
        if (gsim(j,i)==1) then  ! g_ji=1 then compute mutual utility
            do kk=pp(1)+1,pp(2)
                dttemp(kk) = dtfm(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2))
            enddo
        endif
        ! indirect utility and popularity update
        do kk=pp(2)+1,pp(3)
            dttemp(kk) =  dtfv(nn,qq,gsim,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        
       
        
        ! compute odds ration from the model
        do kk=1,pp(3)
            odds = odds + aa(kk)*dttemp(kk)
        enddo
        !odds = odds/temper(ktemp)
        
    
    
        ! update link g_ij
        if (gsim(i,j) ==1 ) then
            if (odds<0) then
                gsim(i,j)=0 
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) - dttemp(kk)
                enddo
                
            else
                if (log(ux)<-odds) then
                   gsim(i,j)=0 
                   do kk=1,pp(3)
                        tsim(kk) = tsim(kk) - dttemp(kk)
                    enddo
                endif ! (log(ux)<-odds)
                
            endif ! (odds<0)
            
        elseif (gsim(i,j)==0) then
            if (odds>=0) then
                gsim(i,j)=1
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) + dttemp(kk)
                enddo
            else
                if (log(ux)<odds) then
                    gsim(i,j)=1
                    do kk=1,pp(3)
                        tsim(kk) = tsim(kk) + dttemp(kk)
                    enddo
                    
                endif !(log(ux)<=odds)
                
            endif !(odds>=0)
        endif !(gsim(i,j) ==1 )
          
    endif !(ux<invp(ktemp))

    
    ! sample every skip iterations
    if (rr<=sample) then
        if (t==skip*rr) then
            !tout(:,rr)=tsim
            tout(rr,:)=tsim
            !kout(rr)=ktemp
            !if (mod(rr,1000)==0) then
            !    print*, 'sample = ', rr
            !    print*, 'tsim = ', tsim
            !endif
            rr=rr+1
        endif
    else
        exit
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! introduce simulated tempering
    
enddo ! t

end subroutine sample_ast_rc



    
    

!************************************
!* subroutine: sample_ast_one
!* version: 1.0
!* date: 11/09/2012
!* last modified:  06/26/2014
!* author: Angelo Mele
!* description: Samples using Accelerated Simulated
!*              Tempering as in Li et al (2004).
!               introduces steps where row or column
!               are modified
!*******************************************

subroutine sample_ast_one(nn,qq,pp,g,x,aa,dtt,sample,&
    &mtemp,qup,tout,seed1,pr,pc)

! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! gsim   = network matrix for simulations
! x      = matrix of exogenous characteristics
! aa     = vector of parameters for simulations
! dtt    = vector of statistics for simulations
! sample = number of networks to sample
! mtemp  = number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g

! output variables
! tout   = vector of sampled network stats
! kout   = vector of temperatures sampled
!***********************************************************

implicit none

! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::sample
real(8), intent(in)::qup  
integer, intent(in)::mtemp 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc


! output vbls and pars
!real(8), dimension(pp(3),sample), intent(inout)::tout
real(8), dimension(pp(3)), intent(inout)::tout
!real(8), dimension(sample), intent(inout)::kout

! working variables
integer, dimension(12):: seedn
integer, dimension(nn,nn)::gsim
integer i,j,k,rr,t, kk  ! indicators for loops
real(8), dimension(pp(3))::tsim, tsim2 ! current simulated network stats
integer ktemp !indicator of current temperature
real(8) invp(mtemp) ! probability of inversion for each temperature
real(8) temper(mtemp) ! vector with temperature levels
real(8) rhotemp(mtemp) ! vector of pseudopriors
real(8) rho0,nn0 ! constants for stochastic approximation of pseudopriors (see Geyer-Thompson 1995 for details)
integer totsim  ! maximum nummber of simulations
real(8) ux,uux ! random number
real(8) odds ! odds of two networks
real(8) dtfu, dtfm, dtfv ! network statistics functions
real(8) dttemp(pp(3))


! set prob of inversion
invp(1) = 0.0 ! .01 !.01 
do k=2,mtemp
invp(k)=invp(k-1)+.05
enddo

! set temperatures
temper(1)=1
do k=2,mtemp
	temper(k)=1 + .1*(k-1)
enddo

! set pseudopriors
rhotemp(:)=100000.0
rho0 = 1.0
nn0 = 10000000.0



! set initial simulated network the one provided inline
gsim = g

! compute observed network stats
tsim(:) = 0.0
call sstat(pp,nn,qq,gsim,x,dtt,tsim)

!tsim = tobs


!** random numbers seed
CALL RANDOM_SEED(size =  k)
seedn = seed1
CALL RANDOM_SEED(put = seedn(1:k) )



!**************************************************!
!******* MAIN SIMULATION LOOP STARTS HERE *********!
!**************************************************!

tout(:)=0.0  ! initialize 
rr=1   ! sample loop indicator
ktemp =1  ! initialize temperature
totsim = sample ! 10000000 ! max number of simulations

! main loop
do t=1, totsim
    
    call random_number(uux)
    !print*, uux
    ! row step
    if (uux<pr) then
        !print*, 'row'
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose row/player i
        gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
		gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii

        endif ! (log(ux)<=min(0.0,odds))
    
    ! column step    
    elseif (pr <= uux .and. uux< pr+pc) then
        !print*, 'column' 
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose column i
        gsim(:,i) = 1 - gsim(:,i) !invert all links of column i
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
			   !print*, 'accepted'
		else
			gsim(:,i) = 1 - gsim(:,i) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
		
        endif ! (log(ux)<=min(0.0,odds))
        
   
    ! one link step
    !elseif (uux <= pr+pc) then
     else   
        !print*, 'link'
        
        
        ! one-link per iteration 
        call random_number(ux)
        i = ceiling(nn*ux) ! choose player i
        call random_number(ux)
        j = ceiling(nn*ux) ! choose agent j
        
        do while (i==j)
            call random_number(ux)
            j = ceiling(nn*ux) ! choose agent j
        enddo ! while loop
        
        call random_number(ux)
        odds = 0.0
        dttemp(:) = 0.0 
        
        ! direct utility update
        do kk=1,pp(1)
            dttemp(kk) = dtfu(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        ! mutual utility update        
        if (gsim(j,i)==1) then  ! g_ji=1 then compute mutual utility
            do kk=pp(1)+1,pp(2)
                dttemp(kk) = dtfm(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2))
            enddo
        endif
        ! indirect utility and popularity update
        do kk=pp(2)+1,pp(3)
            dttemp(kk) =  dtfv(nn,qq,gsim,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        
       
        
        ! compute odds ration from the model
        do kk=1,pp(3)
            odds = odds + aa(kk)*dttemp(kk)
        enddo
        !odds = odds/temper(ktemp)
        
    
    
        ! update link g_ij
        if (gsim(i,j) ==1 ) then
            if (odds<0) then
                gsim(i,j)=0 
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) - dttemp(kk)
                enddo
                
            else
                if (log(ux)<-odds) then
                   gsim(i,j)=0 
                   do kk=1,pp(3)
                        tsim(kk) = tsim(kk) - dttemp(kk)
                    enddo
                endif ! (log(ux)<-odds)
                
            endif ! (odds<0)
            
        elseif (gsim(i,j)==0) then
            if (odds>=0) then
                gsim(i,j)=1
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) + dttemp(kk)
                enddo
            else
                if (log(ux)<odds) then
                    gsim(i,j)=1
                    do kk=1,pp(3)
                        tsim(kk) = tsim(kk) + dttemp(kk)
                    enddo
                    
                endif !(log(ux)<=odds)
                
            endif !(odds>=0)
        endif !(gsim(i,j) ==1 )
          
    endif !(ux<invp(ktemp))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! introduce simulated tempering
    
enddo ! t


tout=tsim
    

    end subroutine sample_ast_one

    

    
    


!************************************
!* subroutine: sample_ast_large
!* version: 1.0
!* date: 07/13/2014
!* last modified: 07/13/2014
!* author: Angelo Mele
!* description: Samples using Accelerated Simulated
!*              Tempering as in Li et al (2004).
!               introduces steps where row or column
!               are modified, and larger steps that are not o(n)
!*******************************************

!subroutine sample_ast_rc(nn,qq,pp,g,x,aa,dtt,skip,sample,&
!    &mtemp,qup,tout,kout,seed1,pr,pc)
subroutine sample_ast_large(nn, qq, pp, g, x, aa, dtt, skip, sample,&
    &mtemp, qup, tout, seed1, pr, pc, pf, size)

! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! gsim   = network matrix for simulations
! x      = matrix of exogenous characteristics
! aa     = vector of parameters for simulations
! dtt    = vector of statistics for simulations
! skip   = number of samples to skip
! sample = number of networks to sample
! mtemp  = number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g
! pf     = prob of large step
! size   = number of links updated in each iteration for large step

! output variables
! tout   = vector of sampled network stats
! kout   = vector of temperatures sampled
!***********************************************************

implicit none

! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::skip
integer, intent(in)::sample
real(8), intent(in)::qup  
integer, intent(in)::mtemp 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc
real(8), intent(in)::pf
integer, intent(in)::size


! output vbls and pars
!real(8), dimension(pp(3),sample), intent(inout)::tout
real(8), dimension(sample, pp(3)), intent(inout)::tout
!real(8), dimension(sample), intent(inout)::kout

! working variables
integer, dimension(12):: seedn
integer, dimension(nn,nn)::gsim
integer i,j,k,rr,t, kk  ! indicators for loops
integer j0, k0, s0 ! indicators for loops
integer trand, tinsamp, diagcheck ! variables for large steps
integer nsq, posrow, poscol
integer, dimension(size)::samplinks
integer, dimension(size,2)::updlinks
real(8), dimension(pp(3))::tsim, tsim2 ! current simulated network stats
integer ktemp !indicator of current temperature
real(8) invp(mtemp) ! probability of inversion for each temperature
real(8) temper(mtemp) ! vector with temperature levels
real(8) rhotemp(mtemp) ! vector of pseudopriors
real(8) rho0,nn0 ! constants for stochastic approximation of pseudopriors (see Geyer-Thompson 1995 for details)
integer totsim  ! maximum nummber of simulations
real(8) ux,uux ! random number
real(8) odds ! odds of two networks
real(8) dtfu, dtfm, dtfv ! network statistics functions
real(8) dttemp(pp(3))


nsq = nn*2

! set prob of inversion
invp(1) = 0.0 ! .01 !.01 
do k=2,mtemp
invp(k)=invp(k-1)+.05
enddo

! set temperatures
temper(1)=1
do k=2,mtemp
	temper(k)=1 + .1*(k-1)
enddo

! set pseudopriors
rhotemp(:)=100000.0
rho0 = 1.0
nn0 = 10000000.0


! set initial simulated network the one provided inline
gsim = g

! compute observed network stats
tsim(:) = 0.0
call sstat(pp,nn,qq,gsim,x,dtt,tsim)

!tsim = tobs


!** random numbers seed
CALL RANDOM_SEED(size =  k)
seedn = seed1
CALL RANDOM_SEED(put = seedn(1:k) )



!**************************************************!
!******* MAIN SIMULATION LOOP STARTS HERE *********!
!**************************************************!

tout(:,:)=0.0  ! initialize 
rr=1   ! sample loop indicator
ktemp =1  ! initialize temperature
totsim = skip*sample ! 10000000 ! max number of simulations

!***************** main loop
do t=1, totsim
    
    call random_number(uux)
    !print*, uux
    
    
    !******************* row step
    
    if (uux<pr) then
        !print*, 'row'
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose row/player i
        gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
		gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii

        endif ! (log(ux)<=min(0.0,odds))
    
    
    !***************************column step   
    
    elseif (pr <= uux .and. uux< pr+pc) then
        !print*, 'column' 
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose column i
        gsim(:,i) = 1 - gsim(:,i) !invert all links of column i
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
			   !print*, 'accepted'
		else
			gsim(:,i) = 1 - gsim(:,i) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
		
        endif ! (log(ux)<=min(0.0,odds))
     
        
        
    !******************************* large step    
    
    elseif (pr+pc <= uux .and. uux <= pr+pc+pf) then
            ! propose a uniform random sample of links
            ! to update (size = number of links)
            k0=1
            do while (k0<=size)
              call random_number(ux)
              trand = ceiling(nsq*ux)
              !trand = 1 + int(nsq*ux)
              tinsamp = 0
              do s0=1,k0-1
                  if (samplinks(s0)==trand) then
                      tinsamp = 1
                      
                  endif
              enddo !s0
              
              if (tinsamp==0) then
                  posrow = ceiling(1.0*trand/nn)
                  poscol = trand - (posrow-1)*nn
                  if (posrow .ne. poscol) then
                      samplinks(k0) = trand
                      k0=k0+1
                      !print*, 'accepted, k0 = ', k0-1
                  endif
              endif !(tinsamp==0)  
            enddo !(k0<=size)
            
            do s0=1,size
                updlinks(s0,1) = ceiling(1.0*samplinks(s0)/nn)
                updlinks(s0,2) = samplinks(s0) - (updlinks(s0,1)-1)*nn
            enddo !k0
            
            do s0=1,size
                gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
            enddo
            call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
            odds = 0.0
            do kk=1,pp(3)
                odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
            enddo ! kk
            call random_number(ux)
            if (log(ux)<=min(0.0,odds)) then
                tsim = tsim2
                !print*, 'accepted'
		    else
                do s0=1,size
                    gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                    &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
                enddo
            endif ! (log(ux)<=min(0.0,odds))
        
    !********************************* one link step
    
    !elseif (uux <= pr+pc) then
     else   
        !print*, 'link'
        
        
        ! one-link per iteration 
        call random_number(ux)
        i = ceiling(nn*ux) ! choose player i
        call random_number(ux)
        j = ceiling(nn*ux) ! choose agent j
        
        do while (i==j)
            call random_number(ux)
            j = ceiling(nn*ux) ! choose agent j
        enddo ! while loop
        
        call random_number(ux)
        odds = 0.0
        dttemp(:) = 0.0 
        
        ! direct utility update
        do kk=1,pp(1)
            dttemp(kk) = dtfu(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        ! mutual utility update        
        if (gsim(j,i)==1) then  ! g_ji=1 then compute mutual utility
            do kk=pp(1)+1,pp(2)
                dttemp(kk) = dtfm(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2))
            enddo
        endif
        ! indirect utility and popularity update
        do kk=pp(2)+1,pp(3)
            dttemp(kk) =  dtfv(nn,qq,gsim,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        
       
        
        ! compute odds ration from the model
        do kk=1,pp(3)
            odds = odds + aa(kk)*dttemp(kk)
        enddo
        !odds = odds/temper(ktemp)
        
    
    
        ! update link g_ij
        if (gsim(i,j) ==1 ) then
            if (odds<0) then
                gsim(i,j)=0 
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) - dttemp(kk)
                enddo
                
            else
                if (log(ux)<-odds) then
                   gsim(i,j)=0 
                   do kk=1,pp(3)
                        tsim(kk) = tsim(kk) - dttemp(kk)
                    enddo
                endif ! (log(ux)<-odds)
                
            endif ! (odds<0)
            
        elseif (gsim(i,j)==0) then
            if (odds>=0) then
                gsim(i,j)=1
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) + dttemp(kk)
                enddo
            else
                if (log(ux)<odds) then
                    gsim(i,j)=1
                    do kk=1,pp(3)
                        tsim(kk) = tsim(kk) + dttemp(kk)
                    enddo
                    
                endif !(log(ux)<=odds)
                
            endif !(odds>=0)
        endif !(gsim(i,j) ==1 )
          
    endif !(ux<invp(ktemp))

    
    ! sample every skip iterations
    if (rr<=sample) then
        if (t==skip*rr) then
            !tout(:,rr)=tsim
            tout(rr,:)=tsim
            !kout(rr)=ktemp
            !if (mod(rr,1000)==0) then
            !    print*, 'sample = ', rr
            !    print*, 'tsim = ', tsim
            !endif
            rr=rr+1
        endif
    else
        exit
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! introduce simulated tempering
    
enddo ! t

    end subroutine sample_ast_large



    
    
    
    
    
!************************************
!* subroutine: sample_ast_bham
!* version: 1.0
!* date: 12/23/2014
!* last modified: 12/23/2014
!* author: Angelo Mele
!* description: Samples using Bhamidi et al (2011)
!*              mean-field approximations
!*******************************************

subroutine sample_ast_bham(nn, qq, pp, g, x, aa, dtt, skip, sample,&
    &mtemp, qup, tout, seed1, pr, pc, pf, pb, pinv, nmus, mus, size)

! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! gsim   = network matrix for simulations
! x      = matrix of exogenous characteristics
! aa     = vector of parameters for simulations
! dtt    = vector of statistics for simulations
! skip   = number of samples to skip
! sample = number of networks to sample
! mtemp  = number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g
! pf     = prob of large step
! pb     = prob of using Bhamidi-like step
! pinv   = prob of matrix inversion
! nmus   = number of modes of the gibbs distribution
! mus    = modes of the gibbs distribution
! size   = number of links updated in each iteration for large step

! output variables
! tout   = vector of sampled network stats
! kout   = vector of temperatures sampled
!***********************************************************

implicit none

! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::skip
integer, intent(in)::sample
real(8), intent(in)::qup  
integer, intent(in)::mtemp 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc
real(8), intent(in)::pf
real(8), intent(in)::pb
real(8), intent(in)::pinv
integer, intent(in)::nmus
real(8), intent(in), dimension(nmus)::mus
integer, intent(in)::size


! output vbls and pars
!real(8), dimension(pp(3),sample), intent(inout)::tout
real(8), dimension(sample, pp(3)), intent(inout)::tout
!real(8), dimension(sample), intent(inout)::kout

! working variables
integer, dimension(12):: seedn
integer, dimension(nn,nn)::gsim, gsim2
integer i,j,k,rr,t, kk  ! indicators for loops
integer j0, k0, s0 ! indicators for loops
integer trand, tinsamp, diagcheck ! variables for large steps
integer nsq, posrow, poscol
integer, dimension(size)::samplinks
integer, dimension(size,2)::updlinks
real(8), dimension(pp(3))::tsim, tsim2 ! current simulated network stats
integer ktemp !indicator of current temperature
real(8) invp(mtemp) ! probability of inversion for each temperature
real(8) temper(mtemp) ! vector with temperature levels
real(8) rhotemp(mtemp) ! vector of pseudopriors
real(8) rho0,nn0 ! constants for stochastic approximation of pseudopriors (see Geyer-Thompson 1995 for details)
integer totsim  ! maximum nummber of simulations
real(8) ux,uux ! random number
real(8) odds ! odds of two networks
real(8) dtfu, dtfm, dtfv ! network statistics functions
real(8) dttemp(pp(3))
integer propmu ! mode proposed for update
integer curmu  ! current mode

curmu = 1

nsq = nn*2

! set prob of inversion
invp(1) = 0.0 ! .01 !.01 
do k=2,mtemp
invp(k)=invp(k-1)+.05
enddo

! set temperatures
temper(1)=1
do k=2,mtemp
	temper(k)=1 + .1*(k-1)
enddo

! set pseudopriors
rhotemp(:)=100000.0
rho0 = 1.0
nn0 = 10000000.0


! set initial simulated network the one provided inline
gsim = g

! compute observed network stats
tsim(:) = 0.0
call sstat(pp,nn,qq,gsim,x,dtt,tsim)

!tsim = tobs


!** random numbers seed
CALL RANDOM_SEED(size =  k)
seedn = seed1
CALL RANDOM_SEED(put = seedn(1:k) )



!**************************************************!
!******* MAIN SIMULATION LOOP STARTS HERE *********!
!**************************************************!

tout(:,:)=0.0  ! initialize 
rr=1   ! sample loop indicator
ktemp =1  ! initialize temperature
totsim = skip*sample ! 10000000 ! max number of simulations

!***************** main loop
do t=1, totsim
    
    call random_number(uux)
    !print*, uux
    
    
    !******************* row step
    
    if (uux<pr) then
        !print*, 'row'
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose row/player i
        gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
		gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii

        endif ! (log(ux)<=min(0.0,odds))
    
    
    !***************************column step   
    
    elseif (pr <= uux .and. uux< pr+pc) then
        !print*, 'column' 
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose column i
        gsim(:,i) = 1 - gsim(:,i) !invert all links of column i
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
			   !print*, 'accepted'
		else
			gsim(:,i) = 1 - gsim(:,i) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
		
        endif ! (log(ux)<=min(0.0,odds))
     
        
        
    !******************************* large step    
    
    elseif (pr+pc <= uux .and. uux < pr+pc+pf) then
            ! propose a uniform random sample of links
            ! to update (size = number of links)
            k0=1
            do while (k0<=size)
              call random_number(ux)
              trand = ceiling(nsq*ux)
              !trand = 1 + int(nsq*ux)
              tinsamp = 0
              do s0=1,k0-1
                  if (samplinks(s0)==trand) then
                      tinsamp = 1
                      
                  endif
              enddo !s0
              
              if (tinsamp==0) then
                  posrow = ceiling(1.0*trand/nn)
                  poscol = trand - (posrow-1)*nn
                  if (posrow .ne. poscol) then
                      samplinks(k0) = trand
                      k0=k0+1
                      !print*, 'accepted, k0 = ', k0-1
                  endif
              endif !(tinsamp==0)  
            enddo !(k0<=size)
            
            do s0=1,size
                updlinks(s0,1) = ceiling(1.0*samplinks(s0)/nn)
                updlinks(s0,2) = samplinks(s0) - (updlinks(s0,1)-1)*nn
            enddo !k0
            
            do s0=1,size
                gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
            enddo
            call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
            odds = 0.0
            do kk=1,pp(3)
                odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
            enddo ! kk
            call random_number(ux)
            if (log(ux)<=min(0.0,odds)) then
                tsim = tsim2
                !print*, 'accepted'
		    else
                do s0=1,size
                    gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                    &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
                enddo
            endif ! (log(ux)<=min(0.0,odds))
    

    !******************* bhamidi-like step
    
    elseif (pr+pc+pf <= uux .and. uux < pr+pc+pf+pb) then
        ! select which mode to propose
        !propmu = 1  ! if only 1 mode
        !if (nmus>1) then
        !    call random_number(ux)
        !    if (ux<0.5) then
        !        propmu = 1
        !    else
        !        propmu = 2
        !    endif !(ux<0.5)
        !endif !(nmus>1) 
        if (curmu == 1) then
            propmu = 2
        else
            propmu = 1
        endif
        
        ! new network is directed Erdos-Renyi 
        ! with density mus(propmu)
        gsim2 = 0
        do i=1,nn
            do j=1,nn
                if (i /= j) then
                    call random_number(ux)
                    if (ux < mus(propmu) ) then
                        gsim2(i,j) = 1
                    endif !(ux < mus(propmu) )
                endif 
            enddo !j 
        enddo  !i 

        ! compute network stats for proposed network
        call sstat(pp,nn,qq,gsim2,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            gsim = gsim2
            tsim = tsim2
            curmu = propmu
            !print*, 'accepted'
        endif ! (log(ux)<=min(0.0,odds))
    

    !******************* inversion step
    
    elseif (pr+pc+pf+pb <= uux .and. uux < pr+pc+pf+pb+pinv) then
        !print*, 'row'
        
        gsim = 1 - gsim !invert all links
		do i=1,nn
            gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        enddo 
        
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim = 1 - gsim !invert all links
    		do i=1,nn
                gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
            enddo 
        endif ! (log(ux)<=min(0.0,odds))
    
        
    !********************************* one link step
    
    !elseif (uux <= pr+pc) then
     else   
        !print*, 'link'
        
        
        ! one-link per iteration 
        call random_number(ux)
        i = ceiling(nn*ux) ! choose player i
        call random_number(ux)
        j = ceiling(nn*ux) ! choose agent j
        
        do while (i==j)
            call random_number(ux)
            j = ceiling(nn*ux) ! choose agent j
        enddo ! while loop
        
        call random_number(ux)
        odds = 0.0
        dttemp(:) = 0.0 
        
        ! direct utility update
        do kk=1,pp(1)
            dttemp(kk) = dtfu(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        ! mutual utility update        
        if (gsim(j,i)==1) then  ! g_ji=1 then compute mutual utility
            do kk=pp(1)+1,pp(2)
                dttemp(kk) = dtfm(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2))
            enddo
        endif
        ! indirect utility and popularity update
        do kk=pp(2)+1,pp(3)
            dttemp(kk) =  dtfv(nn,qq,gsim,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        
       
        
        ! compute odds ration from the model
        do kk=1,pp(3)
            odds = odds + aa(kk)*dttemp(kk)
        enddo
        !odds = odds/temper(ktemp)
        
    
    
        ! update link g_ij
        if (gsim(i,j) ==1 ) then
            if (odds<0) then
                gsim(i,j)=0 
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) - dttemp(kk)
                enddo
                
            else
                if (log(ux)<-odds) then
                   gsim(i,j)=0 
                   do kk=1,pp(3)
                        tsim(kk) = tsim(kk) - dttemp(kk)
                    enddo
                endif ! (log(ux)<-odds)
                
            endif ! (odds<0)
            
        elseif (gsim(i,j)==0) then
            if (odds>=0) then
                gsim(i,j)=1
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) + dttemp(kk)
                enddo
            else
                if (log(ux)<odds) then
                    gsim(i,j)=1
                    do kk=1,pp(3)
                        tsim(kk) = tsim(kk) + dttemp(kk)
                    enddo
                    
                endif !(log(ux)<=odds)
                
            endif !(odds>=0)
        endif !(gsim(i,j) ==1 )
          
    endif !(ux<invp(ktemp))

    
    ! sample every skip iterations
    if (rr<=sample) then
        if (t==skip*rr) then
            !tout(:,rr)=tsim
            tout(rr,:)=tsim
            !kout(rr)=ktemp
            !if (mod(rr,1000)==0) then
            !    print*, 'sample = ', rr
            !    print*, 'tsim = ', tsim
            !endif
            rr=rr+1
        endif
    else
        exit
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! introduce simulated tempering
    
enddo ! t

    end subroutine sample_ast_bham



    
    
    
    
    !************************************
!* subroutine: sample_ast_bhamg
!* version: 1.0
!* date: 12/23/2014
!* last modified: 12/23/2014
!* author: Angelo Mele
!* description: Samples using Bhamidi et al (2011)
!*              mean-field approximations
!*******************************************

subroutine sample_ast_bhamg(nn, qq, pp, g, x, aa, dtt, skip, sample,&
    &mtemp, qup, gout, seed1, pr, pc, pf, pb, pinv, nmus, mus, size)

! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! gsim   = network matrix for simulations
! x      = matrix of exogenous characteristics
! aa     = vector of parameters for simulations
! dtt    = vector of statistics for simulations
! skip   = number of samples to skip
! sample = number of networks to sample
! mtemp  = number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g
! pf     = prob of large step
! pb     = prob of using Bhamidi-like step
! pinv   = prob of matrix inversion
! nmus   = number of modes of the gibbs distribution
! mus    = modes of the gibbs distribution
! size   = number of links updated in each iteration for large step


!***********************************************************

implicit none

! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::skip
integer, intent(in)::sample
real(8), intent(in)::qup  
integer, intent(in)::mtemp 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc
real(8), intent(in)::pf
real(8), intent(in)::pb
real(8), intent(in)::pinv
integer, intent(in)::nmus
real(8), intent(in), dimension(nmus)::mus
integer, intent(in)::size


! output vbls and pars
!real(8), dimension(pp(3),sample), intent(inout)::tout
!real(8), dimension(sample, pp(3)), intent(inout)::tout
!real(8), dimension(sample), intent(inout)::kout
integer, intent(out), dimension(nn,nn)::gout


! working variables
integer, dimension(12):: seedn
integer, dimension(nn,nn)::gsim, gsim2
integer i,j,k,rr,t, kk  ! indicators for loops
integer j0, k0, s0 ! indicators for loops
integer trand, tinsamp, diagcheck ! variables for large steps
integer nsq, posrow, poscol
integer, dimension(size)::samplinks
integer, dimension(size,2)::updlinks
real(8), dimension(pp(3))::tsim, tsim2 ! current simulated network stats
integer ktemp !indicator of current temperature
real(8) invp(mtemp) ! probability of inversion for each temperature
real(8) temper(mtemp) ! vector with temperature levels
real(8) rhotemp(mtemp) ! vector of pseudopriors
real(8) rho0,nn0 ! constants for stochastic approximation of pseudopriors (see Geyer-Thompson 1995 for details)
integer totsim  ! maximum nummber of simulations
real(8) ux,uux ! random number
real(8) odds ! odds of two networks
real(8) dtfu, dtfm, dtfv ! network statistics functions
real(8) dttemp(pp(3))
integer propmu ! mode proposed for update
integer curmu  ! current mode

curmu = 1

nsq = nn*2

! set prob of inversion
invp(1) = 0.0 ! .01 !.01 
do k=2,mtemp
invp(k)=invp(k-1)+.05
enddo

! set temperatures
temper(1)=1
do k=2,mtemp
	temper(k)=1 + .1*(k-1)
enddo

! set pseudopriors
rhotemp(:)=100000.0
rho0 = 1.0
nn0 = 10000000.0


! set initial simulated network the one provided inline
gsim = g
gout = gsim

! compute observed network stats
tsim(:) = 0.0
call sstat(pp,nn,qq,gsim,x,dtt,tsim)

!tsim = tobs


!** random numbers seed
CALL RANDOM_SEED(size =  k)
seedn = seed1
CALL RANDOM_SEED(put = seedn(1:k) )



!**************************************************!
!******* MAIN SIMULATION LOOP STARTS HERE *********!
!**************************************************!

!tout(:,:)=0.0  ! initialize 
rr=1   ! sample loop indicator
ktemp =1  ! initialize temperature
totsim = skip*sample ! 10000000 ! max number of simulations

!***************** main loop
do t=1, totsim
    
    call random_number(uux)
    !print*, uux
    
    
    !******************* row step
    
    if (uux<pr) then
        !print*, 'row'
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose row/player i
        gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
		gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii

        endif ! (log(ux)<=min(0.0,odds))
    
    
    !***************************column step   
    
    elseif (pr <= uux .and. uux< pr+pc) then
        !print*, 'column' 
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose column i
        gsim(:,i) = 1 - gsim(:,i) !invert all links of column i
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
			   !print*, 'accepted'
		else
			gsim(:,i) = 1 - gsim(:,i) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
		
        endif ! (log(ux)<=min(0.0,odds))
     
        
        
    !******************************* large step    
    
    elseif (pr+pc <= uux .and. uux < pr+pc+pf) then
            ! propose a uniform random sample of links
            ! to update (size = number of links)
            k0=1
            do while (k0<=size)
              call random_number(ux)
              trand = ceiling(nsq*ux)
              !trand = 1 + int(nsq*ux)
              tinsamp = 0
              do s0=1,k0-1
                  if (samplinks(s0)==trand) then
                      tinsamp = 1
                      
                  endif
              enddo !s0
              
              if (tinsamp==0) then
                  posrow = ceiling(1.0*trand/nn)
                  poscol = trand - (posrow-1)*nn
                  if (posrow .ne. poscol) then
                      samplinks(k0) = trand
                      k0=k0+1
                      !print*, 'accepted, k0 = ', k0-1
                  endif
              endif !(tinsamp==0)  
            enddo !(k0<=size)
            
            do s0=1,size
                updlinks(s0,1) = ceiling(1.0*samplinks(s0)/nn)
                updlinks(s0,2) = samplinks(s0) - (updlinks(s0,1)-1)*nn
            enddo !k0
            
            do s0=1,size
                gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
            enddo
            call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
            odds = 0.0
            do kk=1,pp(3)
                odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
            enddo ! kk
            call random_number(ux)
            if (log(ux)<=min(0.0,odds)) then
                tsim = tsim2
                !print*, 'accepted'
		    else
                do s0=1,size
                    gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                    &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
                enddo
            endif ! (log(ux)<=min(0.0,odds))
    

    !******************* bhamidi-like step
    
    elseif (pr+pc+pf <= uux .and. uux < pr+pc+pf+pb) then
        ! select which mode to propose
        !propmu = 1  ! if only 1 mode
        !if (nmus>1) then
        !    call random_number(ux)
        !    if (ux<0.5) then
        !        propmu = 1
        !    else
        !        propmu = 2
        !    endif !(ux<0.5)
        !endif !(nmus>1) 
        if (curmu == 1) then
            propmu = 2
        else
            propmu = 1
        endif
        
        ! new network is directed Erdos-Renyi 
        ! with density mus(propmu)
        gsim2 = 0
        do i=1,nn
            do j=1,nn
                if (i /= j) then
                    call random_number(ux)
                    if (ux < mus(propmu) ) then
                        gsim2(i,j) = 1
                    endif !(ux < mus(propmu) )
                endif 
            enddo !j 
        enddo  !i 

        ! compute network stats for proposed network
        call sstat(pp,nn,qq,gsim2,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            gsim = gsim2
            tsim = tsim2
            curmu = propmu
            !print*, 'accepted'
        endif ! (log(ux)<=min(0.0,odds))
    

    !******************* inversion step
    
    elseif (pr+pc+pf+pb <= uux .and. uux < pr+pc+pf+pb+pinv) then
        !print*, 'row'
        
        gsim = 1 - gsim !invert all links
		do i=1,nn
            gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        enddo 
        
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim = 1 - gsim !invert all links
    		do i=1,nn
                gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
            enddo 
        endif ! (log(ux)<=min(0.0,odds))
    
        
    !********************************* one link step
    
    !elseif (uux <= pr+pc) then
     else   
        !print*, 'link'
        
        
        ! one-link per iteration 
        call random_number(ux)
        i = ceiling(nn*ux) ! choose player i
        call random_number(ux)
        j = ceiling(nn*ux) ! choose agent j
        
        do while (i==j)
            call random_number(ux)
            j = ceiling(nn*ux) ! choose agent j
        enddo ! while loop
        
        call random_number(ux)
        odds = 0.0
        dttemp(:) = 0.0 
        
        ! direct utility update
        do kk=1,pp(1)
            dttemp(kk) = dtfu(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        ! mutual utility update        
        if (gsim(j,i)==1) then  ! g_ji=1 then compute mutual utility
            do kk=pp(1)+1,pp(2)
                dttemp(kk) = dtfm(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2))
            enddo
        endif
        ! indirect utility and popularity update
        do kk=pp(2)+1,pp(3)
            dttemp(kk) =  dtfv(nn,qq,gsim,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        
       
        
        ! compute odds ration from the model
        do kk=1,pp(3)
            odds = odds + aa(kk)*dttemp(kk)
        enddo
        !odds = odds/temper(ktemp)
        
    
    
        ! update link g_ij
        if (gsim(i,j) ==1 ) then
            if (odds<0) then
                gsim(i,j)=0 
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) - dttemp(kk)
                enddo
                
            else
                if (log(ux)<-odds) then
                   gsim(i,j)=0 
                   do kk=1,pp(3)
                        tsim(kk) = tsim(kk) - dttemp(kk)
                    enddo
                endif ! (log(ux)<-odds)
                
            endif ! (odds<0)
            
        elseif (gsim(i,j)==0) then
            if (odds>=0) then
                gsim(i,j)=1
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) + dttemp(kk)
                enddo
            else
                if (log(ux)<odds) then
                    gsim(i,j)=1
                    do kk=1,pp(3)
                        tsim(kk) = tsim(kk) + dttemp(kk)
                    enddo
                    
                endif !(log(ux)<=odds)
                
            endif !(odds>=0)
        endif !(gsim(i,j) ==1 )
          
    endif !(ux<invp(ktemp))

    
    ! copy last network
    if (rr<=totsim) then
        gout = gsim
    else
        exit
    endif    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! introduce simulated tempering
    
enddo ! t

    end subroutine sample_ast_bhamg



    
    
    
    
    
    
  
    
!************************************
!* subroutine: sample_ast_bham
!* version: 1.0
!* date: 12/23/2014
!* last modified: 12/23/2014
!* author: Angelo Mele
!* description: Samples using Bhamidi et al (2011)
!*              mean-field approximations
!*******************************************

subroutine sample_ast_bham_est(nn, qq, pp, g, x, aa, dtt, skip, sample,&
    &mtemp, qup, tout, seed1, pr, pc, pf, pb, pinv, nmus, mus, size)

! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! gsim   = network matrix for simulations
! x      = matrix of exogenous characteristics
! aa     = vector of parameters for simulations
! dtt    = vector of statistics for simulations
! skip   = number of samples to skip
! sample = number of networks to sample
! mtemp  = number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g
! pf     = prob of large step
! pb     = prob of using Bhamidi-like step
! pinv   = prob of matrix inversion
! nmus   = number of modes of the gibbs distribution
! mus    = modes of the gibbs distribution
! size   = number of links updated in each iteration for large step

! output variables
! tout   = vector of sampled network stats
! kout   = vector of temperatures sampled
!***********************************************************

implicit none

! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::skip
integer, intent(in)::sample
real(8), intent(in)::qup  
integer, intent(in)::mtemp 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc
real(8), intent(in)::pf
real(8), intent(in)::pb
real(8), intent(in)::pinv
integer, intent(in)::nmus
real(8), intent(in), dimension(nmus)::mus
integer, intent(in)::size


! output vbls and pars
!real(8), dimension(pp(3),sample), intent(inout)::tout
real(8), dimension(pp(3)), intent(out)::tout
!real(8), dimension(sample), intent(inout)::kout

! working variables
integer, dimension(12):: seedn
integer, dimension(nn,nn)::gsim, gsim2
integer i,j,k,rr,t, kk  ! indicators for loops
integer j0, k0, s0 ! indicators for loops
integer trand, tinsamp, diagcheck ! variables for large steps
integer nsq, posrow, poscol
integer, dimension(size)::samplinks
integer, dimension(size,2)::updlinks
real(8), dimension(pp(3))::tsim, tsim2 ! current simulated network stats
integer ktemp !indicator of current temperature
real(8) invp(mtemp) ! probability of inversion for each temperature
real(8) temper(mtemp) ! vector with temperature levels
real(8) rhotemp(mtemp) ! vector of pseudopriors
real(8) rho0,nn0 ! constants for stochastic approximation of pseudopriors (see Geyer-Thompson 1995 for details)
integer totsim  ! maximum nummber of simulations
real(8) ux,uux ! random number
real(8) odds ! odds of two networks
real(8) dtfu, dtfm, dtfv ! network statistics functions
real(8) dttemp(pp(3))
integer propmu ! mode proposed for update
integer curmu  ! current mode

curmu = 1

nsq = nn*2

! set prob of inversion
invp(1) = 0.0 ! .01 !.01 
do k=2,mtemp
invp(k)=invp(k-1)+.05
enddo

! set temperatures
temper(1)=1
do k=2,mtemp
	temper(k)=1 + .1*(k-1)
enddo

! set pseudopriors
rhotemp(:)=100000.0
rho0 = 1.0
nn0 = 10000000.0


! set initial simulated network the one provided inline
gsim = g

! compute observed network stats
tsim(:) = 0.0
call sstat(pp,nn,qq,gsim,x,dtt,tsim)

!tsim = tobs


!** random numbers seed
!CALL RANDOM_SEED(size =  k)
!seedn = seed1
!CALL RANDOM_SEED(put = seedn(1:k) )



!**************************************************!
!******* MAIN SIMULATION LOOP STARTS HERE *********!
!**************************************************!

tout(:)=0.0  ! initialize 
rr=1   ! sample loop indicator
ktemp =1  ! initialize temperature
totsim = skip*sample ! 10000000 ! max number of simulations

!***************** main loop
do t=1, totsim
    
    call random_number(uux)
    !print*, uux
    
    
    !******************* row step
    
    if (uux<pr) then
        !print*, 'row'
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose row/player i
        gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
		gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim(i,:) = 1 - gsim(i,:) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii

        endif ! (log(ux)<=min(0.0,odds))
    
    
    !***************************column step   
    
    elseif (pr <= uux .and. uux< pr+pc) then
        !print*, 'column' 
        
        call random_number(ux)
        i = ceiling(nn*ux) ! choose column i
        gsim(:,i) = 1 - gsim(:,i) !invert all links of column i
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
			   !print*, 'accepted'
		else
			gsim(:,i) = 1 - gsim(:,i) !invert all links of player i
			gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
		
        endif ! (log(ux)<=min(0.0,odds))
     
        
        
    !******************************* large step    
    
    elseif (pr+pc <= uux .and. uux < pr+pc+pf) then
            ! propose a uniform random sample of links
            ! to update (size = number of links)
            k0=1
            do while (k0<=size)
              call random_number(ux)
              trand = ceiling(nsq*ux)
              !trand = 1 + int(nsq*ux)
              tinsamp = 0
              do s0=1,k0-1
                  if (samplinks(s0)==trand) then
                      tinsamp = 1
                      
                  endif
              enddo !s0
              
              if (tinsamp==0) then
                  posrow = ceiling(1.0*trand/nn)
                  poscol = trand - (posrow-1)*nn
                  if (posrow .ne. poscol) then
                      samplinks(k0) = trand
                      k0=k0+1
                      !print*, 'accepted, k0 = ', k0-1
                  endif
              endif !(tinsamp==0)  
            enddo !(k0<=size)
            
            do s0=1,size
                updlinks(s0,1) = ceiling(1.0*samplinks(s0)/nn)
                updlinks(s0,2) = samplinks(s0) - (updlinks(s0,1)-1)*nn
            enddo !k0
            
            do s0=1,size
                gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
            enddo
            call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
            odds = 0.0
            do kk=1,pp(3)
                odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
            enddo ! kk
            call random_number(ux)
            if (log(ux)<=min(0.0,odds)) then
                tsim = tsim2
                !print*, 'accepted'
		    else
                do s0=1,size
                    gsim(updlinks(s0,1),updlinks(s0,2)) = 1 -&
                    &gsim(updlinks(s0,1),updlinks(s0,2)) !invert link
                enddo
            endif ! (log(ux)<=min(0.0,odds))
    

    !******************* bhamidi-like step
    
    elseif (pr+pc+pf <= uux .and. uux < pr+pc+pf+pb) then
        ! select which mode to propose
        !propmu = 1  ! if only 1 mode
        !if (nmus>1) then
        !    call random_number(ux)
        !    if (ux<0.5) then
        !        propmu = 1
        !    else
        !        propmu = 2
        !    endif !(ux<0.5)
        !endif !(nmus>1) 
        if (curmu == 1) then
            propmu = 2
        else
            propmu = 1
        endif
        
        ! new network is directed Erdos-Renyi 
        ! with density mus(propmu)
        gsim2 = 0
        do i=1,nn
            do j=1,nn
                if (i /= j) then
                    call random_number(ux)
                    if (ux < mus(propmu) ) then
                        gsim2(i,j) = 1
                    endif !(ux < mus(propmu) )
                endif 
            enddo !j 
        enddo  !i 

        ! compute network stats for proposed network
        call sstat(pp,nn,qq,gsim2,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            gsim = gsim2
            tsim = tsim2
            curmu = propmu
            !print*, 'accepted'
        endif ! (log(ux)<=min(0.0,odds))
    

    !******************* inversion step
    
    elseif (pr+pc+pf+pb <= uux .and. uux < pr+pc+pf+pb+pinv) then
        !print*, 'row'
        
        gsim = 1 - gsim !invert all links
		do i=1,nn
            gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
        enddo 
        
        call sstat(pp,nn,qq,gsim,x,dtt,tsim2)
        odds = 0.0
        do kk=1,pp(3)
            odds = odds + aa(kk)*(tsim2(kk)-tsim(kk))
        enddo ! kk
        call random_number(ux)
        if (log(ux)<=min(0.0,odds)) then
            tsim = tsim2
            !print*, 'accepted'
		else
			gsim = 1 - gsim !invert all links
    		do i=1,nn
                gsim(i,i) = 0 ! make sure you did not invert the link for g_ii
            enddo 
        endif ! (log(ux)<=min(0.0,odds))
    
        
    !********************************* one link step
    
    !elseif (uux <= pr+pc) then
     else   
        !print*, 'link'
        
        
        ! one-link per iteration 
        call random_number(ux)
        i = ceiling(nn*ux) ! choose player i
        call random_number(ux)
        j = ceiling(nn*ux) ! choose agent j
        
        do while (i==j)
            call random_number(ux)
            j = ceiling(nn*ux) ! choose agent j
        enddo ! while loop
        
        call random_number(ux)
        odds = 0.0
        dttemp(:) = 0.0 
        
        ! direct utility update
        do kk=1,pp(1)
            dttemp(kk) = dtfu(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        ! mutual utility update        
        if (gsim(j,i)==1) then  ! g_ji=1 then compute mutual utility
            do kk=pp(1)+1,pp(2)
                dttemp(kk) = dtfm(nn,qq,x,dtt(kk,1),i,j,dtt(kk,2))
            enddo
        endif
        ! indirect utility and popularity update
        do kk=pp(2)+1,pp(3)
            dttemp(kk) =  dtfv(nn,qq,gsim,x,dtt(kk,1),i,j,dtt(kk,2)) 
        enddo
        
       
        
        ! compute odds ration from the model
        do kk=1,pp(3)
            odds = odds + aa(kk)*dttemp(kk)
        enddo
        !odds = odds/temper(ktemp)
        
    
    
        ! update link g_ij
        if (gsim(i,j) ==1 ) then
            if (odds<0) then
                gsim(i,j)=0 
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) - dttemp(kk)
                enddo
                
            else
                if (log(ux)<-odds) then
                   gsim(i,j)=0 
                   do kk=1,pp(3)
                        tsim(kk) = tsim(kk) - dttemp(kk)
                    enddo
                endif ! (log(ux)<-odds)
                
            endif ! (odds<0)
            
        elseif (gsim(i,j)==0) then
            if (odds>=0) then
                gsim(i,j)=1
                do kk=1,pp(3)
                    tsim(kk) = tsim(kk) + dttemp(kk)
                enddo
            else
                if (log(ux)<odds) then
                    gsim(i,j)=1
                    do kk=1,pp(3)
                        tsim(kk) = tsim(kk) + dttemp(kk)
                    enddo
                    
                endif !(log(ux)<=odds)
                
            endif !(odds>=0)
        endif !(gsim(i,j) ==1 )
          
    endif !(ux<invp(ktemp))

    
    ! sample every skip iterations
    if (rr<=sample) then
        !if (t==skip*rr) then
            !tout(:,rr)=tsim
            tout=tsim
            !kout(rr)=ktemp
            !if (mod(rr,1000)==0) then
            !    print*, 'sample = ', rr
            !    print*, 'tsim = ', tsim
            !endif
    else        
        exit
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! introduce simulated tempering
    
enddo ! t
return
    end subroutine sample_ast_bham_est



  
    
    
!************************************
!* subroutine: bayesest_ex
!* version: 1.0
!* date: 03/05/2010
!* last modified: 01/10/2015
!* author: Angelo Mele
!* description: samples parameters from the model
!*              using an Exchange algorithm 
!*              Accelerated Simulated Tempering
!*              algorithm. 
!************************************

!subroutine bayesest_ex(nn, qq, pp, g, x, aa0, dtt, skip, parsamp, samp,&
!    &step, mtempnet, qup, seed1, pr, pc, pf, pb, pinv, nmus, mus, size,&
!    &mu0, sigma0, parout, accout)
subroutine bayesest_ex(nn, qq, pp, g, x, aa0, dtt, skip, parsamp, samp,&
    &step, mtempnet, qup, seed1, pr, pc, pf, pb, pinv, nmus, mus, size,&
    &mu0, sigma0, parout)


! input variables
! nn     = number of players
! qq     = number of variables
! pp     = number of parameters for u,m,v
! g      = network matrix to start simulations from
! x      = matrix of exogenous characteristics
! aa0    = vector of parameters to start simulations
! dtt    = vector of statistics for simulations
! skip   = number of network samples to skip
! parsamp= number of parameter simulations
! samp   = number of networks to sample
! step   = choleski matrix for RW for parameter sampling
! mtempnet= number of temperatures for simulated tempering
! qup    = probability of increasing temperature in sim tempering
! seed1  = seed for random number generation
! pr     = prob of changing a row of g
! pc     = prob of changing a column of g
! pf     = prob of large step
! pb     = prob of using Bhamidi-like step
! pinv   = prob of matrix inversion
! nmus   = number of modes of the gibbs distribution
! mus    = modes of the gibbs distribution
! size   = number of links updated in each iteration for large step
! mu0    = mean of the prior
! sigma0 = std dev of the prior

! output variables
! parout = matrix of sampled parameters
! accout = acceptance ratio
!***********************************************************


! input vbls and pars
integer, intent(in)::nn 
integer, intent(in)::qq
integer, intent(in)::pp(3)
integer, intent(in), dimension(nn,nn)::g
real(8), intent(in), dimension(nn,qq)::x
real(8), dimension(pp(3)), intent(in)::aa0
integer, dimension(pp(3),2), intent(in)::dtt
integer, intent(in)::skip
integer, intent(in)::parsamp
integer, intent(in)::samp
real(8), dimension(pp(3),pp(3)), intent(in):: step   
real(8), intent(in)::qup  
integer, intent(in)::mtempnet 
integer, intent(in)::seed1
real(8), intent(in)::pr
real(8), intent(in)::pc
real(8), intent(in)::pf
real(8), intent(in)::pb
real(8), intent(in)::pinv
integer, intent(in)::nmus
real(8), intent(in), dimension(nmus)::mus
integer, intent(in)::size
real(8), intent(in), dimension(pp(3))::mu0
real(8), intent(in), dimension(pp(3))::sigma0


! output variables
real(8), intent(inout), dimension(parsamp,pp(3))::parout
!real(8), intent(out)::accout
real(8) accout

! working vbls
real(8), dimension(pp(3)):: tsamp  ! last networks from the network simulation
!real(8), dimension(pp(3),samp)::tout !sample of networks from the network simulation
!real(8), dimension(samp)::kout ! sample temperatures from network simulation
!real(8) kout ! sample temperatures from network simulation
!real(8), dimension(pp(3),pp(3)):: stderr ! var cov from networks samples (not needed)
integer t,i,r,k ! counters for loops
real(8) ux
real(8) pmh
real(8), dimension(pp(3)):: uux ! vector of random variables for parameter update
real(8), dimension(pp(3)):: thetap ! proposed parameter vector 
real(8), dimension(pp(3)):: theta  ! current parameter vector
integer nacc, nrej  ! number of acceptances and rejections for updates
real(8) racc ! ratio of acceptances
!real(8), dimension(mtemp):: rhotemp
!real(8), dimension(mtemp):: temper
!real(8) qdown
!real(8) rst
real(8) normal
!real(8), dimension(samp)::odds
!real(8) maxodds
!real(8) consratio ! ratio of constants
!real(8) mu0,sigma0, mu1,sigma1  ! parameters for update
real(8), dimension(pp(3))::tobs
integer, dimension(12):: seedn
real(8) mu1, sigma1



! parameters for random walk
mu1 = 0.0
sigma1 = 1.0

! compute observed network stats
tobs(:) = 0.0
tsamp(:) = 0.0
call sstat(pp,nn,qq,g,x,dtt,tobs)


!** random numbers seed
CALL RANDOM_SEED(size =  k)
seedn = seed1
CALL RANDOM_SEED(put = seedn(1:k) )


! initialize parameters
theta = aa0


!**** start main loop
t = 1
racc = 0.0
nacc = 0
nrej = 0
thetap(:) = 0.0
do while (t<=parsamp)
	
	
	!**** propose a new theta
	do i=1,pp(3)
		uux(i)=normal(mu1,sigma1)  ! random number for random walk update
	enddo
	do i=1,pp(3)
		thetap(i) = theta(i) + dot_product(step(i,:),uux) ! propose new element i of parameter vector
		!print*, 'dot_product(step(i,:),uux)', dot_product(step(i,:),uux)
    enddo !i  

    
	!**** collect sample to compute constant ratio
!    k = 1
    call sample_ast_bham_est(nn, qq, pp, g, x, thetap, dtt, skip, samp,&
         &mtempnet, qup, tsamp, seed1, pr, pc, pf, pb, pinv, nmus, mus, size)

    !call sample_ast_rc_est(nn,qq,pp,thetap,dt,skip,samp,mtempnet,qup,tsamp,kout,seed1,pr,pc)
    !tsamp = tout(:,samp)

    
    !**** accept new parameter using AST approximate ratio
	!pmh = dot_product(thetap-theta,(tsamp-tobs)/temper(ktemp))  - .5*sum(thetap**2-theta**2)/sigma0
	!pmh = dot_product(theta-thetap,(tsamp-tobs)/temper(ktemp))  - .5*sum(thetap**2-theta**2)/sigma0
	!pmh = dot_product(theta-thetap,(tsamp-tobs))  -.5*sum(  (thetap**2-theta**2)/sigma0 )
	pmh = dot_product(theta-thetap,tsamp-tobs)  
    do i=1,pp(3)
        pmh = pmh - .5*(thetap(i)**2-theta(i)**2)/sigma0(i)
    enddo 
    
	if (pmh>0.0) then
		theta = thetap
	    parout(t,:) = thetap
        !nacc=nacc+1
		!racc = 1.0*nacc
		!racc = racc/t
        !accout = racc
		
	!**** if fp smaller than current fcn, then use MH to update
	else 
		call random_number(ux)
		if (log(ux)<=pmh) then
			theta = thetap
			parout(t,:) = thetap
			!nacc  = nacc + 1
			!racc = 1.0*nacc
			!racc = racc/t
			!accout = racc
		else 
			!nrej  = nrej + 1
            theta = theta
			parout(t,:) = theta
			!racc = 1.0*nacc
			!racc = racc/t
			!accout = racc
		endif
	
	endif

	t = t + 1
	!endif
	



	! update temperature
	


	!! update temperature if AST
	!if(mtemp>1) then
 !
	!if (ktemp==1) then
	!	rst  = min(0.0,log(rhotemp(2)/rhotemp(1))+dot_product(theta,tobs)*(1/temper(2)-1/temper(1))+log(qup/1))
 !
	!	call random_number(ux)
	!	if (log(ux)<=rst) then
	!		ktemp=2
	!	endif !(log(ux)<=rst)
	!elseif (ktemp==mtemp) then
	!	rst  = min(0.0,log(rhotemp(mtemp-1)/rhotemp(mtemp))+dot_product(theta,tobs)*(1/temper(mtemp-1)-1/temper(mtemp))+log(1/qdown))
	!	call random_number(ux)
	!	if (log(ux)<=rst) then
	!		ktemp=mtemp-1
	!	endif !(log(ux)<=rst)
	!else 
	!	call random_number(ux)
	!	if (ux<=qup) then
	!		rst  = min(0.0,log(rhotemp(ktemp+1)/rhotemp(ktemp))+dot_product(theta,tobs)*(1/temper(ktemp+1)-1/temper(ktemp))+log(qdown/qup))
	!		call random_number(ux)
	!		if (log(ux)<=rst) then
	!			ktemp=ktemp+1
	!		endif !(log(ux)<=rst)
	!	else
	!		rst  = min(0.0,log(rhotemp(ktemp-1)/rhotemp(ktemp))+dot_product(theta,tobs)*(1/temper(ktemp-1)-1/temper(ktemp))+log(qup/qdown))
	!		call random_number(ux)
	!		if (log(ux)<=rst) then
	!			ktemp=ktemp-1
	!		endif !(log(ux)<=rst)		
	!	endif !(ux=<qup)
	!endif	
 !
	!endif		
	
enddo ! t main loop
return
end subroutine bayesest_ex


function normal(mean,sigma) !returns a normal distribution 

implicit none 
real(8) normal,tmp 
real(8), intent(in)::mean,sigma 
integer flag 
real(8) fac,gsave,rsq,r1,r2 
save flag,gsave 
data flag /0/ 
	if (flag.eq.0) then 
		rsq=2.0
		do while(rsq.ge.1.0.or.rsq.eq.0.0) ! new from for do 
			call random_number(r1)
			call random_number(r2)
			r1=2.0*r1-1.0
			r2=2.0*r2-1.0 
			rsq=r1*r1+r2*r2 
		enddo 
		fac=sqrt(-2.0*log(rsq)/rsq) 
		gsave=r1*fac 
		tmp=r2*fac 
		flag=1 
	else 
		tmp=gsave 
		flag=0 
	endif 
normal=tmp*sigma+mean 
return 
end function normal
