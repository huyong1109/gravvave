MODULE sol
USE param

CONTAINS

!=======================
SUBROUTINE presol

end SUBROUTINE presol

subroutine pcg(X,B)

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X ! on input, an initial guess for the solution
                         ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------


   real (r8) :: &
      eta0,eta1,rr ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R, &! residual (b-Ax)
      S, &! conjugate direction vector
      Q,WORK0,WORK1 ! various cg intermediate results

   character (char_len) :: &
      noconvrg ! error message for no convergence
    ! hyfile ! editbyhy



!-----------------------------------------------------------------------
!
! compute initial residual and initialize S
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)


      call btrop_operator(S,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - S(:,:,iblock)
      S(:,:,iblock) = c0

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! initialize fields and scalars
!
!-----------------------------------------------------------------------

   eta0 =c1
   solv_sum_iters = solv_max_iters

!-----------------------------------------------------------------------
!
! iterate

!-----------------------------------------------------------------------

   iter_loop: do m = 1, solv_max_iters

!-----------------------------------------------------------------------
!
! calculate (PC)r
! diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         if (lprecond) then
            call preconditioner(WORK1,R,this_block,iblock)
         else
            where (A0(:,:,iblock) /= c0)
               WORK1(:,:,iblock) = R(:,:,iblock)/A0(:,:,iblock)
            elsewhere
               WORK1(:,:,iblock) = c0
            endwhere
         endif

         WORK0(:,:,iblock) = R(:,:,iblock)*WORK1(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! update conjugate direction vector s
!
!-----------------------------------------------------------------------

      !*** (r,(PC)r)
      eta1 = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         S(:,:,iblock) = WORK1(:,:,iblock) + S(:,:,iblock)*(eta1/eta0)

!-----------------------------------------------------------------------
!
! compute As
!
!-----------------------------------------------------------------------

         call btrop_operator(Q,S,this_block,iblock)
         WORK0(:,:,iblock) = Q(:,:,iblock)*S(:,:,iblock)

      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! compute next solution and residual
!
!-----------------------------------------------------------------------


      eta0 = eta1
      eta1 = eta0/global_sum(WORK0, distrb_tropic, &
                             field_loc_center, RCALCT_B)

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         X(:,:,iblock) = X(:,:,iblock) + eta1*S(:,:,iblock)
         R(:,:,iblock) = R(:,:,iblock) - eta1*Q(:,:,iblock)

         if (mod(m,solv_ncheck) == 0) then

            call btrop_operator(R,X,this_block,iblock)
            R(:,:,iblock) = B(:,:,iblock) - R(:,:,iblock)
            WORK0(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! test for convergence
!
!-----------------------------------------------------------------------

      if (mod(m,solv_ncheck) == 0) then

         call update_ghost_cells(R, bndy_tropic, field_loc_center,&
                                                 field_type_scalar)

         rr = global_sum(WORK0, distrb_tropic, &
                         field_loc_center, RCALCT_B) ! (r,r)

         if (rr < solv_convrg) then
            solv_sum_iters = m
            exit iter_loop
         endif

      endif

   enddo iter_loop
      call update_ghost_cells(X, bndy_tropic, field_loc_center, &
                                              field_type_scalar)

   rms_residual = sqrt(rr*resid_norm)

   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
         write(noconvrg,'(a45,i11)') &
           'Barotropic solver not converged at time step ', nsteps_total
         call exit_POP(sigAbort,noconvrg)
      endif
   endif

!-----------------------------------------------------------------------

 end subroutine pcg
 
subroutine btrop_operator(AX,X)

! !DESCRIPTION:
! This routine applies the nine-point stencil operator for the
! barotropic solver. It takes advantage of some 9pt weights being
! shifted versions of others.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      AX ! nine point operator result (Ax)

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j ! dummy counters

!-----------------------------------------------------------------------

   AX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) =abs( A0 (i ,j ,bid)*X(i ,j ,bid)) + &
                   abs( AN (i ,j ,bid)*X(i ,j+1,bid)) + &
                   abs( AN (i ,j-1,bid)*X(i ,j-1,bid)) + &
                   abs( AE (i ,j ,bid)*X(i+1,j ,bid)) + &
                   abs( AE (i-1,j ,bid)*X(i-1,j ,bid)) + &
                   abs( ANE(i ,j ,bid)*X(i+1,j+1,bid)) + &
                   abs( ANE(i ,j-1,bid)*X(i+1,j-1,bid)) + &
                   abs( ANE(i-1,j ,bid)*X(i-1,j+1,bid)) + &
                   abs( ANE(i-1,j-1,bid)*X(i-1,j-1,bid))
   end do
   end do

!-----------------------------------------------------------------------

 end subroutine btrop_operator_abs
