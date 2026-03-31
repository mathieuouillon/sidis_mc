!------------------------------------------------------------------------------
! Initialization of variables to book in the ntuple
!------------------------------------------------------------------------------
subroutine InitKin2Book
    use event_info_module, only: Nu, Q22, XBj, W
    use particles_module, only: TrkGS, Nb_part
    implicit none


    Nu = 0.
    Q22 = 0
    XBj = 0
    W = 0
    TrkGS = 0
    Nb_part = 0

end subroutine InitKin2Book

!------------------------------------------------------------------------------
! Computation of variables to book in the ntuple
!------------------------------------------------------------------------------
subroutine ComputV
    use event_info_module, only: Nu, Q22, XBj, W, y_ele
    use particles_module, only: TrkGS, Nb_part, acc_part, ch_part, id_part, id_mother, &
        px_part, py_part, pz_part, p_part, E_part, m_part, z_part, th_part, tt_part, &
        Pts_part, phih_part, phi_part, Xf_part, Als
    use pythia_commons, only: P, N, K
    use config_module, only: iColl, EColl, iAccept, iAlert
    use acceptance_module, only: pro_acc, deu_acc, tri_acc, he3_acc, he4_acc
    use misc_module, only: ranf
    implicit none


    integer ip
    ! for calculation of Phih
    real A1, A2, A3, AA
    real B1, B2, B3, BB
    real phi_ele
    ! function
    integer clas12_accept
    ! Initial kinematics
    real eip
    real nip, nie
    ! C func
    integer alertaccept

    eip = P(1, 4)
    nip = 0
    nie = P(2, 5)

    if (iColl .eq. 1) then
        nie = EColl
        nip = sqrt(EColl**2 - P(2, 5)**2)
    end if

    ! Booking the hbook
    do ip = 1, N
        if (k(ip, 2) .eq. 22 .and. k(ip, 3) .eq. 1) then
            Nu = (p(ip, 4)*nie + p(ip, 3)*nip)/P(2, 5)
            Q22 = P(ip, 5)**2
            if (Nu .ne. 0) XBj = Q22/2/P(2, 5)/Nu
            if (Nu .ne. 0) y_ele = Nu*P(2, 5)/eip/(nie + nip)
            W = 0.
            W = dsqrt((nie + p(ip, 4))**2 - (-nip + p(ip, 3))**2 - p(ip, 2)**2 - p(ip, 1)**2)
            TrkGS = ip
            phi_ele = atan2(p(ip, 2), p(ip, 1))*57.2958 + 210

        end if

        if (abs(k(ip, 2)) .eq. 211 .or. (k(ip, 2) .eq. 11 .and. k(ip, 1) .eq. 1) &
            .or. abs(k(ip, 2)) .eq. 321 .or. k(ip, 2) .eq. 310 &
            .or. (abs(k(ip, 2)) .eq. 2212 .and. k(ip, 1) .eq. 1) &
            .or. (abs(k(ip, 2)) .eq. 2112 .and. k(ip, 1) .eq. 1) &
            .or. (abs(k(ip, 2)) .gt. 10000 .and. k(ip, 1) .eq. 1) &
            .or. (k(ip, 2) .eq. 22) .or. (k(ip, 2) .eq. 111) &
            ) then

            Nb_part = Nb_part + 1
            acc_part(Nb_part) = 0
            ch_part(Nb_part) = 0
            if (k(ip, 2) .eq. 11 .or. &
                k(ip, 2) .eq. -211 .or. &
                k(ip, 2) .eq. -321 .or. &
                k(ip, 2) .eq. -2212) ch_part(Nb_part) = -1
            if (k(ip, 2) .eq. 211 .or. &
                k(ip, 2) .eq. 321 .or. &
                k(ip, 2) .eq. 1000010030 .or. &
                k(ip, 2) .eq. 1000010020 .or. &
                k(ip, 2) .eq. 2212) ch_part(Nb_part) = 1
            if (k(ip, 2) .eq. 1000020030) ch_part(Nb_part) = 2

            id_part(Nb_part) = k(ip, 2)
            id_mother(Nb_part) = k(k(ip, 3), 2)
            px_part(Nb_part) = p(ip, 1)
            py_part(Nb_part) = p(ip, 2)
            pz_part(Nb_part) = p(ip, 3)
            p_part(Nb_part) = sqrt(p(ip, 4)**2 - p(ip, 5)**2)
            E_part(Nb_part) = p(ip, 4)
            m_part(Nb_part) = p(ip, 5)
            z_part(Nb_part) = (p(ip, 4)*nie + p(ip, 3)*nip)/Nu/P(2, 5)
            th_part(Nb_part) = 57.2957795*acos(pz_part(Nb_part)/sqrt(p(ip, 1)**2 + p(ip, 2)**2 + pz_part(Nb_part)**2))
            phi_part(Nb_part) = atan2(p(ip, 2), p(ip, 1))*57.2958 + 30
            if (phi_part(Nb_part) .lt. 0) phi_part(Nb_part) = phi_part(Nb_part) + 360

            ! TODO update and make frame independent for collider kinematics
            A1 = sin(phi_ele)
            A2 = -cos(phi_ele)
            A3 = 0
            AA = A1**2 + A2**2 + A3**2

            B1 = p(TrkGS, 2)*p(ip, 3) - p(TrkGS, 3)*p(ip, 2)
            B2 = p(TrkGS, 3)*p(ip, 1) - p(TrkGS, 1)*p(ip, 3)
            B3 = p(TrkGS, 1)*p(ip, 2) - p(TrkGS, 2)*p(ip, 1)
            BB = B1**2 + B2**2 + B3**2

            phih_part(Nb_part) = acos((A1*B1 + A2*B2 + A3*B3)/sqrt(AA*BB))*57.2958

     tt_part(Nb_part) = (Nu - p(ip, 4))**2 - (p(TrkGS, 1) - p(ip, 1))**2 - (p(TrkGS, 2) - p(ip, 2))**2 - (p(TrkGS, 3) - p(ip, 3))**2

            Pts_part(Nb_part) = (z_part(Nb_part)*Nu)**2 - m_part(Nb_part)**2 &
                                - ((-p(TrkGS, 4)*p(ip, 4) + p(TrkGS, 1)*p(ip, 1) &
                                    + p(TrkGS, 2)*p(ip, 2) + p(TrkGS, 3)*p(ip, 3) &
                                    + z_part(Nb_part)*Nu**2)**2/(Nu**2 + Q22))

            Xf_part(Nb_part) = &
                (z_part(Nb_part)*P(2, 5)*Nu**2 &
                 - z_part(Nb_part)*Q22*Nu &
                 - (P(2, 5) + Nu)* &
                 (p(TrkGS, 4)*p(ip, 4) - p(TrkGS, 1)*p(ip, 1) &
                  - p(TrkGS, 2)*p(ip, 2) - p(TrkGS, 3)*p(ip, 3))) &
                /sqrt(W**2/4 - p(ip, 5)**2)/W &
                /sqrt(Nu**2 + Q22)
            Als(Nb_part) = (P(ip, 4) - (z_part(Nb_part)*Nu**2 &
                                        - (p(TrkGS, 4)*p(ip, 4) - p(TrkGS, 1)*p(ip, 1) &
                                           - p(TrkGS, 2)*p(ip, 2) - p(TrkGS, 3)*p(ip, 3))) &
                            /sqrt(W**2/4 + Q22))/P(ip, 5)
            acc_part(Nb_part) = 1

            if (iAccept .eq. 1) then
                acc_part(Nb_part) = 0
                if (abs(k(ip, 2)) .eq. 11 .or. abs(k(ip, 2)) .eq. 211 .or. &
                    abs(k(ip, 2)) .eq. 321 .or. k(ip, 2) .eq. 22 .or. &
                    k(ip, 2) .eq. 2112) then
                    acc_part(Nb_part) = clas12_accept(k(ip, 2), p(ip, 1), p(ip, 2), p(ip, 3))
                else if (k(ip, 2) .gt. 2200 .and. p_part(Nb_part) .lt. 0.5 .and. th_part(Nb_part) .lt. 180) then
                    ! if(iAlert.eq.1) then
                    !       acc_part(Nb_part) =             alertaccept(k(ip,2),p_part(Nb_part),th_part(Nb_part))
                    if (iAlert .eq. 0) then
                    if (k(ip, 2) .eq. 2212) then
                        if (ranf(0) .lt. pro_acc(int(p_part(Nb_part)/.02) + 1, &
                                                 int(th_part(Nb_part)/7.2) + 1)) &
                            acc_part(Nb_part) = 1
                    else if (k(ip, 2) .eq. 10002) then
                        acc_part(Nb_part) = 0
                    else if (k(ip, 2) .eq. 10102) then
                        if (ranf(0) .lt. deu_acc(int(p_part(Nb_part)/.02) + 1, &
                                                 int(th_part(Nb_part)/7.2) + 1)) &
                            acc_part(Nb_part) = 1
                    else if (k(ip, 2) .eq. 10103) then
                        if (ranf(0) .lt. tri_acc(int(p_part(Nb_part)/.02) + 1, &
                                                 int(th_part(Nb_part)/7.2) + 1)) &
                            acc_part(Nb_part) = 1
                    else if (k(ip, 2) .eq. 10202) then
                        acc_part(Nb_part) = 0
                    else if (k(ip, 2) .eq. 10203) then
                        if (ranf(0) .lt. he3_acc(int(p_part(Nb_part)/.02) + 1, &
                                                 int(th_part(Nb_part)/7.2) + 1)) &
                            acc_part(Nb_part) = 1
                    else if (k(ip, 2) .eq. 10204) then
                        if (ranf(0) .lt. he4_acc(int(p_part(Nb_part)/.02) + 1, &
                                                 int(th_part(Nb_part)/7.2) + 1)) &
                            acc_part(Nb_part) = 1
                    end if
                    end if
                end if
            end if

        end if
    end do

end subroutine ComputV

subroutine init_recoil
    use acceptance_module, only: pro_acc, deu_acc, tri_acc, he3_acc, he4_acc
    implicit none

    integer i, j
    real a, b, c, d, e


    write (*, *) 'Init recoil acceptance routine'

    open (unit=21, file="datafiles/recoil_acc/acc_proton.dat")
    open (unit=22, file="datafiles/recoil_acc/acc_deuteron.dat")
    open (unit=23, file="datafiles/recoil_acc/acc_trition.dat")
    open (unit=24, file="datafiles/recoil_acc/acc_helion.dat")
    open (unit=25, file="datafiles/recoil_acc/acc_alpha.dat")

    do i = 1, 25
        do j = 1, 25
            read (21, *) a, b, c, d, e
            pro_acc(i, j) = e
            read (22, *) a, b, c, d, e
            deu_acc(i, j) = e
            read (23, *) a, b, c, d, e
            tri_acc(i, j) = e
            read (24, *) a, b, c, d, e
            he3_acc(i, j) = e
            read (25, *) a, b, c, d, e
            he4_acc(i, j) = e
        end do
    end do

    close (21)
    close (22)
    close (23)
    close (24)
    close (25)

end subroutine init_recoil
