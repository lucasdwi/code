ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Main program: computes the p-value of the unconditional test for testing
c one and two-sided hypotheses about the means of two Poisson
c distributions.
c
c INPUT:
c iside = side of the test; 1 for right-sided, 2 for two-sided
c alpha = nominal level of the test
c ki = count of the ith population, i = 1,2
c ni = sample size from the ith population, i=1,2
c d = the difference mean1 - mean2 under the H0
c
c OUTPUT:
c p-value = p-value of the unconditional test
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	implicit doubleprecision (a-h,o-z)
      
1     print*,'Enter 1 for right-tail test or' 
	print*,'      2 for two-tail test'
	print*,' '

	read*, iside
	if(iside .ne. 1 .and. iside .ne. 2) then
	  goto 1
	end if
      
2     print*,'Enter the sample counts, k1 and k2'
	read *, k1, k2
	if(k1 .lt. 0.0d0 .or. k2 .lt. 0.0d0) then
	  print*, 'k1 >= 0, k2 >= 0'
	  goto 2
	end if    
	print*,' '
      
3     print*,'Enter the sample sizes n1 and n2'
	read*, n1, n2
	if(n1 .le. 0.0 .or. n2 .le. 0.0) then 
        print*, 'n1 > 0, n2 > 0'
	  goto 3
	end if
	print*,' '
	
4     print*,'Enter the value of mean1 - mean2 under H0:'
	read*, d
	if(d .lt. 0.0d0) then 
        print*, 'd >= 0'
	  goto 4
	end if
	print*,' '

	elhatk = 1.0d0*(k1+k2)/(n1+n2)-d*n1/(n1+n2)
        print*, 'elhatk is', elhatk
	var = (1.0d0*k1/n1**2 + 1.0d0*k2/n2**2)
        print*, 'var is', var
	t_k1k2 = (1.0d0*k1/n1-1.0d0*k2/n2-d)/sqrt(var)
        print*, 't_k1k2 is', t_k1k2

	call poistest(iside, n1, n2, elhatk, t_k1k2, d, pvalue)
      
	print*, 'The p-value is', pvalue
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Program for computing the p-value of the unconditional test
c In the first subroutine, the sum over i1 is carried out
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine poistest(iside, n1, n2, elhatk, t_k1k2, d, pvalue)
	implicit doubleprecision(a-h,o-z)

	pvalue = 0.0d0

c computing estimates of el1*n1 and el2*n2 under H_0
	elhat1 = n1*(elhatk+d)
        print*, 'elhat1 is', elhat1
	elhat2 = n2*elhatk
		print*, 'elhat2 is', elhat2

c computing the modes 
	i1mode = int(elhat1)
	i2mode = int(elhat2)

c initializing the probability at the i1mode
	pi1mode = poipr(i1mode, elhat1)
	pi1 = pi1mode
        print*, 'pi1 is', pi1

c initializing the probability at the i2mode
        print*, 'initializing the probability at the i2 mode'
	pi2mode = poipr(i2mode, elhat2)
        print*, 'pi2mode is',pi2mode

	do i1 = i1mode, 1000
	  if(pi1 .lt. 1.0d-07)       goto 1
          print*, 'about to call sumi2'
	  call sumi2(iside, n1, n2, elhat2, t_k1k2, i1, pi1, i2mode, 
     &             pi2mode, d, pvalue)
	 print*, 'pvalue after call to sumi2 is ', pvalue
	  pi1 = elhat1*pi1/(i1+1.0d0)
	end do

1	i1 = i1mode-1
        print*, 'in label 1 i1 is now', i1
	pi1 = pi1mode
	pi1 = i1mode*pi1/elhat1
	
	do i1 = i1mode-1, 0, -1
          if(pi1 .lt. 1.0d-07) print*, 'pi1 is lt 1e-7!!'
	  if(pi1 .lt. 1.0d-07) return
	  call sumi2(iside, n1, n2, elhat2, t_k1k2, i1, pi1, i2mode, 
     &             pi2mode, d, pvalue)
	  pi1 = i1*pi1/elhat1
	end do
	print*, 'at the end of poistest pi1 is now', pi1
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Here, we carry out the sum over i2 to compute the p-value of the E-test 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine sumi2(iside, n1, n2, elhat2, t_k1k2, i1, pi1, i2mode,
     &                 pi2mode, d, pvalue)
	implicit doubleprecision(a-h,o-z)
	pi2 = pi2mode
	print*, 'pvalue just inside sumi2 is ', pvalue
	print*,'d just inside sum is', d

	do i2 = i2mode, 1000
	  if(pi2 .lt. 1.0d-07) goto 1
	  elhati1 = 1.0d0*i1/n1
	  elhati2 = 1.0d0*i2/n2
	  diffi = elhati1 - elhati2 - d 
	  print*, 'diffi is ', diffi
	  var = (1.0d0*elhati1/n1 + 1.0d0*elhati2/n2)
	  if(iside .eq. 1) then    
	    if(1.0d0*i1/n1 - 1.0d0*i2/n2 .le. d) then
	      t_i1i2 = 0.0d0
	    else
	      t_i1i2 = diffi/sqrt(var)
	    end if
	    if(t_i1i2 .ge. t_k1k2) pvalue = pvalue + pi1*pi2
	  else if(iside .eq. 2) then
	    print*, 'entering double tail pvalue is', pvalue  
		print*, 'entering double tail pi1 is', pi1
		print*, 'entering double tail pi2 is', pi2
	    if(dabs(1.0d0*i1/n1 - 1.0d0*i2/n2) .le. d) then
	      t_i1i2 = 0.0d0
	    else
	      t_i1i2 = diffi/sqrt(var)
	    end if
	    if(dabs(t_i1i2) .ge. dabs(t_k1k2)) pvalue = pvalue + pi1*pi2
		print*, 'pvalue inside iside2 loop is', pvalue
	  end if  
	  pi2 = elhat2*pi2/(i2+1.0d0)
	end do
c
1     i2 = i2mode-1 
	pi2 = pi2mode
	pi2 = i2mode*pi2/elhat2
        print*, 'pi1 at end of label 1 is', pi1
		print*, 'pi2 at end of label 1 is', pi2
		print*, 't_i1i2 at end of label 1 is', t_i1i2
		print*, 't_k1k2 at end of label 1 is', t_k1k2
		print*, 'pvalue at end of label 1 is', pvalue

	do i2 = i2mode-1, 0, -1
	  if(pi2 .lt. 1.0d-07) return
	  elhati1 = 1.0d0*i1/n1
	  elhati2 = 1.0d0*i2/n2
	  diffi = elhati1 - elhati2 - d
	  var = (1.0d0*elhati1/n1 + 1.0d0*elhati2/n2)
	  if(iside .eq. 1) then    
	    if(1.0d0*i1/n1 - 1.0d0*i2/n2 .le. d) then
	      t_i1i2 = 0.0d0
	    else
	      t_i1i2 = diffi/sqrt(var)
	    end if
	    if(t_i1i2 .ge. t_k1k2) pvalue = pvalue + pi1*pi2
	  else if(iside .eq. 2) then    
	    if(dabs(1.0d0*i1/n1 - 1.0d0*i2/n2) .le. d) then
	      t_i1i2 = 0.0d0
	    else
	      t_i1i2 = diffi/sqrt(var)
	    end if
	    if(dabs(t_i1i2) .ge. dabs(t_k1k2)) pvalue = pvalue + pi1*pi2
	  end if  
	  pi2 = i2*pi2/elhat2
	end do
	print*, 'pi1 at end of sumi2 is', pi1
	print*, 'pi2 at end of sumi2 is', pi2
	print*, 't_i1i2 at end of sumi2 is', t_i1i2
	print*, 't_k1k2 at end of sumi2 is', t_k1k2
	print*, 'pvalue at end of sumi2 is', pvalue
	print*, '============'
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This program computes the P(X = k), where X is a Poisson random
c variable with mean defective rate = el, # of defective items = k
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	double precision function poipr(k, el)
	implicit doubleprecision(a-h,o-z)
	  
	ek = k*1.0d0
	poipr = dexp(-el+ek*dlog(el)-alng(ek+1.0d0))
	end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c Logarithmic gamma function = alng(x),  x > 0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        double precision function alng(x)
        implicit doubleprecision (a-h, o-z)
        double precision b(8)
        logical indx
        data b/8.33333333333333d-2, 3.33333333333333d-2,
     +         2.52380952380952d-1, 5.25606469002695d-1,
     +         1.01152306812684d0,  1.51747364915329d0,
     +         2.26948897420496d0,  3.00991738325940d0/ 
        if(x .lt. 8.0d0) then
        xx = x + 8.0d0
        indx = .true.
        else
        indx = .false.
        xx = x
        end if
c
        fterm = (xx-0.5d0)*dlog(xx) - xx + 9.1893853320467d-1
        sum = b(1)/(xx+b(2)/(xx+b(3)/(xx+b(4)/(xx+b(5)/(xx+b(6)
     +/(xx+b(7)/(xx+b(8))))))))
        alng = sum + fterm
        if(indx)  alng = alng-dlog(x+7.0d0)-dlog(x+6.0d0)-dlog
     +        (x+5.0d0)-dlog(x+4.0d0)-dlog(x+3.0d0)-dlog(x+2.0d0)
     +        -dlog(x+1.0d0)-dlog(x)
        end
