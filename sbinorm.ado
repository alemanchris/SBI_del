program define sbinorm, eclass
// Version 1.0 (This version is not public, please refrain from sharing)
// This version is still preliminary and incomplete, for the complete version refer
// to MATLAB version 3.1 SBI_norm.m
// sbinorm, runs the Stage Based Identification routine descrived in 
// Aleman,  Busch,  Ludwig,  and  Santaeul`alia-Llopis  (2020)  
// Please address any comments to christian.c.aleman[at]gmail.com

        version 14
 
        syntax varlist, tp(integer)
		marksample touse
		fvexpand `varlist' 		
        local cnames `r(varlist)'
        mata tua = `tp'
		mata SBI_parse("`cnames'", "`touse'", tua)

end
//------------------------------------------------------------------------------
mata
// define interpolating function 
//-----------------------------------------------------------------------------
function interp1(x,y,xq)
// Careful the is not interpolating zero, fix that!
{
nel1 = rows(xq)
nel2 = rows(x)
yq = J(nel1,1,-99999)

	for (i=1; i<=nel1; i++) {     // for loop 
    diff_aux = x:-xq[i,1]
	ind_sign = sign(diff_aux)
    aux_del = select(ind_sign, ind_sign[.,1]:<0)
    pos = colsum(abs(aux_del))
    ind = 0
    

		  if (pos > 0) { 
				if (pos<nel2){
				ind = 1
				}
				else{			
				}
		   }
		   else{ 
		   }
			
		  if (ind==1){
				y0 = y[pos,1]
				y1 = y[pos+1,1]
				x0 = x[pos,1]
				x1 = x[pos+1,1]
				aux_del2 = y0 + ((y1 - y0)/(x1 - x0)) * (xq[i,1] - x0)
				yq[i,1] =  aux_del2 //&J(1,1,aux_del2)	  
		   }
		   else{
		   //yq[i,1]=&J(1,1,NULL)
		   }
		}                          // end for loop
   
return(yq)
}
//-----------------------------------------------------------------------------
function func_stage(time,pphi)
{
    tl = rows(time)  // lenght
    y = J(tl,1,pphi[2,1])
	nm = rows(pphi)
    for (i=3; i<=(nm); i++){
        y = y + (time:^(i-2)):*pphi[i,1]
    }
return(y)
}
//-----------------------------------------------------------------------------
function func_dist(todo,pphi,yT,yC,tu,we,fv,g,H){
	pphi2 = pphi'
	dtaT = yT
	dtaC = yC
	svec =  range(1,tu,1)
	tvec = func_stage(svec,pphi2)
	dtaTT = dtaT[1..tu,1]
	dtaCC = dtaC[1..tu,1]
	
	aux = interp1(tvec,dtaCC,svec)
	// Eliminate -99999
	 pred_C_sc = select(aux, aux[.,1]:>-99999)
	 pred_C_sc = pphi2[1,1]:*pred_C_sc 
	 svec = select(svec, aux[.,1]:>-99999)
	 dtaTT = select(dtaTT, aux[.,1]:>-99999)
	 
	 if (rows(dtaTT)<5) {
		dist = 1e10 
	 }
	 else{
		dist = dtaTT:-pred_C_sc;
	 }

    tl = rows(dist)  
	if (we==1){
		wght = sqrt(range(1,tl,1))		
	}
	 else{
		wght = J(tl,1,1)
	}
	
    wdist = wght:*dist
  
    fv = 0.5:*((wdist')*wdist)

   //return(fv)
}

//-----------------------------------------------------------------------------
function kumsum(x)
{
// take column verctors
	lgth = rows(x)
    y = J(lgth,1,0)
	y[1,1] = x[1,1]
    for (i=2; i<=(lgth); i++){
        y[i,1] = y[i-1,1] + x[i,1]
    }
return(y)
}
//-----------------------------------------------------------------------------
//**************************************************************
// CORE ROUTINE STARTS HERE
//**************************************************************
void SBI_parse(string scalar indepvars,  string scalar touse, tua){
		// Fetch data
		real matrix data_mat
        data_mat  = st_data(., indepvars, touse)
		
		// Rename
		yC = data_mat[.,1]
		yT = data_mat[.,2]
		time = data_mat[.,3]
		nm = 2							// choose number of param 2: linear 3 quadratic
		we = 0							// Weighted minimization
		tu =  tua
	
		// Initial guess 
		x0 = (0.91,5,1.05)
		
		 // Perform minimization
		 s=optimize_init()
		 optimize_init_which(s, "min")
		 optimize_init_evaluator(s,&func_dist())
		 optimize_init_argument(s, 1, yT)
		 optimize_init_argument(s, 2, yC)
		 optimize_init_argument(s, 3, tu)
		 optimize_init_argument(s, 4, we)
		 optimize_init_conv_ptol(s, 1e-16)
		 optimize_init_conv_vtol(s, 1e-16)
		 optimize_init_conv_maxiter(s, 150)
		 optimize_init_params(s,x0)
		 cppsi=optimize(s)
		 // dist2 = func_dist((1,2,1.01),yT,yC,tu,we)
		 
		 // Compute gamma (not trapeziodal), add half the day		 
		tu_norm = func_stage(tu,cppsi')	   		 
		svec = time
		tvec = func_stage(svec,cppsi')

		// time_long_y = interp1(par.tvec_c_long,par.time_long,tvec,'linear',NaN); // Convert to time units
		CO_norm = cppsi[1,1]:*interp1(tvec,yC,svec)
		// Eliminate the -99999 for the computation of gamma, but that messes up the order, figure it out
		// I am not taking into account the half days, do it!
		// Trim on indetification window.
		tu1 = floor(tu_norm)
		C_del =  CO_norm[tu..tu1,1]
		T_del =  yT[tu..tu1,1]
		aux_del = rows(T_del)
		
		unr = kumsum(T_del) // rougth approximation for area under red
		unb = kumsum(C_del) // rougth approximation for area under blue
		
		// Policy Effect
		gamma_val = (unr[aux_del,1]-unb[aux_del,1])/unb[aux_del,1]
		
		// Identification window size
		olp = tu_norm-tu
		
		// Print output
		psi0 = strofreal(cppsi[1,1])
		psi1 = strofreal(cppsi[1,2])
		psi2 = strofreal(cppsi[1,3])
		gamma_s = strofreal(gamma_val[1,1])
		olp_s = strofreal(olp[1,1])
		
		s0 = ("psi0 =",psi0)
		s1 = ("psi1 =",psi1)
		s2 = ("psi2 =",psi2)
		s0 = invtokens(s0, " ")
		s1 = invtokens(s1, " ")
		s2 = invtokens(s2, " ")
		
		s3 = ("gamma =",gamma_s)
		s4 = ("olp =",olp_s)
		s3 = invtokens(s3, " ")
		s4 = invtokens(s4, " ")
		
		s0L = ("********")
		s0LL = ("******************************************")
		s1L = ("Results:")
		s2L = ("Normalization Coefficients:")
		printf("\n%s\n",s0L)
		printf("%s\n",s1L)
		printf("%s\n",s0L)
		printf("\n%s\n",s2L)
		printf("%s\n",s0)
		printf("%s\n",s1)
		printf("%s\n",s2)
		s3L = ("Policy Effect:")
		printf("\n%s\n",s3L)
		printf("%s\n",s3)
		s4L = ("Length of ident window:")
		printf("\n%s\n",s4L)
		printf("%s\n",s4)
}
end

