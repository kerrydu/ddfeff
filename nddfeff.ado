
*! version 1.0
* By Kerry Du, 2 Nov 2019 
**
* 
capture program drop nddfeff
program define nddfeff, rclass
    version 16

    gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
		local invars `invars' `word'
		gettoken word 0 : 0, parse("=:,")
    }
    unab invars : `invars'

	gettoken word 0 : 0, parse(" =:,")
    while `"`word'"' != ":" & `"`word'"' != "=" {
        if `"`word'"' == "," | `"`word'"'=="" {
                error 198
        }
		local gopvars `gopvars' `word'
		gettoken word 0 : 0, parse(" =:,")
    }
    unab gopvars : `gopvars'
	
	
    syntax varlist [if] [in], dmu(varname) [gx(varlist) gy(varlist) gb(varlist)  ///
	                                       Time(varname) SEQuential VRS     ///
										   Wmat(string) SAVing(string)           ///
										   maxiter(numlist integer >0 max=1) tol(numlist max=1)]
	
	
	local techtype "contemporaneous production technology"
   
   if `"`time'"'==""{
   		local techtype "global production technology"

   }

	
	
	if "`maxiter'"==""{
		local maxiter=-1
	}
	if "`tol'"==""{
		local tol=-1
	}	
	
	if "`sequential'"!=""{
		if "`time'"==""{
		   disp as error "For sequential model, time() should be specified."
		   error 498
		}
		else{
		   local techflag "<="
		   local techtype "sequential production technology"
		}
	
	}
	


	
	preserve
	marksample touse 
	markout `touse' `invars' `gopvars' `gx' `gy' `gb'

	local bopvars `varlist'
	
	local invars: list uniq invars
	local gopvars: list uniq gopvars
	local bopvars: list uniq bopvars
	
	local ninp: word count `invars'
    local ngo: word count `gopvars'
    local nbo: word count `bopvars'
	local nvar=`ninp'+`ngo'+`nbo'
	
	confirm numeric var `invars' `gopvars' `bopvars'
	
	
	if "`gx'"!=""{
		local ngx: word count `gx'
		if `ngx'!=`ninp'{
		    disp as error "# of input variables != # of variables specified in gx()."
		    error 498
		}
		local gmat `gmat' `gx'
	
	}
	else{
		local invarscopy `invars' 
		forv k=1/`ninp'{
			gettoken word invarscopy:invarscopy
			//di "`word'"
		    tempvar gx_`k'
			qui gen `gx_`k''=-`word'
			local gmat `gmat' `gx_`k''
		}
	
	}
	
	if "`gy'"!=""{
		local ngy: word count `gy'
		if `ngy'!=`ngo'{
		    disp as error "# of desriable output variables != # of variables specified in gy()."
		    error 498
		}
		local gmat `gmat' `gy'
	
	}
	else{
	    local gopvarscopy `gopvars'
		forv k=1/`ngo'{
		    gettoken word gopvarscopy:gopvarscopy
		    tempvar gy_`k'
			qui gen `gy_`k''=`word'
			local gmat `gmat' `gy_`k''
		}
	
	}
		
	if "`gb'"!=""{
		local ngb: word count `gb'
		if `ngb'!=`nbo'{
		    disp as error "# of undesriable output variables != # of variables specified in gb()."
		    error 498
		}
		local gmat `gmat' `gb'
	
	}
	else{
	    local bopvarscopy `bopvars'
		forv k=1/`nbo'{
		    gettoken word bopvarscopy:bopvarscopy
		    tempvar gb_`k'
			qui gen `gb_`k''=-`word'
			local gmat `gmat' `gb_`k''
		}
	
	}	
	
	tempname weightvec
	
	if `"`wmat'"'!=""{
		confirm matrix `wmat'
		local ncol=colsof(`wmat')
		if `ncol'!=`nvar'{
		    dis as error `"# of column of `wmat' != # of input-output variables"'
			exit 498
		}
		mat `weightvec'=`wmat'
		
		
	   }
	else{
		
		mat  `weightvec'=(J(1,`ninp',1)*(1/3/`ninp'),J(1,`ngo',1)*(1/3/`ngo'),J(1,`nbo',1)*(1/3/`nbo'))
			
	}
	
	
	local comvars: list invars & gopvars 
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as input and desriable output simultaneously."
		error 498
	}
	
	local comvars: list invars & bopvars
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as input and undesriable output simultaneously."
		error 498
	}	
	
	local comvars: list gopvars & bopvars
	if !(`"`comvars'"'==""){
		disp as error "`comvars' should not be specified as desriable and undesriable outputs simultaneously."
		error 498
	}	
		
	

	
	local rstype=1
	if "`vrs'"!=""{
	   local rstype=0
	}
	

	

	qui keep   `invars' `gopvars' `bopvars' `dmu' `time' `gmat' `touse'
	qui gen _Row=_n
	qui keep if `touse'
	label var _Row "Row #"
	qui gen double Dval=.
	label var Dval "Value of NDDFs: `techtype'"
	

	
    tempvar tvar dmu2
	

	
	if  `"`time'"'!="" {
	    qui egen `tvar'=group(`time')
	    qui egen `dmu2'=group(`dmu')
	}
	else{
	    qui gen `tvar'=1
		qui gen `dmu2'=_n
	}
	
	sort `dmu2' `tvar' _Row

	foreach v in `invars' `gopvars' `bopvars'{
	   qui gen double B_`v'=.
	   label var B_`v' `"beta:`v'"'
	   local slackvars `slackvars' B_`v'
	
	}
		
	
	qui mata mata mlib index
	mata: _nDDFs(`"`invars'"',`"`gopvars'"',`"`bopvars'"',"`dmu2'","`tvar'", ///
	               "`gmat'","`weightvec'","`techflag'",`rstype',`maxiter',`tol',"Dval",`"`slackvars'"')
	


	order _Row `dmu' `time' Dval  `slackvars'
	keep  _Row `dmu' `time' Dval  `slackvars'
	
	disp _n(2) " Non-raidal Directional Distance Function (NDDF) Results:"
	disp "    (_Row: Row # in the original data; Dval: Estimated value of NDDF.)"
	list, sep(0) 
	di "Note: missing value indicates infeasible problem."
	
	if `"`saving'"'!=""{
	  save `saving'
	  gettoken filenames saving:saving, parse(",")
	  local filenames `filenames'.dta
	  disp _n `"Estimated Results are saved in `filenames'."'
	}

	return local file `filenames'

	
	restore 
	
	end
	


/*	
make ddfeff, replace toc pkg title(Directional Distance Function for Efficiency/Productivity Analysis) ///
             version(1.0) author(Kerry Du) affiliation(Xiamen University) ///
			 email(kerrydu@xmu.edu.cn) install("ddfeff.ado;ddfeff.sthlp;nddfeff.ado;nddfeff.sthlp;lddfeff.mlib")
*/

