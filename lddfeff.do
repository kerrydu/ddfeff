
version 16

cap mata mata drop _DDFmain()
cap mata mata drop trans2()
cap mata mata drop crossDDFs()
cap mata mata drop DDFs()
cap mata mata drop ddf_crs()
cap mata mata drop ddf_vrs()
cap mata mata drop relflag()
cap mata mata drop nddf_crs()
cap mata mata drop nddf_vrs()
cap mata mata drop _nDDFs()
cap mata mata drop _lpres()


mata:
mata clear 

	struct lpres   { real scalar    fval 
					 real matrix    coeff
					 real scalar    converged
					 real scalar    returncode
					  }


 void  function _nDDFs( string scalar     Xvar,  ///
				        string scalar     Yvar,  ///
				        string scalar     Bvar,  ///
						string scalar     idvar, ///
						string scalar     tvar,  ///
						string scalar     gmatname, ///
						string scalar     wmatname, ///						
						string scalar     rel, ///
						real scalar       rstype, ///
						real scalar       maxiter, ///
						real scalar       tol,  ///
						string scalar     vScore, ///
						string scalar     vBeta)
{
	
	
	struct lpres scalar nddfres
	
	X=st_data(.,Xvar)
	Y=st_data(.,Yvar)
	B=st_data(.,Bvar)
	id=st_data(.,idvar)
	t=st_data(.,tvar)
	gmat=st_data(.,gmatname)
	wmat=st_matrix(wmatname)

	nv=cols(X)+cols(Y)+cols(B)
    id2= uniqrows(id)
	t2 = uniqrows(t)
	k  = 1
	fval=J(rows(X),1,.)
	betas=J(rows(X),nv,.)
	for(i=1;i<=length(id2);i++){
		for(j=1;j<=length(t2);j++){
		   //(id:==id2[i]):&(t:==t2[j])
		   flag=(id:==id2[i])+(t:==t2[j])
		   XX   = select(X,flag:==2)
		   YY   = select(Y,flag:==2)
		   BB   = select(B,flag:==2)
		   g    = select(gmat,flag:==2)
		   flag=relflag(rel,t, t2[j])
		   Xref = select(X,flag)
		   Yref = select(Y,flag)
		   Bref = select(B,flag)
		   if(rstype==1){
		     nddfres=nddf_crs(XX,YY,BB,Xref,Yref,Bref,g,wmat,maxiter,tol)		   
		   }
		   else{
             nddfres=nddf_vrs(XX,YY,BB,Xref,Yref,Bref,g,wmat,maxiter,tol)		   
		   }

		   if(nddfres.converged==1){
			   fval[k]=nddfres.fval	
			   betas[k,.]=nddfres.coeff[1,1..length(wmat)]		   
		   }

		   k=k+1
		}
	
	}
	
	
	st_view(Score=.,.,vScore)
	Score[.,.]=fval
	st_view(SSbeta=.,.,vBeta)
	SSbeta[.,.]=betas
		
	
	
	

}


			  

function nddf_crs(  real rowvector    X, ///
			        real rowvector    Y,  ///
			        real rowvector    B, ///
			        real matrix    Xref, ///
			        real matrix    Yref, ///
			        real matrix    Bref, ///
			        real rowvector    g, ///
					real rowvector    wmat, /// 
					real scalar       maxiter, ///
					real scalar       tol)
	{
	
	    class LinearProgram scalar q
		
		nx=length(X)
		ny=length(Y)
		nb=length(B)
		gX=g[1..nx]
		gY=g[(nx+1)..(nx+ny)]
		gB=g[(nx+ny+1)..(nx+ny+nb)]
		c=(wmat,J(1,rows(Xref),0))
		lowerbd=J(1,length(wmat),0),J(1,rows(Xref),0)
		upperbd=J(1,length(c),.)		
		q = LinearProgram()
		q.setCoefficients(c)
		q.setBounds(lowerbd, upperbd)
		if(maxiter!=-1){
		  q.setMaxiter(maxiter)
		}
		if (tol!=-1){
		  q.setTol(tol)
		}	
		//wmat	
		 Aie1=diag(-gX'),J(nx,ny+nb,0),Xref'
		 bie1=X'
		 Aie2=J(ny,nx,0),diag(gY'),J(ny,nb,0),-Yref'
		 bie2=-Y'
		 Aec=J(nb,nx+ny,0),diag(-gB'),Bref'
		 bec=B'
		 q.setEquality(Aec,bec)
		 q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))
		 beta=q.optimize()
		 struct lpres scalar retres
		 retres.fval=q.value()
		 retres.coeff=q.parameters()
         retres.converged=q.converged()
		 retres.returncode=q.returncode()
		 return(retres)	
		

}



/////////////////////////////////////////////////////////////////////
    function nddf_vrs(  real rowvector    X, ///
				        real rowvector    Y,  ///
				        real rowvector    B, ///
				        real matrix    Xref, ///
				        real matrix    Yref, ///
				        real matrix    Bref, ///
				        real rowvector    g, ///
						real rowvector    wmat, /// 
						real scalar       maxiter, ///
						real scalar       tol)
	{
	
	    class LinearProgram scalar q
		
		nx=length(X)
		ny=length(Y)
		nb=length(B)
		gX=g[1..nx]
		gY=g[(nx+1)..(nx+ny)]
		gB=g[(nx+ny+1)..(nx+ny+nb)]
		c=(wmat,J(1,rows(Xref),0))
		lowerbd=J(1,length(wmat),0),J(1,rows(Xref),0)
		upperbd=J(1,length(c),.)		
		q = LinearProgram()
		q.setCoefficients(c)
		q.setBounds(lowerbd, upperbd)
		if(maxiter!=-1){
		  q.setMaxiter(maxiter)
		}
		if (tol!=-1){
		  q.setTol(tol)
		}			
		 Aie1=diag(-gX'),J(nx,ny+nb,0),Xref'
		 bie1=X'
		 Aie2=J(ny,nx,0),diag(gY'),J(ny,nb,0),-Yref'
		 bie2=-Y'
		 Aec=J(nb,nx+ny,0),diag(-gB'),Bref'
		 bec=B'
		 lsum=J(1,length(wmat),0),J(1,rows(Xref),1)
		 q.setEquality(Aec \ lsum,bec \1)
		 q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))
		 beta=q.optimize()
		 struct lpres scalar retres
		 retres.fval=q.value()
		 retres.coeff=q.parameters()
         retres.converged=q.converged()
		 retres.returncode=q.returncode()
		 return(retres)	
		

}







////////////////////////////////////////////////////////////////////
void function _DDFmain(         string scalar       Xname,  ///
						        string scalar       Yname,  ///
						        string scalar       Bname,  ///
								string scalar       idname, ///
								string scalar       tname,  ///
								string scalar     gmatname, /// 
								real scalar       gpt, ///
								real scalar       rstype, ///
								string scalar     rel, ///
								string scalar     Dvalname, ///
								string scalar     mqlname, ///
								string scalar     techname, ///
								string scalar     tecchname, ///
								real scalar       maxiter, ///
								real scalar       tol) 

{

	gmat=st_data(.,gmatname)
	X=st_data(.,Xname)
	Y=st_data(.,Yname)
	B=st_data(.,Bname)
	id=st_data(.,idname)
	t=st_data(.,tname)
    N=length(t)
	//idpos1= (id[2::N]-id[1::(N-1)])	\ 1
	idpos2=1 \ (id[2::N]-id[1::(N-1)])
	retpos=select(1::N,idpos2:==0)
	st_view(Dvalue=.,.,Dvalname)
    if(gpt==1){
		newt=J(N,1,1)
		newid=1::N
		tflag=J(N,1,1)
		Dvalue[.,.]=DDFs(X,Y,B,newid,newt,tflag,newt,gmat,rel,rstype,maxiter,tol)
		if(mqlname!=""){
		   st_view(mql=.,.,mqlname)
		   DD=trans2(id,Dvalue)
	       mql[retpos,.]=(1:+DD[.,1]):/(1:+DD[.,2])

		}

	}
	else{
		//gpt
		tflag=J(N,1,1)
		Dvalue[.,.]=DDFs(X,Y,B,id,t,tflag,t,gmat,rel,rstype,maxiter,tol)		
	
		if(mqlname!=""){
		   st_view(tech=.,.,techname)
		   st_view(tecch=.,.,tecchname)
		   DD=trans2(id,Dvalue)
		   tech[retpos,.]=(1:+DD[.,1]):/(1:+DD[.,2])
		   //length(tech)
		   crD= crossDDFs(X,Y,B,id,t,gmat,rel,rstype,maxiter,tol)
	       tp=(1:+crD[.,2]):*(1:+DD[.,2])
		   tp=tp:/(1:+DD[.,1])
		   tp=tp:/(1:+crD[.,1])
		   tecch[retpos,.]=tp:^0.5
		
		   st_view(mql=.,.,mqlname)
		   mql[.,.]=tech:*tecch

		   
		}	
	
	
	}



}






real matrix function trans2(real colvector id, real colvector beta) 
{
   N=length(id)
   idpos1= (id[2::N]-id[1::(N-1)]) \ 1 
   idpos2=1 \ (id[2::N]-id[1::(N-1)])
   conD=select(beta,idpos1:==0),select(beta,idpos2:==0)
   return(conD)
}

   






real matrix function crossDDFs( real matrix       X,  ///
						        real matrix       Y,  ///
						        real matrix       B,  ///
								real colvector    id, ///
								real colvector    t,  ///
								real matrix       gmat, /// 
								string scalar     rel, ///
								real scalar       rstype, ///
								real scalar       maxiter, ///
								real scalar       tol) 
{


   tflag=(t:<max(t))
   tau0=select(t,t:>min(t))
   //id0=select(id,t:<max(t))
   D21=DDFs(X,Y,B,id,t,tflag,tau0,gmat,rel,rstype,maxiter,tol)
   tflag=(t:>min(t))
   tau0=select(t,t:<max(t))
   D12=DDFs(X,Y,B,id,t,tflag,tau0,gmat,rel,rstype,maxiter,tol)  
   crossD=D12,D21
   return(crossD)
  							
								
}




real matrix function DDFs(      real matrix       X,  ///
						        real matrix       Y,  ///
						        real matrix       B,  ///
								real colvector    id, ///
								real colvector    t,  ///
								real colvector    tflag, ///
								real colvector    tau, ///
								real matrix       gmat, ///
								string scalar     rel, ///
								real scalar       rstype, ///
								real scalar       maxiter, ///
								real scalar       tol)
	{

    id2= uniqrows(id)
	t2 = uniqrows(select(t,tflag))
	tau2= uniqrows(tau)
	//t2,tau2
	k  = 1
	//BETA=J(0,1,.)
	BETA=J(length(select(t,tflag)),1,.)
	for(i=1;i<=length(id2);i++){
		for(j=1;j<=length(tau2);j++){
		   //(id:==id2[i]):&(t:==t2[j])
		   flag=(id:==id2[i])+(t:==t2[j])
		   XX   = select(X,flag:==2)
		   YY   = select(Y,flag:==2)
		   BB   = select(B,flag:==2)
		   g    = select(gmat,flag:==2)
		   flag = relflag(rel,t, tau2[j])
		   Xref = select(X,flag)
		   Yref = select(Y,flag)
		   Bref = select(B,flag)
		   if(rstype==1){
		     BETA[k]=ddf_crs(XX,YY,BB,Xref,Yref,Bref,g,maxiter,tol)		   
		   }
		   else{
             BETA[k]=ddf_vrs(XX,YY,BB,Xref,Yref,Bref,g,maxiter,tol)		   
		   }

		   //BETA= BETA \ btemp
		   //length(BETA)

		   k=k+1
		}
	
	}
	
	return(BETA)
	

}




real scalar function ddf_crs(   real rowvector    X, ///
						        real rowvector    Y,  ///
						        real rowvector    B, ///
						        real matrix    Xref, ///
						        real matrix    Yref, ///
						        real matrix    Bref, ///
						        real rowvector    g, ///
								real scalar       maxiter, ///
								real scalar       tol)
	{
	
	    class LinearProgram scalar q
		
		nx=length(X)
		ny=length(Y)
		nb=length(B)
		gX=g[1..nx]
		gY=g[(nx+1)..(nx+ny)]
		gB=g[(nx+ny+1)..(nx+ny+nb)]
		c=(1,J(1,rows(Xref),0))
		lowerbd=.,J(1,length(c)-1,0)
		upperbd=J(1,length(c),.)		
		q = LinearProgram()
		q.setCoefficients(c)
		q.setBounds(lowerbd, upperbd)
		if(maxiter!=-1){
		  q.setMaxiter(maxiter)
		}
		if (tol!=-1){
		  q.setTol(tol)
		}			
		 Aie1=-gX',Xref'
		 bie1=X'
		 Aie2=gY',-Yref'
		 bie2=-Y'
		 Aec=-gB',Bref'
		 bec=B'
		 q.setEquality(Aec,bec)
		 q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))
		 beta=q.optimize()		
		
		
	     return(beta)

}




real scalar function ddf_vrs(   real rowvector    X, ///
						        real rowvector    Y,  ///
						        real rowvector    B, ///
						        real matrix    Xref, ///
						        real matrix    Yref, ///
						        real matrix    Bref, ///
						        real rowvector    g, ///
								real scalar       maxiter, ///
								real scalar       tol) 
	{
	
	    class LinearProgram scalar q
		
		nx=length(X)
		ny=length(Y)
		nb=length(B)
		gX=g[1..nx]
		gY=g[(nx+1)..(nx+ny)]
		gB=g[(nx+ny+1)..(nx+ny+nb)]
		c=(1,J(1,rows(Xref),0))
		lowerbd=.,J(1,length(c)-1,0)
		upperbd=J(1,length(c),.)		
		q = LinearProgram()
		q.setCoefficients(c)
		q.setBounds(lowerbd, upperbd)
		if(maxiter!=-1){
		  q.setMaxiter(maxiter)
		}
		if (tol!=-1){
		  q.setTol(tol)
		}			
		 Aie1=-gX',Xref'
		 bie1=X'
		 Aie2=gY',-Yref'
		 bie2=-Y'
		 Aec=-gB',Bref'
		 bec=B'
		 lsum=0,J(1,rows(Xref),1)
		 q.setEquality(Aec \ lsum, bec \ 1)
		 q.setInequality((Aie1 \ Aie2 ), (bie1 \ bie2))
		 beta=q.optimize()		
		
		
	     return(beta)

}





real matrix function relflag(string scalar rel, real colvector x, real scalar y)

{
   if(rel=="<="){
      flag=(x:<=y)
   }
   else{
      flag=(x:==y)
   }
   return(flag)

}

mata mlib create lddfeff, replace 
mata mlib add lddfeff *()
mata mlib index

end

