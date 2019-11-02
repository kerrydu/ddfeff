{smcl}
{* *! version 0.3  29 Oct 2019}{...}
{cmd:help nddfeff}
{hline}

{title:Title}

{phang}
{bf:nddfeff} {hline 2} Non-radial Directional Distance Function in Stata

{title:Syntax}

{p 8 21 2}
{cmd:nddfeff} {it:{help varlist:inputvars}} {cmd:=} {it:{help varlist:desirable_outputvars}} {cmd::} {it:{help varlist:undesirable_outputvars}} {ifin} 
{cmd:,} {cmdab:d:mu(}{varname}{cmd:)} [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{cmdab:d:mu:(varname)}}specifies names of DMUs. It is required. 

{synopt:{cmdab:t:ime:(varname)}}specifies time period for contemporaneous production technology. If {opt time:(varname)} is not specified, global production technology is assumed. 
{p_end}

{synopt:{opt wmat(name)}}specifies a weight matrix. 
The default is W=(1/(3*#.inputvars),..,1/(3*#.desirable_outputvars),...,1/(3*#.undesirable_outputvars)).
{p_end}

{synopt:{opt gx(varlist)}}specifies direction components for input adjustment. The default is gx=0.
{p_end}

{synopt:{opt gy(varlist)}}specifies direction components for desirable output adjustment. The default is gy=Y.
{p_end}

{synopt:{opt gb(varlist)}}specifies direction components for undesirable output adjustment. The default is gb=-B. 
{p_end}

{synopt:{cmdab:seq:uential}}specifies sequential production technology.
{p_end}

{synopt:{opt vrs}}specifies production technology with variable returns to scale. By default, production technology with constant returns to scale is assumed.
{p_end}

{synopt:{opt sav:ing(filename[,replace])}}specifies that the results be saved in {it:filename}.dta. 
{p_end}

{synopt:{opt maxiter(#)}}specifies the maximum number of iterations, which must be an integer greater than 0. The default value of maxiter is 16000.
{p_end}

{synopt:{opt tol(real)}}specifies the convergence-criterion tolerance, which must be greater than 0.  The default value of tol is 1e-8.
{p_end}

{synoptline}
{p2colreset}{...}

{title:Description}

{pstd}
{cmd:nddfeff} selects the input and output variables from the user designated data file or in the opened data set and solves directional distance function models by options specified. 

{phang}
The nddfeff program uses the buit-in mata function linearprogram(). Stata 16 or later is required.

{phang}
The nddfeff program requires initial data set that contains the input and output variables for observed units. 

{phang}
Variable names must be identified by inputvars for input variable, by desirable_outputvars for desirable output variable,  
and by undesirable_outputvars for undesirable output variable
 to allow that {cmd:nddfeff} program can identify and handle the multiple input-output data set. The direction vector
  g=(gx,gy,gb) should be specfied by options. The default is g=(-X,Y,-B).



{title:Examples}

{phang}{"use ...\example_ddf.dta"}

{phang}{cmd:. nddfeff labor capital energy= gdp: co2, dmu(id)}

{phang}{cmd:. nddfeff labor capital energy= gdp: co2, dmu(id) time(t) sav(nddf_result)}

{phang}{cmd:. nddfeff labor capital energy= gdp: co2, dmu(id) time(t) seq  sav(nddf_result,replace)}




{title:Saved Results}

{psee}
Macro:

{psee}
{cmd: r(file)} the stored results of {cmd:nddfeff} that have observation rows of DMUs and variable columns with input data, output data, values of DDF, and MLPI.
{p_end}


{marker references}{...}
{title:References}

{phang}
Chung, Y.H., Färe, R., Grosskopf, S. Productivity and Undesirable Outputs: A Directional Distance Function Approach.
 Journal of Environmental Management, 1997, 51:229-240.

{phang}
Oh, D.-h. A global Malmquist-Luenberger productivity index. Journal of Productivity Analysis, 2010, 34:183-197.

{phang}
Oh, D.-h., Heshmati A.  A sequential Malmquist–Luenberger productivity index: Environmentally sensitive productivity growth 
considering the progressive nature of technology. Energy Economics, 2010,3 2:1345-1355.

{phang}
Zhou, P., Ang, B.W., Wang H. Energy and CO2 emission performance in electricity generation: a non-radial directional distance function approach. Eur. J. Oper. Res., 2012, 221:625-635.

{title:Author}

{psee}
Kerry Du

{psee}
Xiamen University

{psee}
Xiamen, China

{psee}
E-mail: kerrydu@xmu.edu.cn
{p_end}
