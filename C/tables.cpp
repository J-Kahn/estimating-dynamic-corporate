#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "tableout.h"

void print_struct_est(string filename, string suffixm, string suffixq)
{
  
  ofstream fileout;
  fileout.open(filename.c_str());
  
  matrix par1 = readcsv("../Results/Estim/par_" + suffixm + ".csv",7,1);
  matrix par2 = readcsv("../Results/Estim/par_" + suffixq + ".csv",7,1);
  matrix jstatm11 = readcsv("../Results/Estim/j_mom_" + suffixm + ".csv",1,1);
  matrix jstatm12 = readcsv("../Results/Estim/j_epfq_" + suffixq + ".csv",1,1);
  matrix jstatm21 = readcsv("../Results/Estim/j_epfq_" + suffixm + ".csv",1,1);
  matrix jstatm22 = readcsv("../Results/Estim/j_mom_" + suffixq + ".csv",1,1);
  matrix vparm1 = readcsv("../Results/Estim/var_par_" + suffixm + ".csv",7,7);
  matrix vparm2 = readcsv("../Results/Estim/var_par_" + suffixq + ".csv",7,7);
  matrix vpar1 = sqrtdiags(vparm1);
  matrix vpar2 = sqrtdiags(vparm2);
  double jstat11 = jstatm11(0);
  double jstat12 = jstatm12(0);
  double jstat21 = jstatm21(0);
  double jstat22 = jstatm22(0);
  
  fileout << "\\begin{table}     \\caption{Structural Estimates of Trade-Off model " << endl;
  fileout << "Parameters with Different Benchmarks \\label{StructuralEstimationTable}}" << endl;
  fileout << "\\begin{center}" << endl;
  fileout << "\\begin{tabular*}{1\\textwidth}{l@{\\extracolsep{\\fill}}cc}" << endl;
  fileout << "Panel A: Parameters    \\\\" << endl;
  fileout << "\\hline" << endl;
  fileout << "Parameter  &   TM     & EPF \\\\ " << endl;
  fileout << "\\hline" << endl;
  fileout << setprecision(3) << "$\\mu$         &    " <<  par1(4) << " &   " <<  par2(4) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(4) << ")   & (" << vpar2(4) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "$\\rho$        &    " <<  par1(5) << " &   " <<  par2(5) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(5) << ")   & (" << vpar2(5) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "$\\sigma$      &    " <<  par1(6) << " &   " <<  par2(6) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(6) << ")   & (" << vpar2(6) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "$\\delta$      &    " <<  par1(0) << " &   " <<  par2(0) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(0) << ")   & (" << vpar2(0) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "$\\psi$        &    " <<  par1(2) << " &   " <<  par2(2) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(2) << ")   & (" << vpar2(2) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "$s$            &    " <<  par1(3) << " &   " <<  par2(3) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(3) << ")   & (" << vpar2(3) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "$\\lambda$     &    " <<  par1(1) << " &   " <<  par2(1) << "      \\\\ " << endl;
  fileout << setprecision(3) << "               &    (" << vpar1(1) << ")   & (" << vpar2(1) << ")         \\\\ " << endl;
  fileout << setprecision(3) << "Overidentifying       &  " << jstat11  << " & " << jstat12 << "     \\\\" << endl;
  fileout << setprecision(3) << "restrictions $\\chi^2$ &    &      \\\\ " << endl;
  fileout << setprecision(3) << "Out-of-sample         &   " <<  jstat21   << " & " << jstat22 << "    \\\\ " << endl;
  fileout << setprecision(3) << "$\\chi^2$              &     &     \\\\ " << endl;
  fileout << "\\hline " << endl;
  fileout << "\\\\ " << endl;
  fileout << "\\end{tabular*} " << endl;
  
  fileout.close();
}

void print_mom_est(string filename, string suffixm, string suffixq)
{
  
  ofstream fileout;
  fileout.open(filename.c_str());

  matrix mom = readcsv("../Data/mom_comp.csv",8,1);
  matrix mom1 = readcsv("../Results/Estim/mom_" + suffixm + ".csv",8,1);
  matrix mom2 = readcsv("../Results/Estim/mom_" + suffixq + ".csv",8,1);
  matrix vmmom1 = readcsv("../Results/Estim/var_mom_" + suffixm + ".csv",8,8);
  matrix vmmom2 = readcsv("../Results/Estim/var_mom_" + suffixq + ".csv",8,8);
  matrix vmom1 = sqrtdiags(vmmom1);
  matrix vmom2 = sqrtdiags(vmmom2);
  matrix tmom1 = mom1/vmom1;
  matrix tmom2 = mom2/vmom2;
  
fileout << "\\begin{tabular*}{1\\textwidth}{l@{\\extracolsep{\\fill}}ccccccc} " << endl;
fileout << "%\\begin{tabular}{lccccccc} " << endl;
fileout << "Panel B: Moments \\\\ " << endl;
fileout << "\\hline " << endl;
fileout << "                & Data &&\\multicolumn{2}{c}{Traditional moments} && \\multicolumn{2}{c}{EPF polynomial} \\\\ " << endl;
fileout << "                &      && Simulated & {\\em t}-statistic  && Simulated & {\\em t}-statistic  \\\\ " << endl;
fileout << "\\hline \\\\[2pt] " << endl;
fileout << setprecision(3) << "Mean leverage          &  " << mom(0) << " && " << mom1(0) << " &  " << tmom1(0) << "   &&  " << mom2(0) << " &  " << tmom2(0) << "  \\\\ " << endl;
fileout << setprecision(3) << "Variance leverage      &  " << mom(1) << " && " << mom1(1) << " &  " << tmom1(1) << "   &&  " << mom2(1) << " &  " << tmom2(1) << "  \\\\ " << endl;
fileout << setprecision(3) << "Mean investment        &  " << mom(2) << " && " << mom1(2) << " &  " << tmom1(2) << "   &&  " << mom2(2) << " &  " << tmom2(2) << "  \\\\ " << endl;
fileout << setprecision(3) << "Variance investment    &  " << mom(3) << " && " << mom1(3) << " &  " << tmom1(3) << "   &&  " << mom2(3) << " &  " << tmom2(3) << "  \\\\ " << endl;
fileout << setprecision(3) << "Mean distributions     &  " << mom(4) << " && " << mom1(4) << " &  " << tmom1(4) << "   &&  " << mom2(4) << " &  " << tmom2(4) << "  \\\\ " << endl;
fileout << setprecision(3) << "Variance distributions &  " << mom(5) << " && " << mom1(5) << " &  " << tmom1(5) << "   &&  " << mom2(5) << " &  " << tmom2(5) << "  \\\\ " << endl;
fileout << setprecision(3) << "Mean profits           &  " << mom(6) << " && " << mom1(6) << " &  " << tmom1(6) << "   &&  " << mom2(6) << " &  " << tmom2(6) << "  \\\\ " << endl;
fileout << setprecision(3) << "Variance profits       &  " << mom(7) << " && " << mom1(7) << " &  " << tmom1(7) << "   &&  " << mom2(7) << " &  " << tmom2(7) << "  \\\\ " << endl;
fileout << "\\hline " << endl;
fileout << "%\\end{tabular} " << endl;
fileout << "\\end{tabular*} " << endl;
fileout << "\\parbox[c]{6.5in}{\\footnotesize Panel A contains parameter estimates parameters from a trade-off " << endl;
fileout << "model using two differen benchmarks: traditional moments (TM) and empirical policy functions (EPF). Indirect inference is performed by " << endl;
fileout << "minimizing the (inverse covariance matrix weighted) distance between the " << endl;
fileout << "simulated values of each set of benchmarks and the corresponding values " << endl;
fileout << "found in Compustat. " << endl;
fileout << "$\\delta$ is the depreciation rate of capital, $\\lambda$ is the linear " << endl;
fileout << "equity issuance cost, $s$ is the collateral value of capital, and $\\psi$ " << endl;
fileout << "is the quadratic adjustment cost.  The parameters in the second pane are " << endl;
fileout << "estimated in a seperate stage from data on productivity and leverage. " << endl;
fileout << "$\\alpha$ is the curvature of the production function, $\\rho$ is the " << endl;
fileout << "persistence of the shock to productivity, and $\\sigma$ is the variance " << endl;
fileout << "of the autoregressive error in productivity. Clustered standard errors are in parentheses under the parameter estimates. In Panel B," << endl;
fileout << "the first column present the data estimates used in the traditional moments (TM) estimation. The second and third contain the simulated moments from the traditional moments estimation, along" << endl;
fileout << "with the {\\em t}-statistic for the difference between the actual and simulated moments. The next two pairs of columns contain analogous results for the two empirical policy function (EPF) estimations.}" << endl;

fileout << "\\end{center} " << endl;
fileout << "\\end{table} " << endl;
  
  fileout.close();
}



void print_struct_mc(string filename, matrix diff1, matrix p1, matrix parm1, matrix diff2, matrix p2, matrix parm2, matrix diff3, matrix p3, matrix parm3, matrix diff4, matrix p4, matrix parm4)
{
  ofstream fileout;
  fileout.open(filename.c_str());
  
  matrix bias1 = cmean(diff1), bias2 = cmean(diff2), bias3 = cmean(diff3), bias4 = cmean(diff4);
  matrix rmse1 = csqrtmeansq(diff1), rmse2 = csqrtmeansq(diff2), rmse3 = csqrtmeansq(diff3), rmse4 = csqrtmeansq(diff4);
  matrix pt1 = cmean(p1), pt2 = cmean(p2), pt3 = cmean(p3), pt4 = cmean(p4);
  matrix par1 = cmean(parm1), par2 = cmean(parm2), par3 = cmean(parm3), par4 = cmean(parm4);
  
  fileout << "\\begin{table}[hhh] " << endl;
  fileout << "  \\caption{Monte Carlo Comparison of Simulated Moments Estimators with Different Matching Objectives\\label{tab:parameters_small}} " << endl;
  fileout << "  \\begin{center} " << endl;
  fileout << "  \\begin{tabular*}{1\\textwidth}{r@{\\extracolsep{\\fill}}ccccc} " << endl;
  fileout << "    \\hline \\hline " << endl;
  fileout << "            &  \\multicolumn{2}{c}{Traditional Moments} & &\\multicolumn{2}{c}{Empirical Policy Functions} \\\\ " << endl;
  fileout << "            &  Diagonal  & Clustered &$\\qquad$ & Diagonal & Clustered \\\\ " << endl;
  fileout << "\\hline " << endl;
  fileout << "\\multicolumn{1}{l}{$\\delta$ (depreciation rate)} \\\\ " << endl;
  fileout << "   Mean Estimate       & " << par1(0) << " & " <<  par2(0) <<  " & " <<  par3(0) << " & " <<  par4(0) << "                     \\\\ " << endl;
  fileout << "   Average $\\%$ Bias & " << bias1(0) << " & " << bias2(0)  << " & " << bias3(0) << " & " << bias4(0) << "                    \\\\ " << endl;
  fileout << "   RMSE $\\%$         & " << rmse1(0) << " & " << rmse2(0) << " & " <<  rmse3(0) << " & " << rmse4(0) << "                        \\\\ " << endl;
  fileout << "   $\\Pr(t)$          & " <<   pt1(0) << " & " <<   pt2(0) << " & " <<    pt3(0) << " & " <<   pt4(0) << "                     \\\\  [10pt] " << endl;
  fileout << "\\multicolumn{1}{l}{$\\lambda$ (equity issuance cost)} \\\\ " << endl;
  fileout << "   Mean Estimate       & " << par1(1) << " & " <<  par2(1) <<  " & " <<  par3(1) << " & " <<  par4(1) << "                   \\\\ " << endl;
  fileout << "   Average $\\%$ Bias & " << bias1(1) << " & " << bias2(1)  << " & " << bias3(1) << " & " << bias4(1) << "                    \\\\ " << endl;
  fileout << "   RMSE $\\%$         & " << rmse1(1) << " & " << rmse2(1) << " & " <<  rmse3(1) << " & " << rmse4(1) << "                        \\\\ " << endl;
  fileout << "   $\\Pr(t)$          & " <<   pt1(1) << " & " <<   pt2(1) << " & " <<    pt3(1) << " & " <<   pt4(1) << "                     \\\\  [10pt] " << endl;
  fileout << "\\multicolumn{1}{l}{$s$ (collateral parameter)} \\\\ " << endl;
  fileout << "   Mean Estimate       & " << par1(3) << " & " <<  par2(3) <<  " & " <<  par3(3) << " & " <<  par4(3) << "                   \\\\ " << endl;
  fileout << "   Average $\\%$ Bias & " << bias1(3) << " & " << bias2(3)  << " & " << bias3(3) << " & " << bias4(3) << "                    \\\\ " << endl;
  fileout << "   RMSE $\\%$         & " << rmse1(3) << " & " << rmse2(3) << " & " <<  rmse3(3) << " & " << rmse4(3) << "                        \\\\ " << endl;
  fileout << "   $\\Pr(t)$          & " <<   pt1(3) << " & " <<   pt2(3) << " & " <<    pt3(3) << " & " <<   pt4(3) << "                     \\\\  [10pt] " << endl;
  fileout << "\\multicolumn{1}{l}{$\\delta$ (investment adjustment cost)} \\\\ " << endl;
  fileout << "   Mean Estimate       & " << par1(2) << " & " <<  par2(2) <<  " & " <<  par3(2) << " & " <<  par4(2) << "                   \\\\ " << endl;
  fileout << "   Average $\\%$ Bias & " << bias1(2) << " & " << bias2(2)  << " & " << bias3(2) << " & " << bias4(2) << "                    \\\\ " << endl;
  fileout << "   RMSE $\\%$         & " << rmse1(2) << " & " << rmse2(2) << " & " <<  rmse3(2) << " & " << rmse4(2) << "                        \\\\ " << endl;
  fileout << "   $\\Pr(t)$          & " <<   pt1(2) << " & " <<   pt2(2) << " & " <<    pt3(2) << " & " <<   pt4(2) << "                     \\\\  [10pt] " << endl;
  fileout << "\\hline " << endl;
  fileout << "\\end{tabular*} " << endl;
  fileout << "\\parbox[c]{6.5in}{Indicated expectations and probabilities are estimates " << endl;
  fileout << "based on 100 Monte Carlo samples of size 75,000.  The samples are " << endl;
  fileout << "generated from the model in Section \\ref{sec:model}.  Bias is " << endl;
  fileout << "expressed as a percent of the true coefficient value.  RMSE indicates " << endl;
  fileout << "root mean squared error and is also expressed as a percent of the true " << endl;
  fileout << "coefficient.$\\Pr(t)$ is the fraction of the time we observe a rejection of the null hypothesis that a parameter equals  its true value using " << endl;
  fileout << "a {\\em t}-test." << endl;
  fileout << "}" << endl;
  fileout << "\\end{center}" << endl;
  fileout << "\\end{table}" << endl;
  
  fileout.close();
}

/*void print_reject_mc(matrix p1, matrix p2, matrix p3, matrix p4)
{
  
  ofstream fileout;
  fileout.open(filename.c_str());
  
  
  
  fileout << " \\begin{table} " << endl;
  fileout << "   \\caption{Monte Carlo Comparison of Specification Tests}\\label{tab:Jstats} " << endl;
  fileout << " \\begin{center} " << endl;
  fileout << " % \\begin{tabular*}{1\\textwidth}{r@{\\extracolsep{\\fill}}ccccc} " << endl;
  fileout << "   \\begin{tabular}{lccccc} " << endl;
  fileout << "             &  \\multicolumn{2}{c}{Traditional Moments} & &\\multicolumn{2}{c}{Empirical Policy Functions} \\\\ " << endl;
  fileout << "             &  Diagonal  & Clustered &$\\qquad$ & Diagonal & Clustered \\\\ " << endl;
  fileout << "     \\hline \\hline " << endl;
  fileout << " \\multicolumn{1}{l}{Panel A: Minimized objective function} \\\\  \\\\ " << endl;
  fileout << "  Sample size = ???                                                                                     \\\\ " << endl;
  fileout << "     Overidentification test              & " <<  jst11[ ] << " & " <<  jst12[ ] << " & " <<  jst13[ ] << " & " << jst14[ ] << "\\\\ " << endl;
  fileout << "     Moment {\\em t}-statistics:                                                                         \\\\ " << endl;
  fileout << "       \\mbox{  } maximum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } median  rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } minimum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ \\\\ " << endl;
  fileout << "  Sample size = ???                                                                                     \\\\ " << endl;
  fileout << "     Overidentification test              & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "     Moment {\\em t}-statistics:                                                                         \\\\ " << endl;
  fileout << "       \\mbox{  } maximum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } median  rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } minimum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ \\\\ " << endl;
  fileout << " \\multicolumn{1}{l}{Panel B: True-parameter objective function}  &   &  &  &                \\\\     \\\\ " << endl;
  fileout << "  Sample size = ???                                                                                     \\\\ " << endl;
  fileout << "     Overidentification test              & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "     Moment {\\em t}-statistics:                                                                         \\\\ " << endl;
  fileout << "       \\mbox{  } maximum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } median  rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } minimum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\  \\\\ " << endl;
  fileout << "  Sample size = ???                                                                                    \\\\ " << endl;
  fileout << "     Overidentification test              & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "     Moment {\\em t}-statistics:                                                                        \\\\ " << endl;
  fileout << "       \\mbox{  } maximum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } median  rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "       \\mbox{  } minimum rejection rate  & " << maxt11[ ] << " & " << maxt12[ ] << " & " << maxt13[ ] << " & " << maxt14[ ] << "\\\\ " << endl;
  fileout << "     \\hline \\hline " << endl;
  fileout << "   \\end{tabular} " << endl;
  fileout << "   \\medskip " << endl;
  
  fileout << " \\parbox{5.2in}{Indicated expectations and probabilities are estimates " << endl;
  fileout << " based on 100 Monte Carlo samples of size 1,000.  The samples are " << endl;
  fileout << " generated from the model in Section \\ref{sec:model}.  We report the fraction of trials that produce a rejection of a nominal 5\\% test of fileout << " the model overidentifying restrictions. The first " << endl;
  fileout << " column is for a model that is specified exactly as in Section \\ref{sec:model}, and the second is for a model with a 0.001 fixed cost of fileout << " adjusting the capital stock.} " << endl;
  fileout << " \\end{center} " << endl;
  fileout << " \\end{table} " << endl;
}*/