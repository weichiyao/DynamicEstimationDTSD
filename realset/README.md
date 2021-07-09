### Description of the real dataset 

We revisit the bankruptcy data in Section 3 of Bou-Hamad et al. (2011a) from United States firms that conducted IPOs (Initial Public Offerings) between 1990 and 1999. The data come from COMPUSTAT and each firm is followed yearly starting from its IPO until 2006. The target variable is the bankruptcy, measured in number of years after the IPO. All firms that filed for bankruptcy under Chapter 7 or 11 are considered bankrupt. The data set includes 1143 firms, of which 189 went bankrupt during the study. The time-varying covariates are the following financial ratios collected yearly:
\begin{itemize}
  \item $X_1$ = Current Assets/Current Liabilities
  \item $X_2$ = Net Working Capital/Total Assets
  \item $X_3$ = Sales/Net Working Capital
  \item $X_4$ = Sales/Net Fixed Assets
  \item $X_5$ = Sales/Total Assets
  \item $X_6$ = Total Debt/Total Assets
  \item $X_7$ = Long--Term Debt/(Long--Term Debt+Total Equity)
  \item $X_8$ = Net Income/Total Assets
  \item $X_9$ = Net Income/Total Equity
  \item $X_{10}$ = Market Value of Equity/Book Value of Total Debt.\\
  In addition, the following time-invariant covariate is included:
  \item $X_{11}$ = Category based on Standard Industrial Classification (SIC) Division Structure.
\end{itemize}