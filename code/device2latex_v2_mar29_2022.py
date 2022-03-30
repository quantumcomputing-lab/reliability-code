#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 20:41:22 2022
@author: samudra
Automated tex generation for input to latex2html
"""

import re
import os

output_file="/home/samudra/Desktop/sjv_website_dev/tex_subfile.txt"

device_list = {
'ibm_washington': 127,
'ibmq_brooklyn': 65,
'ibmq_toronto': 27,
'ibmq_mumbai': 27,
'ibm_cairo': 27,
'ibmq_montreal': 27,
'ibm_auckland': 27,
'ibm_hanoi': 27,
'ibmq_guadalupe': 16,
'ibm_perth': 7,
'ibm_lagos': 7,
'ibmq_jakarta': 7,
'ibmq_manila': 5,
'ibmq_bogota': 5,
'ibmq_santiago': 5,
'ibmq_quito': 5,
'ibmq_lima': 5,
}

text_before = r"""\subsection{washington}
\subsubsection{Metric 1: Initialization fidelity}
\begin{table}
\centering
%\caption{Caption text}
%\vspace{10mm}
\label{tbl:FI}
\begin{center}
\begin{tabular}{|c|c|}
\multicolumn{2}{c}{METRIC 1: INITIALIZATION FIDELITY}\\
\hline
\hline
Confidence interval (across time) & Confidence interval (point-in-time)\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_across-time.png}}
%\caption{A figure}
\label{fig:FI_across-time}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_point-in-time.png}}
%\caption{A figure}
\label{fig:FI_point-in-time}
\end{figure}
\\
\multicolumn{2}{c}{RAW DATA: POINT-IN-TIME}\\
\multicolumn{2}{c}{
\begin{figure}
\centering
\includegraphics[width=1.2\linewidth]{\htmladdimg{figures/washington/FI_raw_lattice.png}}
%\caption{A figure}
\label{fig:FI_raw_lattice}
\end{figure}
}
\\
\multicolumn{2}{c}{RAW DATA: ACROSS TIME}\\
\hline
\hline
Qubits # 0-15 & Qubits # 16-31\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_0.png}}
%\caption{A figure}
\label{fig:FI_raw_time_0}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_1.png}}
%\caption{A figure}
\label{fig:FI_raw_time_1}
\end{figure}
\\
\hline
\hline
Qubits # 32-47 & Qubits # 48-63\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_2.png}}
%\caption{A figure}
\label{fig:FI_raw_time_2}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_3.png}}
%\caption{A figure}
\label{fig:FI_raw_time_3}
\end{figure}
\\
\hline
\hline
Qubits # 64-79 & Qubits # 80-95\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_4.png}}
%\caption{A figure}
\label{fig:FI_raw_time_4}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_5.png}}
%\caption{A figure}
\label{fig:FI_raw_time_5}
\end{figure}
\\
\hline
\hline
Qubits # 96-111 & Qubits # 112-126\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_6.png}}
%\caption{A figure}
\label{fig:FI_raw_time_6}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FI_raw_time_7.png}}
%\caption{A figure}
\label{fig:FI_raw_time_7}
\end{figure}
\\
\hline
\end{tabular}
\end{center}
\end{table}

\subsubsection{Metric 2: Gate fidelity}
\begin{table}
\centering
%\caption{Caption text}
%\vspace{10mm}
\label{tbl:FG}
\begin{center}
\begin{tabular}{|c|c|}
\multicolumn{2}{c}{METRIC 2: GATE FIDELITY}\\
\hline
\hline
Confidence interval (across time) & Confidence interval (point-in-time)\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_across-time.png}}
%\caption{A figure}
\label{fig:FG_across-time}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_point-in-time.png}}
%\caption{A figure}
\label{fig:FG_point-in-time}
\end{figure}
\\
\multicolumn{2}{c}{RAW DATA: POINT-IN-TIME}\\
\multicolumn{2}{c}{
\begin{figure}
\centering
\includegraphics[width=1.2\linewidth]{\htmladdimg{figures/washington/FG_raw_lattice.png}}
%\caption{A figure}
\label{fig:FG_raw_lattice}
\end{figure}
}
\\
\multicolumn{2}{c}{RAW DATA: ACROSS TIME}\\
\hline
\hline
Set 1 & Set 2\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_0.png}}
%\caption{A figure}
\label{fig:FG_raw_time_0}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_1.png}}
%\caption{A figure}
\label{fig:FG_raw_time_1}
\end{figure}
\\
\hline
\hline
Set 3 & Set 4\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_2.png}}
%\caption{A figure}
\label{fig:FG_raw_time_2}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_3.png}}
%\caption{A figure}
\label{fig:FG_raw_time_3}
\end{figure}
\\
\hline
\hline
Set 5 & Set 6\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_4.png}}
%\caption{A figure}
\label{fig:FG_raw_time_4}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_5.png}}
%\caption{A figure}
\label{fig:FG_raw_time_5}
\end{figure}
\\
\hline
\hline
Set 7 & Set 8\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_6.png}}
%\caption{A figure}
\label{fig:FG_raw_time_6}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/FG_raw_time_7.png}}
%\caption{A figure}
\label{fig:FG_raw_time_7}
\end{figure}
\\
\hline
\end{tabular}
\end{center}
\end{table}

\subsubsection{Metric 3: Duty cycle}
\begin{table}
\centering
%\caption{Caption text}
%\vspace{10mm}
\label{tbl:tau}
\begin{center}
\begin{tabular}{|c|c|}
\multicolumn{2}{c}{METRIC 3: DUTY CYCLE}\\
\hline
\hline
Confidence interval (across time) & Confidence interval (point-in-time)\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_across-time.png}}
%\caption{A figure}
\label{fig:tau_across-time}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_point-in-time.png}}
%\caption{A figure}
\label{fig:tau_point-in-time}
\end{figure}
\\
\multicolumn{2}{c}{RAW DATA: ACROSS TIME}\\
\hline
\hline
Set 1 & Set 2\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_0.png}}
%\caption{A figure}
\label{fig:tau_raw_time_0}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_1.png}}
%\caption{A figure}
\label{fig:tau_raw_time_1}
\end{figure}
\\
\hline
\hline
Set 3 & Set 4\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_2.png}}
%\caption{A figure}
\label{fig:tau_raw_time_2}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_3.png}}
%\caption{A figure}
\label{fig:tau_raw_time_3}
\end{figure}
\\
\hline
\hline
Set 5 & Set 6\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_4.png}}
%\caption{A figure}
\label{fig:tau_raw_time_4}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_5.png}}
%\caption{A figure}
\label{fig:tau_raw_time_5}
\end{figure}
\\
\hline
\hline
Set 7 & Set 8\\
\hline
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_6.png}}
%\caption{A figure}
\label{fig:tau_raw_time_6}
\end{figure}
&
\begin{figure}
\centering
\includegraphics[width=0.8\linewidth]{\htmladdimg{figures/washington/tau_raw_time_7.png}}
%\caption{A figure}
\label{fig:tau_raw_time_7}
\end{figure}
\\
\hline
\end{tabular}
\end{center}
\end{table}

\subsubsection{Metric 4: Addressability}
\begin{table}
\centering
%\caption{Caption text}
%\vspace{10mm}
\label{tbl:nmi}
\begin{center}
\begin{tabular}{|c|c|}
\multicolumn{2}{c}{METRIC 4: 	ADDRESSABILITY}\\
\hline
\hline
\multicolumn{2}{c}{POINT-IN-TIME: NMI LATTICE}\\
\multicolumn{2}{c}{
\begin{figure}
\centering
\includegraphics[width=1.2\linewidth]{\htmladdimg{figures/washington/nmi.png}}
%\caption{A figure}
\label{fig:nmi}
\end{figure}
}
\\
\hline
\end{tabular}
\end{center}
\end{table}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


"""

def device2latex(text_before, before_devicename, after_devicename, output_file):

    text_before_rawed = r'{}'.format(text_before)
    pattern = re.compile(re.escape('washington'), re.IGNORECASE)
    text_after = pattern.sub(after_devicename, text_before_rawed )
    text_after_rawed = r'{}'.format(text_after)

    #output_file=after_devicename+'.txt'
    f = open(output_file, 'a')
    #f.write(text_after_rawed.encode('utf-16'))
    f.write(text_after_rawed)
    f.close()


try:
    os.remove(output_file)
except OSError:
    pass

for after_devicename in device_list:
    ind = after_devicename.find('_')
    device2latex(text_before=text_before, 
                 before_devicename = "washington", 
                 after_devicename = after_devicename[ind+1:],
                 output_file=output_file)
