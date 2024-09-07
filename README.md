# CnCovsimJe
Code to reproduce the results in the article submission to Proceedings of the Royal Society A journal.
It is neccesary a computer with linux or MacosX and install gcc and tclsh.

To reproduce the results of the Fig4, 

In first, compile the main files, as follow

gcc SIQR_COV_PAPER_RSTA_SEPT24_fig4.c -lm -o SIQR_COV_PAPER_RSTA_SEPT24_fig4

gcc SIQR_COV_PAPER_RSTA_SEPT24_fig5.c -lm -o SIQR_COV_PAPER_RSTA_SEPT24_fig5

Then, run the following simulations

nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 0.3 0 1 0 0.2213 5 5 10 1 100 &  

nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 0.7 0 1 0 0.2213 5 5 10 1 100 & 

nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 1 0 1 0 0.2213 5 5 10 1 100 &

nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig4 10000 1000 2 1 2 1 0 0.2213 5 5 10 1 100 &

nohup ./SIQR_COV_PAPER_RSTA_SEPT24_fig5 10000 1000 2 0.81 0 1 0 0.2213 5 10 15 1 100 &



