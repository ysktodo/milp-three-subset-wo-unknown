This souce code is the verification code for the new-type cube attack introduced in CRYPTO 2018. 

This code is written by C++ with Gurobi API. 

Therefore, to compile and run this code, you need to install Gurobi Optimizer in advance. 

If you already install the Gurobi Optimizer version 8.1, you just type 
+++
 make 
+++

If your Gurobi Optimizer is not version 8.1, please change LIB option in makefile. 



After comple, please type
+++
 ./a.out -r 855 -t [core number]
+++

You can change the core number that Gurobi Optimizer uses. 


