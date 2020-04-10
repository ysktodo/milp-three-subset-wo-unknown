This source code recovers superpolys for Trivium and Grain-128AEAD. 

This code is written by C++ with Gurobi API. 

Therefore, to compile and run this code, you need to install Gurobi Optimizer in advance. 

If you already install the Gurobi Optimizer version 8.1, you just run
+++
 make
+++

If your Gurobi Optimizer is not version 8.1, please change LIB option in makefile. 

After compile, if you want to try the superpoly recovery for 840- or 841-round Trivium, you just run
+++
	./a.out -r [840 or 841] -trivium -t [option : thread number]
+++

Note that this code does not return the answer quickly. 
It depends on the performance of your computer, and if you execute this code in a cheap computer, you need to wait a few days. 
We highly recommend that this code is executed on the computer with good performance. 

If you want to try the superpoly recovery for 190-round Grain-128AEAD, you just run
+++
	./a.out -r 190 -grain -t [option : thread number]
+++

Moreover, you want to try 15 superpolies that are used in the key-recovery attack against Grain-128AEAD, you just run
+++
	./a.out -r 190 -grain -subcube -t [option : thread number]
+++

Similarly to the case of Trivium, this code does not return the answer quickly. 
Therefore, we highly recommend that this code is executed on the computer with good performance. 


This source code also provides the practical verification, where the superpoly is recovered under the randomly chosen cube whose size is at most several bits and the correctness of the recovered superpoly is experimentally verified by using 100 randomly generated secret key bits and non-cube IV bits. 
If you want to try this verification, you just run
+++
	\tt{./a.out -trivium -practical}
+++
for Trivium and 
+++
	\tt{./a.out -grain -practical}
+++
for Grain-128AEAD. 

