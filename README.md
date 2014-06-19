sg-lite
=======

This is intended a light weight Sparse Grid implementation focused mostly on density estimation. This is a WIP and this readme will be filled with up to date comments about the process.

19/06/14: I have made the following design choices for now
		  - C++11 : For the time being the library will be written for C++11. The reason to use C++ is so that I 
		  	have access to both Eigen (a C++ template library for linear algebra) and std library containers. The
		  	reason to use C++11 is that this allows me to use shared_ptr without requiring boost. Eventually if it
		  	seems feasible I will try and copy in my own implementation of the required elements to try and move 
		  	closer to C.
		  - RegularGrid: I have decided that I will first focus on getting algoirthms working for the regular grid,
		    then for the combination grid, then finally for a full SparseGrid. This may make a regular Sparse Grid
		    difficult to implement if it becomes needed, but there are other libraries for this, and starting from
		    that point seems top heavy when I hope for this to be a light weight implementation
		  - Avoiding Recursion: I have decided to try and avoid recursion in algorithms where possible. I do not 
		    have a good reason for this, other than that it might be faster, and seems like a fun thing to try...


Comments: I appreciate any comments on the current state of the code, but do keep in mind it is a WIP. I will try
          and keep minimial build instructions in this readme as it becomes realistic, so that people can see what 
          should be working.

-------
INSTALL
-------

A minimal example should be workable by simply calling 'make' in bash. This will create a 'test' executable. 
Further details to foolow. (19/06/14)