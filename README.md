wsp3dovm
========

Weighted Shortest Path with OpenVolumeMesh (wsp3dovm).

This is a Visual Studio 2019 solution but it should work in newer versions of Visual Studio as well.

wsp3dovm depends on two external libraries:

* Boost (tested with version 1.70)
* OpenVolumeMesh (tested with version 2.0)

These libraries can be installed in separate folders which are parallel to wsp3dovm. 
This is the easiest way because the wsp3dovm.sln file is already configured for this.
Depending on your machine you might install these libraries into other folders and then 
edit the includes and library paths in wsp3dovm accordingly.

In the following, the common root folder for all code is C:\Users\Frank\Source


Download and install Boost
--------------------------
Download the prebuilt Windows binaries (installer) boost_1_70_0-unsupported-msvc-14.2-64.exe from http://www.boost.org/
(1´4.2 indicates the Visual Studio tools version) and install it to a new empty folder, C:\local\boost_1_70_0.


Download, build and install OpenVolumeMesh
------------------------------------------
Download a zipped snapshot of the source tree from http://www.openvolumemesh.org/ under Downloads
and extract the source tree to C:\Users\Frank\Source\OpenVolumeMesh.

Open and read the included README which contains instructions to build OpenVolumeMesh by using CMake.

Create the build folder C:\Users\Frank\Source\OpenVolumeMesh\build for CMake output.
The CMake generator should be set to "Visual Studio 12 2013 Win64".
Now, choose Configure, set the Configuration variable CMAKE_INSTALL_PREFIX 
to the folder C:\Users\Frank\Source\OpenVolumeMesh itself
and the generate a Visual Studio 12 solution from CMake.

CMake may complain some missing features, but this can be ignored.

CMake creates a Visual Studio Solution OpenVolumeMesh.sln, open it.
The build configuration should be x64 Debug.
In project OpenVolumeMesh navigate to header file TopologyKernel.hh and open it.
Go to line 580, it should read "private:", change that to "public:" and save it.

Now build ALL_BUILD, then INSTALL. This will create two new subfolders include and lib in the OpenVolumeMesh folder.
(Do not build OpenVolumeMesh directly, this will fail.)


Build wsp3dovm
--------------
Open Visual Studio and use the Team Explorer to clone the git repo from github 
to a local repository in folder C:\Users\Frank\Source\wsp3dovm.

See blogs.msdn.com/b/visualstudioalm/archive/2013/02/06/set-up-connect-and-publish-using-visual-studio-with-git.aspx

Open wsp3dovm.sln. It might complain that project Zurich1 is not found, this can be ignored.
The main C++ project is wsp3dovm, set it as the startup project. (Right click on it in the Solution Explorer window.)

Set the build configuration to x64 Debug and build the solution. 
The build has some errors related to private class members like outgoing_hes_per_vertex_.
Go to such an error location and click on Go To Definition. This will open TopologyKernel.hh
from OpenVolumeMesh. Change private to public and rebuild. This time there should be no more errors but some warnings.


Dry Run
-------
wsp3dovm needs some Boost DLLs for execution. In Visual Studio edit the project properties. 
Under Configuration Properties > Debugging edit Environment. Set the following PATH:

	PATH=$(SolutionDir)\..\..\boost_1_55_0\lib64-msvc-12.0
	
which is automatically added before wsp3dovm is executed from within Visual Studio. 
You may also add the folder C:\Users\Frank\Source\Repos\boost_1_55_0\lib64-msvc-12.0
to the computers PATH variable to run wsp3dovm.exe without Visual Studio.

In Visual Studio press Ctrl+F5 to execute wsp3dovm. 
Because there are no commandline parameters given, a usage message will be displayed and the program exits immediately.


First Test (unweighted case)
----------------------------

Input: A 3D tetrahydralization consisting of two text files:

* .node file containing 3D coordinates of all vertices (3D points) of the tetrahydralization, one per line.
* .ele file containing all tetrahedra (4 vertex indices), one per line.

The first lines contain header information, for details see http://wias-berlin.de/software/tetgen/1.5/doc/manual/manual006.html.

Example test1:

test1.node

	5	3 0 0
	0	1.0	2.0	 1.0
	1	4.0	0.0	 0.0
	2	0.0	4.0	 0.0
	3	0.0	0.0	 0.0
	4	2.0	1.0	-1.0

test1.ele

	2	4 0
	0	0 1 2 3
	1	1 2 3 4

Run wsp3dovm with the commandline arguments

	--start_vertex 0 --termination_vertex 4 --yardstick 0.3 --write_mesh_vtk 1 test1

This will output some statistics and generate the following output files:

* test1.vtk
* test1._wsp_path_s0_t4.vtk
* test1._wsp_tree_s0.vtk
* test1._wsp_path_cells_s0_t4.*

The .vtk files can be visualized by ParaView from kitware http://www.kitware.com/opensource/opensource.html

The file test1.vtk contains the input tetrahedralization in .vtk format for visualization
The file test1._wsp_path_s0_t4.vtk contains the approximated shortest path between start vertex 0 and termination vertex 4
The file test1._wsp_tree_s0.vtk contains the tree of all single source shortest paths starting at vertex 0
The files test1._wsp_path_cells_s0_t4 represent the "3D buffer" around the approximated shortest path, a subset of the input tetrahedralization
which can be used for refinements.


Second Test (weighted case)
---------------------------

Weights can be added easily to the .ele file. The .node file remains the same. 
In this example, we use the refractive index of air and water:

test2.node

	5	3 0 0
	0	1.0	2.0	 1.0
	1	4.0	0.0	 0.0
	2	0.0	4.0	 0.0
	3	0.0	0.0	 0.0
	4	2.0	1.0	-1.0

test2.ele

	2	4 1
	0	0 1 2 3   	1.000293
	1	1 2 3 4		1.3330

Rerun wsp3dovm and compare the approximated shortest paths.


Switching to Release Mode
-------------------------
When the firsts tests run successfully, it is a good idea to switch from Debug mode to Release mode for
wsp3dovm and OpenVolumeMesh to increase execution speed. Rerun the tests to make sure that everything works okay.

Since the resulting data structures become huge, 64-bit mode 8x64) is highly recommended.


Tetgen
------
The tetgen software by Hang Si (http://wias-berlin.de/software/tetgen/) can generate 3D tetrahydralizations 
from several input formats and can be used as a preprocessing step. It comes with a makefile but can
easily be compiled using visual Studio into an executable file tetgen.exe.


A Larger Example (twolayermdl)
------------------------------
Here are some examples of 3D structures:

http://wias-berlin.de/software/tetgen/1.5/fformats.examples.html

We use twolayermdl which is a triangular surface mesh of some geographical sub-surface structure and consists of two files:

* twolayermdl.node
* twolayermdl.smesh

Download these files and use tetgen.exe to create a tetrahedralization:

	tetgen.exe -p -q1.4/18 -V twolayermdl.smesh

The commandline parameters  
* specify the input as an boundary description (surface mesh) of a 3D piecewise linear complex (-p)
* request a maximum radius/edge ratio bound of 1.4 and a minimum dihedral angle bound of 18° (-q1.4/18) and
* print a mesh quality report (-V)

Consult the tetgen manual by Hang Si for details: http://wias-berlin.de/software/tetgen/1.5/doc/manual/index.html

Tetgen creates 4 output files:

* twolayermdl.1.edge
* twolayermdl.1.ele
* twolayermdl.1.face
* twolayermdl.1.node

We use the .ele and .node file as in the above tests. Now we execute a batch script (run_twolayermdl.cmd) of shortest path approximations in a Windows command shell (cmd).

run_twolayermdl.cmd:

	set W=x64\Release\wsp3dovm.exe

	for /L %%Y in (200,-10,100) do (
	  %W% --random_s_t_vertices 100 --spanner_stretch 0.0 --yardstick %%Y twolayermdl.1 >> logfile
	)

The yardstick defines the length of subdivision of each tetrahedron, i.e. Steiner points are created  with a max. distance of yardstick.

Decreasing values between 200 and 100 and a step size of -10 are used in this script.

Since script execution will take several minutes, I often use Process Explorer from https://technet.microsoft.com/en-us/sysinternals
for supervision of CPU and memory usage while those scripts are running. It shows memory usage of up to 6.7 GB and 100% use of 1 CPU core.

Now extract the average shortest path approximation ratio for each value of yardstick. 
(Here it is assumed that the 3D structure is convex and hence the Euclidean distance between the end points is the true length of the shortest path.)

	find "avg shortest path approximation ratio" logfile

The output should look like

	avg shortest path approximation ratio: 1.03114
	avg shortest path approximation ratio: 1.02966
	avg shortest path approximation ratio: 1.02717
	avg shortest path approximation ratio: 1.02563
	avg shortest path approximation ratio: 1.02419
	avg shortest path approximation ratio: 1.02197
	avg shortest path approximation ratio: 1.01988
	avg shortest path approximation ratio: 1.01832
	avg shortest path approximation ratio: 1.01604
	avg shortest path approximation ratio: 1.01437
	avg shortest path approximation ratio: 1.01261

showing how the approximation improves when the yardstick decreases.

It could be favorable to either switch to Linux or use cygwin under Windows. This would allow bash scripts and more powerful tools like grep, awk, perl etc.. to extract statistics data.
