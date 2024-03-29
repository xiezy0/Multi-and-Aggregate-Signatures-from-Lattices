PALISADE Lattice Cryptography Library Experimental Repository 
=====================================

PALISADE is a general lattice cryptography library that currently
includes efficient implementations of the following lattice
cryptography capabilities:
	
* Fully Homomorphic Encryption (FHE)
   * Brakerski/Fan-Vercauteren (BFV) scheme for integer arithmetic
   * Brakerski-Gentry-Vaikuntanathan (BGV) scheme for integer arithmetic
   * Cheon-Kim-Kim-Song (CKKS) scheme for real-number arithmetic
   * Ducas-Micciancio (FHEW) and Chillotti-Gama-Georgieva-Izabachene (TFHE) schemes for Boolean circuit evaluation
* Multi-Party Extensions of FHE (to support multi-key FHE)
   * Threshold FHE for BGV, BFV, and CKKS schemes
   * Proxy Re-Encryption for BGV, BFV, and CKKS schemes
   
This experimental repository holds the implementation of the following lattice-based cryptography capabilities (as of version 1.11, this implementations will no longer be included in the main PALISADE release).

* Digital Signature
	
Users should build and install PALISADE on their system. The PALISADE include and library files are not included in this repostiory.

PALISADE is a cross-platform C++11 library supporting Linux, Windows, and macOS. The supported compilers are g++ v6.1 or later and clang++ v6.0 or later.

PALISADE is available under the BSD 2-clause license.

Further information about PALISADE:

[License Information](License.md)

[Library Wiki with documentation](https://gitlab.com/palisade/palisade-development/wikis/home)

[Code of Conduct](Code-of-conduct.md)

[Governance](Governance.md)

[Contributing to PALISADE](Contributing.md)


Build Instructions
=====================================
This repository has been tested to run with PALISADE development release 1.11.

* Install PALISADE from that release on your system. Full instructions
  for this are to be found in the [README.md](https://gitlab.com/palisade/palisade-release/-/blob/master/README.md) file in the PALISADE
  repo.

Run `make install` at the end to install the system to the default
location (you can change this location, but then you will have to
change the Makefile in this repo to reflect the new location).

Note you may have to execute the following on your system to
automatically find the installed libraries and include files:

> `sudo ldconfig`

  found in that repo. 

* Clone this repo on your system 

We use CMake to build signatures. The high-level (platform-independent)
procedure for building PALISADE is as follows (for OS-specific
instructions, see the section "Detailed information about building
PALISADE" at the bottom of this page). Note PALISADE has similar
requirements, so if that builds on your system you are all set :

* Install system prerequisites (if not already installed), including a C++ compiler with OMP support, cmake, make, and autoconf.

* Create a directory where the binaries will be built. The typical choice is a subfolder "build". In this case, the commands are:
```
mkdir build
cd build
cmake ..
```

Note that CMake will check for any system dependencies that are needed
for the build process. If the CMake build does not complete
successfully, please review the error CMake shows at the end. If the
error does not go away (even though you installed the dependency), try
running "make clean" to clear the CMake cache.

* Build the executables by running the following command (this will take few minutes; using the -j make command-line flag is suggested to speed up the build)
```
make
```
If you want to build only library files or some other subset of signatures, please review the last paragraph of this page.

After the "make" completes, you should see the signature library file in the lib folder, binaries of examples in bin/examples and binaries for unit tests in the bin/unittest folder.

* Optionally, the library, `PALISADEsignature`,  can be installed.
```bash
make install
```
**Note** - The default installation path is `/usr/local/` and likely requires admin priviledges.

If this library is installed it must be found uniquely but with Palisade as well. See [CMakeLists.Users.txt](CMakeLists.User.txt) for how to use this library.

Testing and cleaning the build
-------------------

Run unit tests to make sure all capabilities operate as expected
```
make testsignature
```

Run sample code to test, e.g.,
```
bin/examples/signature/gpv
```

To remove the files built by make, you can execute
```
make clean
```

Supported Operating Systems
--------------------------
PALISADE CI continually tests our builds on the following operating systems:

* Ubuntu [18.04]
* macOS [Mojave]
* Centos 7
* NVIDIA Xavier [Linux for Tegra 4.2.2]
* MinGW (64-bit) on Windows 10

PALISADE users have reported successful operation on the following systems:

* FreeBSD
* Ubuntu [16.04]

Please let us know the results if you have run PALISADE any additional systems not listed above.

### Windows

PALISADE extensions are encapsulated in independent repositories and are therefore installed in unique paths. The following snippet gives an example of how to include the extensions when using MinGW:

```bash
palisadeLibs=(
    '/c/Program Files (x86)/PALISADE/lib' 
    '/c/Program Files (x86)/PALISADEabe/lib' 
    '/c/Program Files (x86)/PALISADEsignature/lib' 
    '/c/Program Files (x86)/PALISADEtrapdoor/lib')

for ((i=0; i < ${#palisadeLibs[*]}; i++))
do
    lib=${palisadeLibs[$i]}
    PATH="${PATH}:${lib}"
done
```
**Note** This can be pasted into `.bashrc`, `.profile`, etc. to make it perminent
