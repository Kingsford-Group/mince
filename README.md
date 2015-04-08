
Mince
=====

Mince is a technique for encoding collections of short reads. The Mince software is written in C++11 and is availble under an open-source license.

If you use this software, please cite:

* Carl Kingsford and Rob Patro. Data-dependent Bucketing Improves Reference-free Compression of Sequencing Reads. Under review (2014).

Installation
============

If you can, use the binaries available on github. These should work for MacOS
and most Linux variants. If you want to build from sources, see the relevant section below.

Usage
=====

To compress:

        mince -e -l U -r IN.fastq -o OUTBASENAME
    
or with a paired-end library

        mince -e -l IU -1 IN1.fastq -2 IN2.fastq -o OUTBASENAME
    
Mince will create several output files.

To decompress:

    	mince -d -i BASENAME -o OUT
    
where BASENAME is the basename for the file to decompress (same as the OUTBASENAME used when encoding). This will write the sequences in FASTA format to OUT.fa. If you originally encoded a paired-end set of reads, this will write the sequences in FASTA format to the files OUT1.fa and OUT2.fa.

The reads will NOT be in the same order as in the original file.

You can compress multiple files together using standard process subsitution syntax:

	mince [options] ­1 <(cat f1_1.fq f2_1.fq f3_1.fq) ­2 <(cat f1_2.fq f2_2.fq f3_2.fq)

Other Options
```
      -v [ --version ]           print version information
      -h [ --help ]              print help message
      -p [ --threads ] arg (=20) number of concurrent threads to use for
                                 compression / decompression; the minimum value is
                                 4.
      -l [ --libtype ] arg       library format string [only for encoding]
      -1 [ --mates1 ] arg        mate file 1 [only for encoding]
      -2 [ --mates2 ] arg        mate file 2 [only for encoding]
      -r [ --unmated_reads ] arg unmated reads [only for encoding]
      -i [ --input ] arg         input base file [only for decoding]
      -o [ --output ] arg        output base file [for both encoding / decoding]
      -b [ --blength ] arg (=15) length of the bucket string [1,16]
      -n [ --norc ]              don't allow reverse complementing when encoding /
                                 bucketing
      -e [ --encode ]            encode the input file
      -d [ --decode ]            decode the input file
```

The specification of the library format string is the same as is expected in
the sailfish and salmon software, and is described in detail
[here](http://sailfish.readthedocs.org/en/develop/salmon.html#what-s-this-libtype)
. A visual depiction of the different possible paired-end library types is
available [here](http://sailfish.readthedocs.org/en/develop/library_type.html).

Building from sources
=====================

### Dependencies

Building Mince requires:

1. A C++11 compatible compiler (it has been tested with GCC and clang)
2. The [CMake](http://www.cmake.org) build system 
3. Intel TBB library (The build system can download this if you don't have it)
4. PLZip program
5. Boost (The build system can download this if you don't have it)

### How to build

If you want to build Mince from source, follow these instructions. This is
fairly easy if you have the right dependencies.

First, Mince uses CMake to configure the build. To start, download the latest
source code (either using git or by downloading the source release from
github). 

Next, create a directory inside the mince directory:

```
    > mkdir build
    > cd build
```

Now, run cmake to configure Mince for your system:

```
    >  cmake [FLAGS] ..
```
The ".." should point to the root directory of the mince source tree that you
downloaded.  You should replace [FLAGS] as described below.

* -DFETCH\_BOOST=TRUE -- If you don't have Boost installed (or have an older
  version of it), you can provide the FETCH_BOOST flag instead of the
  BOOST_ROOT variable, which will cause CMake to fetch and build Boost locally.

* -DBOOST_ROOT= -- Tells CMake where an existing installation of Boost resides,
  and looks for the appropriate version in . This is the top-level directory
  where Boost is installed (e.g. /opt/local).

* -DTBB\_INSTALL_DIR= -- Tells CMake where an existing installation of Intel's
  TBB is installed (), and looks for the appropriate headers and libraries
  there. This is the top-level directory where TBB is installed (e.g.
  /opt/local).
  
* -DCMAKE\_INSTALL_PREFIX= -- is the directory to which you wish Sailfish to be installed. If you don't specify this option, it will be installed locally in the top-level directory (i.e. the directory directly above "build").

Mince may build if you omit any flags, but probably you will have to specify
either the -DFETCH\_BOOST (if you don't have boost installed), or -DBOOST\_ROOT
(if you have boost installed someplace non-standard). 

If the cmake command completes, then you can actually build the software via:

```
    > make
```

This will produce the executables in the "src" directory under your build directory. You can copy the binaries to your MAKE\_INSTALL\_PREFIX via:

```
    > make install
```

