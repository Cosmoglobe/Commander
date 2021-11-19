#-------------------------------------------------------------
# Bourne/Korn shell functions required to configure Healpix packages
# ------------------------------------------------------------
# 2008-03-26, IAP, EH
# 2008-07-24  ---> version 2.10
# 2008-09-22, EH: remove need for IDL_DIR
#             no long term effects (on IDL_PATH and IDL_STARTUP) of hidl or hidlde
#             do not expand HEALPIX variable in F90 config file
# 2008-09-26: more stringent tests on F90 compiler (checkF90Compilation) and
#                linkage to fitsio (checkF90FitsioLink)
#             clean-up macros creating 'to_be_removed' files
# 2008-11-13  ---> version 2.11
# 2008-11-17: corrected typo in C 'c-shared' target
#             detects correctly gfortran for version >= 4.3
# 2008-11-21: solved potential problem with multiple cfitio*.tar.gz for C++
#             introduced ${HEAD}
#             replaced ~ with ${HOME}
# 2009-02-25: added Fortran compilation flag for 64 bit INTEGER
#            still unknown for Fujitsu, Lahey, Portland, Absoft
# 2009-06-18: added WLRPATH so that F90 codes can be linked to shared cfitsio library
# 2009-06-26: replace echoLn with printf
# 2009-07-10: debugging on MacOS, libgif -> libhpxgif
# 2009-10-12: replace gfortran -dumpversion -> gfortran --version
# 2009-10-21: removed bashisms:
#               replaced ' == ' tests with ' = ' or ' -eq '
#               got rid of arrays in pickCppCompilation
# 2010-06-22: supports zsh (M. Tomasi)
# 2010-12-09: added IdentifyCParallCompiler to compile C libpsht with OpenMP
# 2011-01-28: C++ configuration ask for preinstalled cfitsio library
# 2011-01-31: keep track of previous choice of FITSDIR and FITSINC (within same session)
#           : propose OpenMP by default
# 2011-03-07: allow linking with shared libcfitsio for the C++ port
# 2012-02-27: better parsing of config.* files in C++ configuration
# 2012-05-30:    and ignore healpy specific config.* files.
# 2012-11-05: supports python (healpy) configuration
#             proposes -fPIC compilation of F90 code
# 2013-04-18: work-around for GCC 4.4 bug
# 2013-07-26: F90: add output location of modules ($MODDIR). Hacked from CMake.
# 2014-11-25: propose cfitsio-free compilation of C package
# 2015-05-12: correct bashism (==) introduced above (problematic for dash and zsh)
# 2015-07-31: improved g95 support; updated support address
# 2016-04-28: tentatively added MINGW for Windows
# 2016-06-02: debugged gcc detection in IdentifyCCompiler
# 2016-08-10: 1st attempt to better detect python version
# 2017-06-21: deal with FL (Fawlty Language) IDL clone (generateConfFlFile)
# 2019-01-21: call autotool configure for C++;
#             read system variables
#                CC,            (used in C/sharp/C++/F90/healpy) 
#                C_FITS, C_SHARED       (C)                      
#                CXX,                           (C++/    healpy) 
#                CXXFLAGS                       (C++)            
#                SHARP_COPT               (sharp)                
#                FC, F_DIRSUFF, F_OPT, F_PARAL, F_SHARED   (F90) 
#                FITSDIR,               (C/      C++/F90)        
#                FITSINC,               (C/      C++)            
#                PYTHON                                 (healpy) 
#                PROFILE_EDIT           (profile)                
#                papersize, ps_com, pdf_com, gif_com       (IDL) 
#             to define proposed default values.
# 2019-10-01: addition and documentation of --auto=list   mode
# 2020-01-06: added missing space in whereisCmd
# 2020-01-14: removed spurious \r in a some comments
# 2020-01-23: removed bashism, improved FITSDIR behavior
# 2021-01-26: better handling of libcfitsio without libcurl
#=====================================
#=========== General usage ===========
#=====================================
#   checkDir: search for installation directories and create them
#             is necessary
#   echoLn:
#   findFITSLib: search for FITSIO library
#   fullPath: convert relative to absolute directory names
#
#
#-------------
checkDir () {
    l=""
    for d in $*; do
	[ ! -d $d ] && l="$l $d"
    done
    if [ "x$l" != "x" ]; then
	echo "Warning: The following directories could not be found:"
	for d in $l; do
	    echo "$d"
	done
	if [ $INTERACTIVE -eq 1 ] ; then 
	    echoLn "Should I attempt to create these directories (y|n) [y]? "
	    read answer
	    if [ "x$answer" != "xn"  -a  "x$answer" != "xN"  ]; then
		for d in $l; do
		    ${MKDIR} $d 1>${DEVNULL} 2>&1
		    if [ $? -gt 0 ]; then
			echo "Error: Could not create directory $d"
			crashAndBurn
		    fi
		done
	    else
		echo "Create installation directories first."
		crashAndBurn
	    fi
	else 
	    # remove extra question in automated mode, to make script more predictable
	    echoLn "They now will be created ... "
	    for d in $l; do
		${MKDIR} $d 1>${DEVNULL} 2>&1
		if [ $? -gt 0 ]; then
		    echo "Error: Could not create directory $d"
		    crashAndBurn
		fi
	    done
	    echo "done."
	fi

    fi
}
#-------------
echoLn () {
#     if [ "${OS}" = "Linux" -o "${OS}" = "Darwin" ]; then
# 	echo -n "$*"
#     else
# 	echo "$*\c"
#     fi
    ${PRINTF} "$*"
}
#-------------
findFITSLib () {
    for dir in $* /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 \
    /usr/local/lib/cfitsio /usr/local/lib64/cftisio /usr/local/src/cfitsio \
    ${HOME}/lib ${HOME}/lib64 \
    `${LS} -dr /softs/cfitsio*/lib 2> ${DEVNULL}` \
    `${LS} -dr /softs/cfitsio/*/lib 2> ${DEVNULL}` \
    `${LS} -dr /usr/common/usg/cfitsio/*/lib 2> ${DEVNULL}` \
    `${LS} -dr ${HOME}/?oft*/cfitsio*/ 2> ${DEVNULL}` \
    `${LS} -dr ${HOME}/?oft*/cfitsio*/lib 2> ${DEVNULL}` \
    `${LS} -dr ${HEADAS}/../heacore/*-*-*/lib 2> ${DEVNULL}` ; do
	if [ -r "${dir}/lib${LIBFITS}.a" -o -r "${dir}/lib${LIBFITS}.so" -o -r "${dir}/lib${LIBFITS}.dylib" ] ; then
	    FITSDIR=$dir
	    break
	fi
    done
}
#-------------
findFITSInclude () {
    for dir in $* /usr/include /usr/local/include \
    /usr/local/src/cfitsio \
    ${HOME}/include ${HOME}/include64 \
    `${LS} -dr /softs/cfitsio*/include 2> ${DEVNULL}` \
    `${LS} -dr /softs/cfitsio/*/include 2> ${DEVNULL}` \
    `${LS} -dr /usr/common/usg/cfitsio/*/include 2> ${DEVNULL}` \
    `${LS} -dr ${HOME}/?oft*/cfitsio*/ 2> ${DEVNULL}` \
    `${LS} -dr ${HOME}/?oft*/cfitsio*/include 2> ${DEVNULL}` \
    `${LS} -dr ${HEADAS}/../heacore/*-*-*/include 2> ${DEVNULL}` ; do
	if [ -r "${dir}/fitsio.h" ] ; then
	    FITSINC=$dir
	    break
	fi
    done
}
#-------------
findFITSPrefix () {
    for dir in $* /usr /usr/local /usr/local/lib/cfitsio /usr/local/cfitsio /usr/local/lib64/cftisio /usr/local/src/cfitsio ${HOME}/softs/cfitsio/3.24 /usr/common/usg/cfitsio/3.26 ; do
	testlib="${dir}/lib/lib${LIBFITS}"
	if ( ([ -r "${testlib}.a" ] || [ -r "${testlib}.so" ] || [ -r "${testlib}.dylib" ]) && [ -r "${dir}/include/fitsio.h" ] ) ; then
	    FITSPREFIX=$dir
	    break
	fi
    done
}
#-------------
fullPath () {
    t='TEMP=`cd $TEMP; pwd`'
    for d in $*; do
	eval `echo $t | sed 's/TEMP/'$d'/g'`
    done
}
#-------------
askMake () {
if [ "${MAKESET}" = 0 ] ; then
    echoLn "Enter make command ($MAKE): "
    read answer
    [ "x$answer" != "x" ] && MAKE=$answer
    MAKESET=1
fi
}
#-------------
goodBye () {
    echo
    if [ -s Makefile -a ${edited_makefile} -eq 1 ] ; then
	echo
	echo "You can run \"(GNU)make\" to build all the packages configured so far,"
	echo "        and \"(GNU)make test\" to test them."
	echo " If supported, \"make -j\" will parallelize the tasks"
	echo "and speed up the compilation and installation"
	echo
    fi
    echo "Good Bye !"
    echo
    exit 0
}
#-------------
crashAndBurn () {
    echo
    echo "Something went wrong ..."
    echo "Quitting configuration script !"
    echo
    exit -1
}
#-------------
whereisCmd () {
    for d in $*; do
	CMD=`${WHEREIS} $d`
	if [ "x${CMD}" != "x"  -a "x${CMD}" != "x${d}:" ] ; then
	    break
	fi
    done
}
#-------------
isTrue () {
    [ "$1" = "1"  -o "$1" = "y"  -o "$1" = "Y" -o "$1" = "t"  -o "$1" = "T" ] && echo 1 || echo 0
}
#-------------
isFalse () {
    [ "$1" = "0"  -o "$1" = "n"  -o "$1" = "N" -o "$1" = "f"  -o "$1" = "F" ] && echo 1 || echo 0
}
#=================
# automatic mode
#=================
# fillFile () { #$1=leading line, $2=number of blank lines, $3=output file
#     lead=$1
#     rep=$2
#     multicr='XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
#     string=$(echo $multicr | ${AWK} -v var=$rep '{ string=substr($0, 1, var); print string; }' )
#     printf `echo ${lead}X${string} | sed 's|X|\\\n|g'` >> $3
# }
fillFile () { #$1=leading line, $2=number of blank lines, $3=output file
    echo $1 $2 | ${AWK} '{print $1} {while ($2--) print ""}' >> $3
}

processAutoList(){

    autofile=$1
    autolist=$2

    do_profile=0
    do_c=0
    do_cxx=0
    do_f90=0
    do_idl=0
    do_healpy=0
    do_sharp=0
    sharp_configured=0
    test_Sharp

    # ---- parse list ----
    #autolist=${autolist/,/ } # not in dash
    #autolist=`${SED} 's|,| |g' <<< $autolist` # not in dash
    autolist=`echo $autolist | ${SED} 's/,/ /g'`
    for item in $autolist; do
	case $item in
	    all)
            do_profile=1
            do_c=1
            do_sharp=1
            do_cxx=1
            do_f90=1
            do_idl=1
            do_healpy=1
            ;;
	    profile)
            do_profile=1
            ;;
	    c)
	    do_profile=1
	    do_c=1
            ;;
	    cxx|cpp)
	    do_profile=1
	    [ "$sharp_configured" = "0" ] && do_sharp=1 # only if not already configured
	    do_cxx=1
            ;;
	    f90)
	    do_profile=1
	    [ "$sharp_configured" = "0" ] && do_sharp=1 # only if not already configured
	    do_f90=1
            ;;
	    idl)
	    do_profile=1
	    do_idl=1
            ;;
	    healpy)
	    do_healpy=1
            ;;
	    sharp)
	    do_sharp=1
            ;;
	    *)
	    echo "Unknown item listed in automatic mode: $item"
	    echo "must be among: ${usage4}"
	    crashAndBurn
	esac
    done

    # ---- fill file ----
    ${RM} -f ${autofile}
    touch    ${autofile}
    #
    [ $do_profile -eq 1 ] && fillFile 9 1 ${autofile}
    [ $do_sharp   -eq 1 ] && fillFile 7 2 ${autofile}
    if [ $do_c    -eq 1 ]; then
	if [ `isTrue $C_FITS` -eq 1 ]; then
	    fillFile 2 8 ${autofile}
	else
	    fillFile 2 5 ${autofile}
	fi
    fi
    [ $do_cxx -eq 1 ]    && fillFile 4 6  ${autofile}
    [ $do_f90 -eq 1 ]    && fillFile 3 13 ${autofile}
    [ $do_idl -eq 1 ]    && fillFile 1 4  ${autofile}
    [ $do_healpy -eq 1 ] && fillFile 5 3  ${autofile}
    fillFile 0 0 ${autofile} # to exit gracefully

    do_profile=0
    do_c=0
    do_cxx=0
    do_f90=0
    do_idl=0
    do_healpy=0
    do_sharp=0

}
#=====================================
#=========== C package ===========
#=====================================
# setCDefaults
# add64bitCFlags
# askCUserMisc
# editCMakefile
# writeCpkgconfigFile
# C_config

setCDefaults () {

    OPT="-O2 -Wall"
    C_AR="ar -rsv"
    PIC="-fPIC" # works with gcc, icc and clang
    C_WLRPATH=""

    case $OS in
        AIX)
	    CC="xlc"
            OPT="-O -DRS6000"
	    PIC="-G"
            CF64="-q64"
            AR64="-X64"
	;;
	Linux)
	    C_WLRPATH="-Wl,-R"
	;;
    esac

    LIBDIR=$HEALPIX/lib
    INCDIR=$HEALPIX/include
    CDIR="src/C"

    FITSINC="/usr/local/include"

    #DO_C_SHARED=0
}
# -------------
add64bitCFlags () {

    if [ "x$CF64$AR64" != "x" ]; then
	echo "Do you want to make a 64 bit compilation ? (y|n) [n]: "
	read answer
	if [ "x$answer" = "xy" -o "x$answer" = "xY" ]; then
	    OPT="$OPT $CF64"
	    C_AR="$C_AR $AR64"
	fi
    fi
}
#-------------
askCUserMisc () {

    checkDir $INCDIR $LIBDIR
    fullPath INCDIR LIBDIR

    echoLn "enter C compiler you want to use [$CC]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CC=$answer

    add64bitCFlags

    echoLn "enter options for C compiler [$OPT]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && OPT=$answer

    echoLn "enter archive creation (and indexing) command [$C_AR]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && C_AR=$answer

    echoLn "do you want the HEALPix/C library to include CFITSIO-related functions ? (y|n)"
    if [ `isTrue ${C_FITS}` -eq 1 ] ; then
	echoLn " [y]: "
	read answer
	[ `isFalse ${answer}` -eq 1 ] && C_WITHOUT_CFITSIO=1 || C_WITHOUT_CFITSIO=0
    else
	echoLn " [n]: "
	read answer
	[ `isTrue ${answer}` -eq 1 ] && C_WITHOUT_CFITSIO=0 || C_WITHOUT_CFITSIO=1
    fi
    [ $INTERACTIVE -eq 0 ] && echo $answer
    if [ ${C_WITHOUT_CFITSIO} -eq 0 ]; then
	echoLn "enter full name of cfitsio library [lib${LIBFITS}.a]: "
	read answer
	[ $INTERACTIVE -eq 0 ] && echo $answer
	[ "x$answer" != "x" ] && LIBFITS=`${BASENAME} $answer ".a" | ${SED} "s/^lib//"`

	findFITSLib $LIBDIR $FITSDIR
	echoLn "enter location of cfitsio library [$FITSDIR]: "
	read answer
	[ $INTERACTIVE -eq 0 ] && echo $answer
	[ "x$answer" != "x" ] && FITSDIR=$answer
	fullPath FITSDIR

	[ "x$C_WLRPATH" != "x" ] && C_WLRPATH="${C_WLRPATH}${FITSDIR}"

	lib="${FITSDIR}/lib${LIBFITS}.a"
	if [ ! -r $lib ]; then
	    echo "error: fits library $lib not found"
	    crashAndBurn
	fi
	checkFitsioCurl $lib

	guess1=${FITSDIR}
	guess2=`${DIRNAME} ${guess1}`
	guess3="${guess2}/include"

	findFITSInclude $INCDIR ${guess1} ${guess2} ${guess3} $FITSINC
	echoLn "enter location of cfitsio header fitsio.h [$FITSINC]: "
	read answer
	[ $INTERACTIVE -eq 0 ] && echo $answer
	[ "x$answer" != "x" ] && FITSINC=$answer
	fullPath FITSINC

	inc="${FITSINC}/fitsio.h"
	if [ ! -r $inc ]; then
	    echo "error: cfitsio include file $inc not found"
	    crashAndBurn
	fi
    else
	C_WITHOUT_CFITSIO=1
    fi

    echoLn "A static library is produced by default. Do you also want a shared library ? "
    if [ `isTrue $C_SHARED` -eq 1 ]; then
	echoLn "(y|n) [y]: "
	read answer
	[ `isFalse $answer` -eq 1 ] && DO_C_SHARED=0 || DO_C_SHARED=1
    else
	echoLn "(y|n) [n]: "
	read answer
	[ `isTrue $answer` -eq 1 ] && DO_C_SHARED=1 || DO_C_SHARED=0
    fi
    [ $INTERACTIVE -eq 0 ] && echo $answer
}
#-------------
editCMakefile () {

    echoLn "Editing top Makefile  for C ..."

    clibtypes='c-static'
    if [ ${DO_C_SHARED} -eq 1 ]; then
	if [ "${OS}" = "Darwin" ]; then
	    clibtypes="${clibtypes} c-dynamic"
	else
	    ###clibtypes="${clibtypes} shared" # corrected 2008-11-17
	    clibtypes="${clibtypes} c-shared"
	fi
    fi

#    [ -r Makefile ] && ${CP} Makefile Makefile.bkup

    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^C_CC.*$|C_CC        = $CC|"   |\
	${SED} "s|^C_PIC.*$|C_PIC       = $PIC|"   |\
	${SED} "s|^C_OPT.*$|C_OPT       = $OPT|"   |\
	${SED} "s|^C_INCDIR.*$|C_INCDIR      = $INCDIR|"   |\
	${SED} "s|^C_LIBDIR.*$|C_LIBDIR      = $LIBDIR|"   |\
	${SED} "s|^C_AR.*$|C_AR        = $C_AR|"   |\
	${SED} "s|^C_WITHOUT_CFITSIO.*$|C_WITHOUT_CFITSIO = $C_WITHOUT_CFITSIO|" |\
	${SED} "s|^C_CFITSIO_INCDIR.*$|C_CFITSIO_INCDIR = $FITSINC|" |\
	${SED} "s|^C_CFITSIO_LIBDIR.*$|C_CFITSIO_LIBDIR = $FITSDIR|" |\
	${SED} "s|^C_EXTRA_LIB.*$|C_EXTRA_LIB = $CFITSIOCURL|" |\
	${SED} "s|^C_WLRPATH.*$|C_WLRPATH = $C_WLRPATH|" |\
	${SED} "s|^C_ALL.*|C_ALL     = ${clibtypes}|" |\
	${SED} "s|^HEALPIX=.*$|HEALPIX	= $HEALPIX|" |\
	${SED} "s|^ALL\(.*\) c-void \(.*\)|ALL\1 c-all \2|" |\
	${SED} "s|^TESTS\(.*\) c-void \(.*\)|TESTS\1 c-test \2|" |\
	${SED} "s|^CLEAN\(.*\) c-void \(.*\)|CLEAN\1 c-clean \2|" |\
	${SED} "s|^DISTCLEAN\(.*\) c-void \(.*\)|DISTCLEAN\1 c-distclean \2|" |\
	${SED} "s|^TIDY\(.*\) c-void \(.*\)|TIDY\1 c-tidy \2|" > Makefile

    echo " done."
    edited_makefile=1

}
# -----------------------------------------------------------------
writeCpkgconfigFile (){
    #pkgconfigFile=${HEALPIX}/lib/chealpix.pc
    pkgconfigFile=${HEALPIX}/lib/pkgconfig/chealpix.pc
    ${MKDIR} -p `${DIRNAME} ${pkgconfigFile}`
    echo
    echo "Writing pkgconfig file: ${pkgconfigFile}"

    echo "# HEALPix/C pkg-config file"         > ${pkgconfigFile}
    echo " "                                   >> ${pkgconfigFile}
    echo "prefix=${HEALPIX}"                   >> ${pkgconfigFile}
    echo "libdir=\${prefix}/lib"               >> ${pkgconfigFile}
    echo "includedir=\${prefix}/include"       >> ${pkgconfigFile}
    echo " "                                   >> ${pkgconfigFile}
    echo "Name: chealpix"                      >> ${pkgconfigFile}
    echo "Description: C library for HEALPix (Hierarchical Equal-Area iso-Latitude) pixelisation of the sphere" >> ${pkgconfigFile}
    echo "Version: ${HPXVERSION}"             >> ${pkgconfigFile}
    echo "URL: https://healpix.sourceforge.io" >> ${pkgconfigFile}
    if [ ${C_WITHOUT_CFITSIO} -eq 0 ] ; then
	echo "Requires: cfitsio"               >> ${pkgconfigFile}
    fi
    echo "Libs: -L\${libdir} -lchealpix"       >> ${pkgconfigFile}
    echo "Cflags: -I\${includedir} ${PIC}"     >> ${pkgconfigFile}

    ${CAT} ${pkgconfigFile}
    echo "           -------------- "
    echo

}

#-------------
C_config () {

    setCDefaults
    askCUserMisc
#    makeCInstall
    editCMakefile
    [ $NOPROFILEYET = 1 ] && installProfile
    writeCpkgconfigFile

}


#=====================================
#========= libsharp package ==========
#=====================================
#   setSharpDefaults: default variables for libsharp
#-------------
setSharpDefaults () {

    SHARPDIR=${HEALPIX}/src/common_libraries/libsharp
    SHARPBLD=${SHARPDIR}/build
    npgc=`$CC -V 2>&1          | ${GREP} PGI | ${WC} -l` # PGI C
    if [ $npgc != 0 ] ; then
	SHARP_COPT="${SHARP_COPT--O3 -fast -mp}" # -O3 -fast -mp  unless already defined
    else
	#SHARP_COPT="-DMULTIARCH -O3 -ffast-math -fopenmp" # reserve for multiplatform compilation
	SHARP_COPT="${SHARP_COPT--O3 -ffast-math -march=native -fopenmp}" # -O3 -ffast-math -march=native -fopenmp unless already defined
    fi
    SHARP_LDFLAGS=""
    if [ "${OS}" = "Darwin" -a `isTrue ${USE_ATRPATH}` -eq 1 ]; then
	#SHARP_LDFLAGS="-Wl,-install_name,@rpath/libsharp.0.dylib"
	SHARP_LDFLAGS="-Wl,-install_name,@rpath/libsharp.1.dylib"
	#SHARP_LDFLAGS="-Wl,-install_name,@rpath/libsharp.1.dylib,-compatibility_version,1,-current_version,1"
	# requires to add -Wl,-rpath,${HEALPIX}/lib to C++ flags ($CXXFLAGS)
    fi

}
#-------------
askSharpUserMisc () {

    echoLn "enter C compiler you want to use [$CC]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CC=$answer

    echo "enter options for C compiler: "
    echo "For optimal performance, this should include '-ffast-math'"
    echo "and '-march=native' (in gcc or clang, or your compiler's equivalent options)."
    echo "(If you are using gcc or clang and you want to produce a portable,"
    echo "high-performance library, you can also replace '-march=native'"
    echo "by '-DMULTIARCH'.)"
    #echoLn "[$CFLAGS]: "
    echoLn "[$SHARP_COPT]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && SHARP_COPT=$answer
}
#-------------
editSharpMakefile () {

    SHARPCDIR=${HEALPIX}/src/common_libraries/libsharp
    SHARPBLD=${SHARPCDIR}/build
    SHARPLIB=${SHARPLDIR}/libsharp.a
    SHARPMKF=${SHARPBLD}/Makefile
    echo
    echo "Running configure in ${SHARPBLD} ... "
    (\rm -rf ${SHARPBLD}; \
    mkdir ${SHARPBLD}; \
    cd ${SHARPBLD}; \
    CC="${CC}" CFLAGS="${SHARP_COPT}" LDFLAGS="${SHARP_LDFLAGS}" ${SHARPCDIR}/configure --prefix=${SHARPPREFIX} || crashAndBurn) || crashAndBurn
    echo "done."
    echo 
    echo "edit top Makefile for libsharp ..."
    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^ALL\(.*\) sharp-void \(.*\)|ALL\1 sharp-all \2|" |\
	${SED} "s|^TESTS\(.*\) sharp-void \(.*\)|TESTS\1 sharp-test \2|" |\
	${SED} "s|^CLEAN\(.*\) sharp-void \(.*\)|CLEAN\1 sharp-clean \2|" |\
	${SED} "s|^DISTCLEAN\(.*\) sharp-void \(.*\)|DISTCLEAN\1 sharp-distclean \2|" |\
	${SED} "s|^TIDY\(.*\) sharp-void \(.*\)|TIDY\1 sharp-tidy \2|" |\
	${SED} "s|^HEALPIX=.*$|HEALPIX	= $HEALPIX|" |\
        ${SED} "s|^SHARPLDIR.*|SHARPLDIR=${SHARPLDIR}|" |\
	${SED} "s|^SHARPBLD.*|SHARPBLD =${SHARPBLD}|" |\
	${SED} "s|^SHARPLIB.*|SHARPLIB =${SHARPLIB}|" |\
	${SED} "s|^SHARPMKF.*|SHARPMKF =${SHARPMKF}|" |\
	${SED} "s|^# sharp configuration.*|# sharp configuration: (rm -rf ${SHARPBLD}; mkdir ${SHARPBLD}; cd ${SHARPBLD}; CC=\"${CC}\" CFLAGS=\"${SHARP_COPT}\" LDFLAGS=\"${SHARP_LDFLAGS}\" ${SHARPCDIR}/configure --prefix=${SHARPPREFIX})|" > Makefile

#     mv -f Makefile Makefile_tmp
#     ${CAT} Makefile_tmp |\
#  	${SED} 's|^$(SHARPMKF).*$|$(SHARPMKF):\\
# 	mkdir -p $(SHARPBLD); \\\ \\
# 	\$(RM) -r  $(SHARPBLD)/*; \\\ \\
# 	cd $(SHARPBLD); \\\\ \\
# 	CC=${CC} CFLAGS=${SHARP_COPT} LDFLAGS=${SHARP_LDFLAGS} ${SHARPCDIR}/configure --prefix=${SHARPPREFIX}; \\\ \\
# 	cd \$(HEALPIX)|' > Makefile

    echo " done."
    edited_makefile=1
}
#-------------
Sharp_config () {

    setSharpDefaults
    askSharpUserMisc
    editSharpMakefile
}
#-------------
Sharp_install () {
    echo
    echo "Compiling and installing libsharp ..."
    makejobs=''
    ngnumake=`${MAKE} -v 2>&1   | ${GREP} GNU | ${WC} -l`
    [ $ngnumake != 0 ] && makejobs='-j'  # parallel make for GNUmake
    ( cd ${SHARPBLD} && ${MAKE} ${makejobs} install || crashAndBurn )
    echo "done."
}

#-------------
test_Sharp () {
    sharp_configured=0
    if [ -s Makefile ] ; then
	sharp_configured=`${CAT} Makefile | ${GREP} -c "^ALL\(.*\) sharp-all"`
    fi
}

#=====================================
#=========== C++ package ===========
#=====================================
#   setCppDefaults: defaults variables for C++
#-------------
setCppDefaults () {

    CFLAGS="-O3 -fopenmp"

    CXXBLD=${CXXDIR}/build
    CXX_LDFLAGS=""
    CXXFLAGS="${CXXFLAGS--O3 -fopenmp}" # -O3 -fopenmp unless already defined
    if [ "${OS}" = "Darwin" -a `isTrue ${USE_ATRPATH}` -eq 1 ]; then
	CXX_LDFLAGS="-Wl,-install_name,@rpath/libhealpix_cxx.2.dylib"
	# requires to add -Wl,-rpath,${HEALPIX}/lib to C++ flags ($CXXFLAGS)
	CXXFLAGS="${CXXFLAGS} -Wl,-rpath,${HEALPIX}/lib"
    fi

}
#-------------
askCppUserMisc () {

    echoLn "enter C compiler you want to use [$CC]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CC=$answer

    echoLn "enter options for C compiler [$CFLAGS]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CFLAGS=$answer

    echoLn "enter C++ compiler you want to use [$CXX]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CXX=$answer

    echoLn "enter options for C++ compiler [$CXXFLAGS]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CXXFLAGS=$answer


    #echo "${FITSINC} ${FITSDIR}"
    CFITSIO_CFLAGS=""
    echoLn "enter path to fitsio.h"
    if [ "$FITSINC" = "/usr/local/include" ]; then
	echoLn ", in case of blank answer, autotools will look in standard locations []: "
	read answer
    else
	echoLn " [$FITSINC]: "
	read answer
	[ "x$answer" = "x" ] && CFITSIO_CFLAGS="-I${FITSINC}"
    fi
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CFITSIO_CFLAGS="-I$answer"

    CFITSIO_LIBS=""
    echoLn "enter path to cfitsio library"
    if [ "$FITSDIR" = "/usr/local/lib" ]; then
	echoLn ", in case of blank answer, autotools will look in standard locations []: "
	read answer
    else
	echoLn "[$FITSDIR]: "
	read answer
	[ "x$answer" = "x" ] && CFITSIO_LIBS="-L${FITSDIR} -lcfitsio"
    fi
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CFITSIO_LIBS="-L$answer -lcfitsio"
    echo $CFITSIO_CFLAGS  $CFITSIO_LIBS
}
#-------------
editCppMakefile () {

    echo
    echo "Running configure in ${CXXBLD} ... "
    (\rm -rf ${CXXBLD}; \
    mkdir ${CXXBLD}; \
    cd ${CXXBLD}; \
    CC="${CC}" CFLAGS="${CFLAGS}" CXX="${CXX}" CXXFLAGS="${CXXFLAGS}" LDFLAGS="${CXX_LDFLAGS}" SHARP_CFLAGS="-I${HEALPIX}/include" SHARP_LIBS="-L${HEALPIX}/lib -lsharp" CFITSIO_CFLAGS="${CFITSIO_CFLAGS}" CFITSIO_LIBS="${CFITSIO_LIBS}" ${CXXDIR}/configure --prefix=${CXXPREFIX} || crashAndBurn)
    echo "done"

    echo 
    echo "edit top Makefile for C++ ..."
    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^ALL\(.*\) cpp-void \(.*\)|ALL\1 cpp-all \2|" |\
	${SED} "s|^TESTS\(.*\) cpp-void \(.*\)|TESTS\1 cpp-test \2|" |\
	${SED} "s|^CLEAN\(.*\) cpp-void \(.*\)|CLEAN\1 cpp-clean \2|" |\
	${SED} "s|^DISTCLEAN\(.*\) cpp-void \(.*\)|DISTCLEAN\1 cpp-distclean \2|" |\
	${SED} "s|^TIDY\(.*\) cpp-void \(.*\)|TIDY\1 cpp-tidy \2|" |\
	${SED} "s|^HEALPIX=.*$|HEALPIX	= $HEALPIX|" |\
	${SED} "s|^CXXBLD.*|CXXBLD=${CXXBLD}|" > Makefile

    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^# C++ configuration.*$|# C++ configuration: (rm -rf ${CXXBLD} ; mkdir ${CXXBLD} ; cd ${CXXBLD}; CC=\"${CC}\" CFLAGS=\"${CFLAGS}\" CXX=\"${CXX}\" CXXFLAGS=\"${CXXFLAGS}\"  LDFLAGS=\"${CXX_LDFLAGS}\" SHARP_CFLAGS=\"-I${HEALPIX}/include\" SHARP_LIBS=\"-L${HEALPIX}/lib -lsharp\" CFITSIO_CFLAGS=\"${CFITSIO_CFLAGS}\" CFITSIO_LIBS=\"${CFITSIO_LIBS}\" ${CXXDIR}/configure --prefix=${CXXPREFIX})|"> Makefile

#     mv -f Makefile Makefile_tmp
#     ${CAT} Makefile_tmp |\
# 	${SED} "s|^cpp-config.*$|cpp-config: sharp-all\\
# 	mkdir -p \$(CXXBLD); \\\ \\
# 	\$(RM) -r  \$(CXXBLD)/*; \\\ \\
# 	cd \$(CXXBLD); \\\\ \\
# 	CC=\"${CC}\" CFLAGS=\"${CFLAGS}\" CXX=\"${CXX}\" CXXFLAGS=\"${CXXFLAGS}\"  SHARP_CFLAGS=\"-I${HEALPIX}/include\" SHARP_LIBS=\"-L${HEALPIX}/lib -lsharp\" CFITSIO_CFLAGS=\"${CFITSIO_CFLAGS}\" CFITSIO_LIBS=\"${CFITSIO_LIBS}\" ${CXXDIR}/configure --prefix=${CXXPREFIX}; \\\ \\
# 	cd \$(HEALPIX)|" > Makefile
    echo " done."
    edited_makefile=1
}
#-------------
generateConfCppFile () {
	echo "Generating $HPX_CONF_CPP"

    	echo "# C++ configuration for HEALPix `date`" > $HPX_CONF_CPP

    case $SHELL in
    csh|tcsh)
	${CAT} <<EOF >>$HPX_CONF_CPP
setenv PATH \${HEALPIX}/bin:\${PATH}
EOF
	;;
    sh|ksh|bash|zsh)
	${CAT} <<EOF >>$HPX_CONF_CPP
PATH="\${HEALPIX}/bin:\${PATH}"
export PATH
EOF
	;;
    *)
	echo "Shell $SHELL not supported yet."
	${RM} $HPX_CONF_CPP
	;;
    esac

}
#-------------
Cpp_config () {
    test_Sharp
    if [ $sharp_configured = "0" ] ; then
	echo "Configuring the libsharp library first, since this is a dependency:"
	Sharp_config
	Sharp_install
	echo
	echo "Now configuring Healpix C++ itself:"
    fi
    HPX_CONF_CPP=$1
    setCppDefaults
    askCppUserMisc
    generateConfCppFile
    #installCppPackage
    editCppMakefile
    [ $NOPROFILEYET = 1 ] && installProfile
}

#=====================================
#=========== healpy Python package ===========
#=====================================
# #-------------
# Healpy_config () {
#     # CFITSIO: make a first guess
#     LIBFITS=cfitsio
#     fullPath FITSDIR
#     guess2=`${DIRNAME} ${FITSDIR}`
#     guess3=`${DIRNAME} ${FITSINC}`
#     findFITSPrefix $FITSDIR $FITSINC ${guess2} ${guess3}
#     # ask user
#     echo "Enter directory prefix for CFitsio"
#     echoLn " ie containing lib/libcfitsio.* and include/fitsio.h ($FITSPREFIX): "
#     read answer
#     [ "x$answer" != "x" ] && FITSPREFIX=$answer
#     # double check
#     inc="${FITSPREFIX}/include/fitsio.h"
#     if [ ! -r $inc ]; then
# 	echo "error: cfitsio include file $inc not found"
# 	crashAndBurn
#     fi
#     # apply
#     editHealpyMakefile

#     # update paths for C and C++
#     FITSDIR=${FITSPREFIX}/lib
#     FITSINC=${FITSPREFIX}/include

# }
# #-------------
# editHealpyMakefile () {


#     echoLn "edit top Makefile for Python (healpy) ..."

#     mv -f Makefile Makefile_tmp
#     ${CAT} Makefile_tmp |\
# 	${SED} "s|^ALL\(.*\) healpy-void\(.*\)|ALL\1 healpy-all \2|" |\
# 	${SED} "s|^TESTS\(.*\) healpy-void\(.*\)|TESTS\1 healpy-test \2|" |\
# 	${SED} "s|^CLEAN\(.*\) healpy-void\(.*\)|CLEAN\1 healpy-clean \2|" |\
# 	${SED} "s|^DISTCLEAN\(.*\) healpy-void\(.*\)|DISTCLEAN\1 healpy-distclean \2|" |\
# 	${SED} "s|^TIDY\(.*\) healpy-void\(.*\)|TIDY\1 healpy-tidy \2|" > Makefile

#     echo " done."
#     edited_makefile=1

# }



#-------------
Healpy_config () {  # for healpy 1.7.0

    HPY_PYTHON="$PYTHON"
    tmpfile=to_be_removed
    HPY_SETUP='setup.py' # default setup
    HPY_SETUP2='setup2.py' # backup setup
    HPY_DIR='src/healpy/'

    # ask for python command
    echoLn "Enter python command [$HPY_PYTHON] "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && HPY_PYTHON="$answer"

    # test python version number
    ${HPY_PYTHON} --version 1> ${tmpfile} 2>&1
    #python_version=`${CAT} ${tmpfile} | ${AWK} '{print \$NF}'` # current version (last field)
    python_version=`${CAT} ${tmpfile} | ${AWK} '{print \$2}'` # current version (2nd field)
    #python_reqrd="2.4" # minimal version supported
    #python_reqrd="2.6" # minimal version supported
    python_reqrd="2.7" # minimal version supported (1.12.8)
    p_v1=`echo ${python_version} | ${AWK} '{print $1*10}'`
    p_v2=`echo ${python_reqrd}   | ${AWK} '{print $1*10}'`
    ${RM} ${tmpfile}
    if [ ${p_v1} -lt ${p_v2} ]; then
	echo
	echo "Python version (${python_version}) must be >= ${python_reqrd}"
	echo
	crashAndBurn
    fi


#     # special treatement for MacOSX
#     if [ "${OS}" = "Darwin" ]; then
# 	# find out compiler and options used by python (and therefore healpy in setup.py)
# 	HPY_CC=`${HPY_PYTHON}   -c "from distutils.sysconfig import get_config_var ; print get_config_var('CC')"`
# 	HPY_OPTS=`${HPY_PYTHON} -c "from distutils.sysconfig import get_config_var ; print get_config_var('CFLAGS')"`

# 	# test these options on a C code
# ${CAT} > ${tmpfile}.c <<EOF
# int main(){
# }
# EOF
# 	${HPY_CC} -E ${HPY_OPTS} ${tmpfile}.c -o ${tmpfile}.p  1>${DEVNULL} 2>&1
# 	# if test fails, create back up setup2.py and make sure it will be used
# 	if [ ! -e ${tmpfile}.p ] ; then
# 	    cat ${HPY_DIR}${HPY_SETUP} | \
# 		sed "s|'--disable-shared',|'--disable-shared', '--disable-dependency-tracking',|g" > \
# 		${HPY_DIR}${HPY_SETUP2}
# 	    HPY_SETUP=${HPY_SETUP2}
# 	fi
# 	# clean up
# 	${RM}  ${tmpfile}.*
#     fi

    echoLn "enter C compiler you want to use [$CC]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CC=$answer
    HPY_CC="${CC}"

    echoLn "enter C++ compiler you want to use [$CXX]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CXX=$answer
    HPY_CXX="${CXX}"

    # apply
    editHealpyMakefile

}
#-------------
editHealpyMakefile () {

    echo 
    echo "edit top Makefile for Python (healpy) ..."
    #echo ${HPY_SETUP}
    HPY_VARS="CC=\"${HPY_CC}\" CXX=\"${HPY_CXX}\" "
    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^HPY_SETUP.*$|HPY_SETUP    = ${HPY_SETUP}|" |\
	${SED} "s|^HPY_PYTHON.*$|HPY_PYTHON   = ${HPY_PYTHON}|" |\
	${SED} "s|^HPY_VARS.*$|HPY_VARS     = ${HPY_VARS}|" |\
	${SED} "s|^ALL\(.*\) healpy-void\(.*\)|ALL\1 healpy-all \2|" |\
	${SED} "s|^TESTS\(.*\) healpy-void\(.*\)|TESTS\1 healpy-test \2|" |\
	${SED} "s|^CLEAN\(.*\) healpy-void\(.*\)|CLEAN\1 healpy-clean \2|" |\
	${SED} "s|^DISTCLEAN\(.*\) healpy-void\(.*\)|DISTCLEAN\1 healpy-distclean \2|" |\
	${SED} "s|^HEALPIX=.*$|HEALPIX	= $HEALPIX|" |\
	${SED} "s|^TIDY\(.*\) healpy-void\(.*\)|TIDY\1 healpy-tidy \2|" > Makefile

    echo " done."
    edited_makefile=1

}
#=====================================
#=========== IDL package ===========
#=====================================
#-------------
askPaperSize () {

    echo
    echo "Please indicate what size of paper you will be mainly "
    echo "using for your Postscript printouts (eg: a4/letter)"
    echoLn "Enter choice [$papersize]       "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && papersize="$answer"
    echo "Note: the plot bounding box (BBox) rather than the paper size "
    echo "will be used when previewing a Postcript file"
}
#-------------
askPS () {

    echo
    echo  "Please indicate the Postscript previewer you want to use "
    echo  " (eg: gs, ghostview, gv, ggv, kghostview, evince)"
    echoLn "Enter choice [$ps_com]        "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && ps_com="$answer"
    ngs=`$ps_com -v 2>&1 | ${GREP} -i ghostscript | ${WC} -l`
    ngv=`$ps_com -v 2>&1 | ${GREP} -i gv | ${WC} -l`
    nghostview=`echo $ps_com | ${GREP} -i ghostview | ${WC} -l`

    if [ $ngs != 0 ] ; then
	ps_scom="gs"
    elif [ $ngv != 0 ] ; then
	ps_scom="gv"
    elif [ $nghostview != 0 ] ; then
	ps_scom="ghostview"
    else
	ps_scom=$ps_com
    fi
}
#-------------
askPDF () {

    echo
    echo  "Please indicate the PDF previewer you want to use "
    echo  " (eg: gv, xpdf, kpdf, open)"
    echoLn "Enter choice [$pdf_com]        "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && pdf_com="$answer"
    pdf_scom=$pdf_com
}
#-------------
askGif () {

    echo
    echo "Please indicate the program to be used to view "
    echo "the GIF and PNG files generated by Healpix (eg: display, open)"
    echo "Note that xv may not be able to deal with PNG files"
    echoLn "Enter choice [$gif_com]       "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && gif_com="$answer"
    nnet=`$gif_com -version 2>&1 | ${GREP} -i netscape | ${WC} -l`

    if [ $nnet != 0 ] ; then
	gif_scom="netscape"
    else
	gif_scom=$gif_com
    fi

}
#-------------
generateProIdlFile () {
    echo
    echo "* Generating $previewfile"
    ${CAT} <<EOF > $previewfile
; ----------------------------------------
;
; File generated automatically by $0
; in $HEALPIX
; on `date`
; Used by preview_file.pro
;
; ----------------------------------------
papersize = '$papersize'
media = '$media'
ps_com = '$ps_com'
ps_scom = '$ps_scom'
pdf_com = '$pdf_com'
pdf_scom = '$pdf_scom'
gif_com  = '$gif_com'
gif_scom = '$gif_scom'
settings = 'user'
;-----------------------------------------
EOF
}
#-------------
generateConfIdlFile () {
    #SHELL=`${BASENAME} ${SHELL-/bin/sh}`
    echo
    echo "* Generating $HPX_CONF_IDL"
    echo "containing:"


    echo "# IDL configuration for HEALPix `date`" > $HPX_CONF_IDL

    case $SHELL in
    csh|tcsh)
	${CAT} <<EOF >>$HPX_CONF_IDL
# back up original IDL config, or give default value
if (\$?IDL_PATH) then
    setenv OIDL_PATH    "\${IDL_PATH}"
else
    setenv OIDL_PATH    "<IDL_DEFAULT>"
endif
if (\$?IDL_STARTUP) then
    setenv OIDL_STARTUP "\${IDL_STARTUP}"
else
    setenv OIDL_STARTUP
endif
# create Healpix IDL config, and return to original config after running Healpix-enhanced IDL
setenv HIDL_PATH  "+\${HEALPIX}/src/idl:\${OIDL_PATH}"
setenv HIDL_STARTUP \${HEALPIX}/src/idl/HEALPix_startup
alias  hidl    'setenv IDL_PATH \${HIDL_PATH} ; setenv IDL_STARTUP \${HIDL_STARTUP} ; idl   ; setenv IDL_PATH \${OIDL_PATH} ; setenv IDL_STARTUP \${OIDL_STARTUP}'
alias  hidlde  'setenv IDL_PATH \${HIDL_PATH} ; setenv IDL_STARTUP \${HIDL_STARTUP} ; idlde ; setenv IDL_PATH \${OIDL_PATH} ; setenv IDL_STARTUP \${OIDL_STARTUP}'
EOF
	;;
#    sh)
#	${CAT} <<EOF >>$HPX_CONF_IDL
#HIDL_PATH=+$HEAPIX/src/idl:$IDL_PATH
#export HIDL_PATH
#EOF
#	;;
    sh|ksh|bash|zsh)
	${CAT} <<EOF >>$HPX_CONF_IDL
# make sure IDL related variables are global
export IDL_PATH IDL_STARTUP
# back up original IDL config, or give default value
OIDL_PATH="\${IDL_PATH-<IDL_DEFAULT>}"
OIDL_STARTUP="\${IDL_STARTUP}"
# create Healpix IDL config, and return to original config after running Healpix-enhanced IDL
HIDL_PATH="+\${HEALPIX}/src/idl:\${OIDL_PATH}"
HIDL_STARTUP="\${HEALPIX}/src/idl/HEALPix_startup"
alias hidl="IDL_PATH=\"\${HIDL_PATH}\"   ; IDL_STARTUP=\${HIDL_STARTUP} ; idl   ; IDL_PATH=\"\${OIDL_PATH}\" ; IDL_STARTUP=\${OIDL_STARTUP} "
alias hidlde="IDL_PATH=\"\${HIDL_PATH}\" ; IDL_STARTUP=\${HIDL_STARTUP} ; idlde ; IDL_PATH=\"\${OIDL_PATH}\" ; IDL_STARTUP=\${OIDL_STARTUP} "
EOF
	;;
    *)
	echo "Shell $SHELL not supported yet."
	${RM} $HPX_CONF_IDL
	crashAndBurn
	;;
    esac

    ${CAT} $HPX_CONF_IDL
    echo
    echo
}
#-------------
generateConfGdlFile () {
    echo
    echo "* Generating $HPX_CONF_GDL"
    echo "containing:"


    echo "# GDL configuration for HEALPix `date`" > $HPX_CONF_GDL

    case $SHELL in
    csh|tcsh)
	${CAT} <<EOF >>$HPX_CONF_GDL
# back up original GDL config, or give default value
if (\$?GDL_PATH) then
    setenv OGDL_PATH    "\${GDL_PATH}"
else
    setenv OGDL_PATH#    "<GDL_DEFAULT>"
endif
if (\$?GDL_STARTUP) then
    setenv OGDL_STARTUP "\${GDL_STARTUP}"
else
    setenv OGDL_STARTUP
endif
# create Healpix GDL config, and return to original config after running Healpix-enhanced GDL
setenv HGDL_PATH  "+\${HEALPIX}/src/idl:\${OGDL_PATH}"
setenv HGDL_STARTUP \${HEALPIX}/src/idl/HEALPix_startup
alias  hgdl    'setenv GDL_PATH \${HGDL_PATH} ; setenv GDL_STARTUP \${HGDL_STARTUP} ; gdl   ; setenv GDL_PATH \${OGDL_PATH} ; setenv GDL_STARTUP \${OGDL_STARTUP}'
alias  hgdlde  'setenv GDL_PATH \${HGDL_PATH} ; setenv GDL_STARTUP \${HGDL_STARTUP} ; gdlde ; setenv GDL_PATH \${OGDL_PATH} ; setenv GDL_STARTUP \${OGDL_STARTUP}'
EOF
	;;
    sh|ksh|bash|zsh)
	${CAT} <<EOF >>$HPX_CONF_GDL
# make sure GDL related variables are global
export GDL_PATH GDL_STARTUP
# back up original GDL config, or give default value
OGDL_PATH="\${GDL_PATH-<GDL_DEFAULT>}"
OGDL_STARTUP="\${GDL_STARTUP}"
# create Healpix GDL config, and return to original config after running Healpix-enhanced GDL
HGDL_PATH="+\${HEALPIX}/src/idl:\${OGDL_PATH}"
HGDL_STARTUP="\${HEALPIX}/src/idl/HEALPix_startup"
alias hgdl="GDL_PATH=\"\${HGDL_PATH}\"   ; GDL_STARTUP=\${HGDL_STARTUP} ; gdl   ; GDL_PATH=\"\${OGDL_PATH}\" ; GDL_STARTUP=\${OGDL_STARTUP} "
alias hgdlde="GDL_PATH=\"\${HGDL_PATH}\" ; GDL_STARTUP=\${HGDL_STARTUP} ; gdlde ; GDL_PATH=\"\${OGDL_PATH}\" ; GDL_STARTUP=\${OGDL_STARTUP} "
EOF
	;;
    *)
	echo "Shell $SHELL not supported yet."
	${RM} $HPX_CONF_GDL
	crashAndBurn
	;;
    esac

    ${CAT} $HPX_CONF_GDL
    echo
    echo
}
#-------------
generateConfFlFile () {
    echo
    echo "* Generating $HPX_CONF_FL"
    echo "containing:"

    echo "# FL configuration for HEALPix `date`" > $HPX_CONF_FL

    if [ "x$FL_DIR" != "x" ] ; then
	if [ "${OS}" = "Darwin" ]; then
	    fl_path="\${FL_DIR}/Contents/MacOS/fl64_cmd"
	    flde_path="\${FL_DIR}/Contents/MacOS/fl64_gui"
	else
	    fl_path=`${LS} ${FL_DIR}/bin/fl??_cmd`
	    flde_path=`${LS} ${FL_DIR}/bin/fl??_gui`
	fi
    else
	whereisCmd fl64_cmd fl32_cmd
	fl_path=${CMD}
	whereisCmd fl64_gui fl32_gui
	flde_path=${CMD}
	if [ "${OS}" = "Darwin" -a "x${fl_path}" = "x" ]; then
	    fl_path='/Applications/fl.app/Contents/MacOS/fl64_cmd'
	    flde_path='/Applications/fl.app/Contents/MacOS/fl64_gui'
	fi
    fi

    case $SHELL in
    csh|tcsh)
	${CAT} <<EOF >>$HPX_CONF_FL
# back up original FL config, or give default value
if (\$?IDL_PATH) then
    setenv OFL_PATH    "\${IDL_PATH}"
else
    setenv OFL_PATH
endif
if (\$?IDL_STARTUP) then
    setenv OFL_STARTUP "\${IDL_STARTUP}"
else
    setenv OFL_STARTUP
endif
# create Healpix FL config, and return to original config after running Healpix-enhanced FL
setenv HFL_PATH  "+\${HEALPIX}/src/idl:\${OFL_PATH}"
setenv HFL_STARTUP \${HEALPIX}/src/idl/HEALPix_startup
alias  hfl    'setenv IDL_PATH \${HFL_PATH} ; setenv IDL_STARTUP \${HFL_STARTUP} ; ${fl_path}   ; setenv IDL_PATH \${OFL_PATH} ; setenv IDL_STARTUP \${OFL_STARTUP}'
alias  hflde  'setenv IDL_PATH \${HFL_PATH} ; setenv IDL_STARTUP \${HFL_STARTUP} ; ${flde_path} ; setenv IDL_PATH \${OFL_PATH} ; setenv IDL_STARTUP \${OFL_STARTUP}'
EOF
	;;
    sh|ksh|bash|zsh)
	${CAT} <<EOF >>$HPX_CONF_FL
# make sure FL related variables are global
export FL_PATH FL_STARTUP
# back up original FL config, or give default value
OFL_PATH="\${FL_PATH-<FL_DEFAULT>}"
OFL_STARTUP="\${FL_STARTUP}"
# create Healpix FL config, and return to original config after running Healpix-enhanced FL
HFL_PATH="+\${HEALPIX}/src/idl:\${OFL_PATH}"
HFL_STARTUP="\${HEALPIX}/src/idl/HEALPix_startup"
alias hfl="FL_PATH=\"\${HFL_PATH}\"   ; FL_STARTUP=\${HFL_STARTUP} ; ${fl_path}   ; FL_PATH=\"\${OFL_PATH}\" ; FL_STARTUP=\${OFL_STARTUP} "
alias hfde="FL_PATH=\"\${HFL_PATH}\" ; FL_STARTUP=\${HFL_STARTUP} ; {flde_path} ; FL_PATH=\"\${OFL_PATH}\" ; FL_STARTUP=\${OFL_STARTUP} "
EOF
	;;
    *)
	echo "Shell $SHELL not supported yet."
	${RM} $HPX_CONF_FL
	crashAndBurn
	;;
    esac

    ${CAT} $HPX_CONF_FL
    echo
    echo
}
#-------------
setIdlDefaults () {

    papersize="${papersize-letter}" # letter   unless already defined
    media=" "
    ps_com="${ps_com-gv}"           # gv       unless already defined
    [ "${OS}" = "Darwin" ]  && pdf_com="${pdf_com-open}"  || pdf_com="${pdf_com-xpdf}" # open or xpdf unless already defined
    #gif_com="${gif_com-netscape}"   # netscape unless already defined
    [ "${OS}" = "Darwin" ]  && gif_com="${gif_com-open}"  || gif_com="${gif_com-display}" # open or display unless already defined
    previewfile=${HEALPIX}/src/idl/visu/idl_default_previewer.pro

#     # if IDL_PATH is undefined, then set it to +IDL_DIR
#     if [ -z "${IDL_PATH}" ] ; then
# 	if [ -z "${IDL_DIR}" ] ; then
# 	    echo
# 	    echo "IDL_DIR is undefined. Make sure it is defined before proceeding."
# 	    echo
# 	    crashAndBurn
# 	fi
# 	IDL_PATH="+${IDL_DIR}"
#     fi
}
#-------------
idl_config () {

    HPX_CONF_IDL=$1
    HPX_CONF_GDL=$2
    HPX_CONF_FL=$3
    setIdlDefaults
    askPaperSize
    askPS
    askGif
    askPDF
    generateProIdlFile
    generateConfIdlFile
    generateConfGdlFile
    generateConfFlFile
    [ $NOPROFILEYET = 1 ] && installProfile
}

#=====================================
#=========== F90 package ===========
#=====================================
#
#   setF90Defaults: set default values of variables
#   sun_modules : test weither the Sun compiler creates modules ending with .M or .mod
#   ifc_modules : test weither the IFC compiler creates .d or .mod (version7) modules
#   checkF90Fitsio:        check that CFITSIO library contains Fortran wrapper
#   checkFitsioCurl:       check if   CFITSIO library is linked to libcurl
#   checkF90FitsioLink:    check that CFITSIO library links to Fortran test code
#   checkF90FitsioVersion: check that CFITSIO library is recent enough
#   checkCParall: check that C compiler supports OpenMP
#   GuessF90Compiler: tries to guess compiler from operating system
#   askOpenMP: ask user for compilation of OpenMP source files
#   askF90PIC: ask user for -fPIC compilation of code
#   patchF90: all patches to apply to F90 and/or C compilers
#   IdentifyCParallCompiler: identify C compiler used for parallel compilation of SHT routines
#   IdentifyF90Compiler : identify Non native f90 compiler
#   add64bitF90Flags: add 64 bit flags to F90 (and C) compiler
#   countF90Bits: count number of addressing bits in code produced by F90 compiler
#   countCBits:   count number of addressing bits in code produced by C   compiler
#   checkF90Compilation: check that F90 compiler actually works
#   checkF90LongLong: check that F90 support 8 byte integers
#   askUserF90:  ask user the f90 compiler command
#   showDefaultDirs: show default directories
#   updateDirs: update those directories
#   showActualDirs: show actual directories
#   askUserMisc:  ask user to confirm or change various defaults
#   askPgplot: ask if user wants to link with PGPlot
#   editF90Makefile: create makefile from template
#   generateConfF90File: generates configuration file for F90
#   offerF90Compilation: propose to perform F90 compilation
#   f90_shared: deal with shared F90 library
#   writeF90pkgconfigFile: writes pkg-config (.pc) file for F90 library
#   f90_config: top routine for F90
#
#-------------
setF90Defaults () {
    FFLAGS="-I\$(F90_INCDIR)"
    CFLAGS="-O3 -std=c99"  # OK for gcc, icc and clang
    #f90_LDFLAGS="-L\$(F90_LIBDIR) -L\$(FITSDIR) -lhealpix -lhpxgif -lsharp -l\$(LIBFITS)"
    F90_LDFLAGS="-L\$(F90_LIBDIR) -L\$(FITSDIR) -L\$(SHARPLDIR) -lhealpix -lhpxgif -lsharp -l\$(LIBFITS)"
    F90_BINDIR="./bin"
    F90_INCDIR="./include"
    F90_LIBDIR="./lib"
    F90_BUILDDIR="./build"
    F_DIRSUFF="${F_DIRSUFF}"
    MOD="mod"
    FPP="-D"
    PARALL=""
    PRFLAGS=""
    F90_AR="ar rv"
    FTYPE=""
    PPFLAGS=""
    FF64=""
    CF64=""
    AR64=""
    PGFLAG=""
    PGLIBS=""
    PGLIBSDEF="-L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11"
    F90_WLRPATH="" # to add a directory to the (linker) runtime library search path
    F90PIC="-fPIC"
    F90_LIBSUFFIX=".a" # static library by default

    echo "you seem to be running $OS"

    case $OS in
	AIX)
	    FC="${FC-xlf90_r}";; # xlf90_r unless already defined
	Linux)
	    FC="${FC-gfortran}";;
	Darwin)
	    FC="${FC-gfortran}";;
	SUPER-UX)
	    FC="${FC-f90}";; # f90 unless already defined
	*)
	    FC="${FC-f90}";; # f90 unless already defined
    esac

    FCNAME="$OS Native compiler"
}

# -----------------------------------------------------------------

sun_modules () {
tmpfile=to_be_removed
suffix=.f90
${CAT} > ${tmpfile}${suffix} << EOF
   module ${tmpfile}
       integer :: i
   end module ${tmpfile}
EOF
   $FC -c ${tmpfile}${suffix} -o ${tmpfile}.o

   if test -s ${tmpfile}.M  ; then
       MOD="M"
   else
       MOD="mod"
   fi

   ${RM} ${tmpfile}.*
}

# -----------------------------------------------------------------
ifc_modules () {
tmpfile=to_be_removed
suffix=.f90
${CAT} > ${tmpfile}${suffix} << EOF
   module ${tmpfile}
       integer :: i
   end module ${tmpfile}
EOF
   $FC -c ${tmpfile}${suffix} -o ${tmpfile}.o 2> ${DEVNULL}

   if test -s ${tmpfile}.d  ; then
    # version 5 and 6 of ifc
	echo "This version of ifc is no longer supported"
	echo "use a more recent version (7 or higher)"
	crashAndBurn
#         IFCMOD="d"
#         IFCINC="-cl,\$(HEALPIX)/include/list.pcl"
#         IFCVERSION="ifcold"
   else
       IFCMOD="mod"
#         IFCINC="-I\$(HEALPIX)/include"
       IFCINC="-I\$(F90_INCDIR)"
#         IFCVERSION="ifnew"
       IFCVERSION="ifcnew"
   fi

   ${RM}  ${tmpfile}.*
   ${RM}  TO_BE_REMOVED.*
}

#----------
par2curly () {
    # replace () with {}
    echo `echo "$*" | ${SED} "s|(|{|g" | sed "s|)|}|g"`
}
#----------
checkF90Fitsio () {
    cfitsiolib=$1
    sanity=`${NM} ${cfitsiolib} 2> ${DEVNULL} | ${GREP} read | ${GREP} T | ${WC} -l` # make sure that nm, grep and wc are understood
    if [ $sanity -gt 0 ] ; then
	check=`${NM} ${cfitsiolib} 2> ${DEVNULL} | ${GREP} ftgkey | ${GREP} T | ${WC} -l` # count ftgkey definition
	if [ $check -eq 0 ] ; then
	    echo
	    echo "*WARNING*: the cfitsio library ${cfitsiolib}"
	    echo "does not seem to include the Fortran interface;"
	    echo "this will prevent compilation of the Healpix F90 package."
	    echo
	    echo "When installing cfitsio, make sure that a Fortran compiler is known"
	    echo "to the cfitsio configure script."
	    echo
	fi
    fi

}
#----------
checkFitsioCurl () {
    cfitsiolib=$1
    if [ -z ${CFITSIOCURL} ] ; then
	CFITSIOCURL=" "
	sanity=`${NM} ${cfitsiolib} 2> ${DEVNULL} | ${GREP} read | ${GREP} T | ${WC} -l` # make sure that nm, grep and wc are understood
	if [ $sanity -gt 0 ] ; then
	    check=`${NM} ${cfitsiolib} 2> ${DEVNULL} | ${GREP} curl_ | ${WC} -l` # count curl calls
	    if [ $check -gt 1 ] ; then # curl_verbose does not count (2021-01-26)
		CFITSIOCURL="-lcurl"
 		CURLCONFIG='curl-config' # use curl-config to find non /usr/lib*/ installations
 		testcurl=`${WHEREIS} ${CURLCONFIG}`
		#echo "TESTCURL:  $testcurl"
		if [ "x${testcurl}" != "x" -a "x{testcurl}" != "x${CURLCONFIG}:" ] ; then
		    #pathcurl=`${CURLCONFIG} --prefix`
		    #pcloc=${pathcurl:0:2} # bashism !!
		    pathcurl=`${CURLCONFIG} --libs`
		    pcloc=`echo $pathcurl | awk '{ print substr($0, 1, 2); }'`
		    if [ "pcloc}" = "-L" ] ; then
 			PATHCURL="${pathcurl}"
 			CFITSIOCURL="${PATHCURL} ${CFITSIOCURL}"
 		    fi
		fi

		# adding -lcurl is necessary if the cfitsio library is static
		#echo "WARNING: his version of CFITSIO must be linked with libcurl (flag -lcurl will be added)"
	    fi
	fi
    fi
}
# ----------------
checkF90FitsioLink () {
# check that F90 routines can link with F90-fitsio wrappers
# requires compilation of F90 code
    tmpfile=to_be_removed
    suffix=.f90
    # write simple program to link with fitsio
cat > ${tmpfile}${suffix} << EOF
    program needs_fitsio
	character(len=6) :: string='abcdef'
	call ftupch(string)
    end program needs_fitsio
EOF
    # compile and link
    FFLAGS_=`par2curly  ${FFLAGS}`  # $() -> ${}
    FFLAGS_=`eval echo ${FFLAGS_}` # ${} -> value
    ${FC} ${FFLAGS_}  ${tmpfile}${suffix} -o ${tmpfile}.x -L${FITSDIR} -l${LIBFITS} ${CFITSIOCURL} ${F90_WLRPATH}
    # test
    if [ ! -s ${tmpfile}.x ]; then
	echo
	echo "F90 codes do not link correctly with ${FITSDIR}/lib${LIBFITS}.a"
	echo "Check that in the cfitsio library:"
	echo " - the Fortran wrappers were correctly compiled, and"
	echo " - the library (C routines and F90 wrappers) was compiled "
	echo "   with a number of bits compatible with ${FC} ${FFLAGS}"
	if [ "x${CFITSIOCURL}x" != "x x" ] ; then
	    echo "Make sure that libcurl is properly installed"
	    echo "(see HEALPix installation documentation for more details)"
	fi
	crashAndBurn
    fi

    # clean up
    ${RM} ${tmpfile}.*


}
# ----------------
checkF90FitsioVersion () {
# check that FITSIO version is recent enough
# requires compilation of F90 code
    tmpfile=./to_be_removed # do not forget ./ to allow execution
    suffix=.f90
    # write simple test program
cat > ${tmpfile}${suffix} << EOF
    program date_fitsio
	real:: version
	call ftvers(version)
	write(*,'(f5.3)') version
    end program date_fitsio
EOF
    # compile and link
    FFLAGS_=`par2curly  ${FFLAGS}`  # $() -> ${}
    FFLAGS_=`eval echo ${FFLAGS_}` # ${} -> value
    ${FC} ${FFLAGS_}  ${tmpfile}${suffix} -o ${tmpfile}.x -L${FITSDIR} -l${LIBFITS} ${CFITSIOCURL} ${F90_WLRPATH_}
    #CFITSIOVREQ="3.14"            # required  version of CFITSIO (in Healpix 3.00)
    CFITSIOVREQ="3.20"            # required  version of CFITSIO (in Healpix 3.30)
    CFITSIOVREC="3.44"            # recommended  version of CFITSIO (according to NASA)
    # run if executable
    if [ -x ${tmpfile}.x ]; then
	CFITSIOVERSION=`${tmpfile}.x` || CFITSIOVERSION=-1 # available version of CFITSIO
	v1=`echo ${CFITSIOVERSION} | ${AWK} '{print $1*1000}'` # multiply by 1000 to get integer
	v2=`echo ${CFITSIOVREQ}    | ${AWK} '{print $1*1000}'`
	v3=`echo ${CFITSIOVREC}    | ${AWK} '{print $1*1000}'`
	if [ $v1 -lt 0   ]; then
	    echo
	    echo "The code compiled with"
	    echo "${FC} ${FFLAGS_}  ${tmpfile}${suffix} -o ${tmpfile}.x -L${FITSDIR} -l${LIBFITS} ${CFITSIOCURL} ${F90_WLRPATH_}"
	    echo "can not be executed."
	    echo "(missing shared/dynamic library ?)"
	    crashAndBurn
	fi
	if [ $v1 -lt $v2 ]; then
	    echo
	    echo "CFITSIO version in ${FITSDIR}/lib${LIBFITS}.a  is  $CFITSIOVERSION "
	    echo "CFITSIO >= ${CFITSIOVREQ} is expected for Healpix-F90"
	    echo
	    ${RM} ${tmpfile}.*
	    crashAndBurn
	fi
	if [ $v1 -lt $v3 ]; then
	    echo
	    echo "WARNING: CFITSIO version in ${FITSDIR}/lib${LIBFITS}.a  is  $CFITSIOVERSION "
	    echo "         CFITSIO >= ${CFITSIOVREC} is recommended by NASA security team"
	    echo
	    ${RM} ${tmpfile}.*
	fi
    else
	echo "Warning: unable to check that CFITSIO is recent enough (>= ${CFITSIOVREQ})"
	${RM} ${tmpfile}.*
    fi

    # clean up
    ${RM} ${tmpfile}.*


}

# -----------------------------------------------------------------

GuessF90Compiler () {
    case $OS in
	AIX)
	    IdentifyF90Compiler;;
	SunOS)
	    sun_modules
	    FFLAGS=`echo $FFLAGS | ${SED} "s/-I/-M/g"`
	    F90_LDFLAGS="$F90_LDFLAGS -lm -lnsl -lsocket"
	    OFLAGS="-fast";;
	IRIX*)
	    OS="IRIX"
	    F90_LDFLAGS="$F90_LDFLAGS -lm"
	    OFLAGS="-fast-O"
	    PRFLAGS="-mp";;
	Linux)
	    F90_AR="ar -rsv" # archive with index table
  	    OFLAGS="-O"
	    IdentifyF90Compiler;;
	Darwin)
  	    OFLAGS="-O"
	    F90_AR="libtool -static -s -o"  # archive with index table
	    IdentifyF90Compiler;;
	OSF*)
	    OS="OSF"
	    OFLAGS="-O5 -fast"
	    PRFLAGS="-omp";;
	SUPER-UX)
	    FC="f90"
	    FFLAGS="$FFLAGS"
	    OFLAGS="-C vopt"
	    CFLAGS="-C vopt"
	    FPP="-D"
	    PRFLAGS="-P openmp";;
	CYGWIN*)
	    OFLAGS="-O"
	    IdentifyF90Compiler;;
	MINGW*)
	    OFLAGS="-O"
	    IdentifyF90Compiler;;
	*)
	    echo "\"$OS\" is not supported yet"
	    crashAndBurn;;
    esac
    # override final OFLAGS with user provided F_OPT, if any
    [  -z "${F_OPT}"  ] && F_OPT="${OFLAGS}" || OFLAGS="${F_OPT}"
}

# -----------------------------------------------------------------

askOpenMP () {
    OpenMP=`isTrue ${F_PARAL}`
    echo " The Spherical Harmonics Transform (C and F90) routines used by "
    echo "synfast/anafast/smoothing/plmgen"
    echo "and some routines used by ud_grade and alteralm respectively"
    echo "have a parallel implementation (based on OpenMP)."
#     echo " It has been successfully tested on xlf (IBM), "
#     echo "gcc and ifort (Linux and/or MacOSX) compilers (systems)"
    echo "Do you want to use :"
    echo " 0) the standard serial implementation ?"
    echo " 1) the parallel implementation "
    echoLn "Enter choice                                      [$OpenMP]: "
    read answer
#    [ "x$answer" != "x" ] && OpenMP="$answer"
    [ {$OPenMP} = 0 -a `isTrue ${answer}`  = 1 ] && OpenMP=1
    [ {$OPenMP} = 1 -a `isFalse ${answer}` = 1 ] && OpenMP=0
    if [ $OpenMP = 1 ] ; then

	# deal with C and F90 flags
	IdentifyCParallCompiler
	if [ "x$PRFLAGS" != "x" ] ; then
	    FFLAGS="$FFLAGS $PRFLAGS"
	    if [ "x$PRCFLAGS" != "x"  ] ; then
		CFLAGS="$CFLAGS $PRCFLAGS"
	    fi
	else
	    echo "WARNING: Healpix+OpenMP not tested for  \"$FCNAME\" under \"$OS\" "
	    echo "Contact us (https://healpix.sourceforge.io/support.php) "
	    echo "if you already used OpenMP in this configuration."
	    echo "Will perform serial implementation of C and F90 routines instead."
	fi

    fi
}
# -----------------------------------------------------------------

askF90PIC () {
    DoF90PIC="1"
    echo " "
    echo " Do you want a Position Independent Compilation  (option  \"$F90PIC\") "
    echoLn "(recommended if the Healpix-F90 library is to be linked to external codes) (y|n) [y]: "
    read answer
    if [ "x$answer" = "xy"  -o "x$answer" = "xY"  -o "x$answer" = "x" ]; then
	if [ "x$F90PIC" != "x" ] ; then
	    # update FFLAGS
	    FFLAGS="$FFLAGS $F90PIC"
 	    # update CFLAGS
 	    CPIC="-fPIC" # hacked from setCDefaults
	    case $OS in
		AIX)
		    CPIC="-G"
		;;
	    esac
 	    CFLAGS="$CFLAGS $CPIC"
	else
	    echo "PIC compilation flag not known for  \"$FCNAME\" under \"$OS\" "
	    echo "standard static compilation will be performed"
	fi
    fi
}

# -----------------------------------------------------------------
patchF90 (){
# all patches to apply to F90 and/or C compilers

# F90 compiler: nothing!

# C compiler:
#  *  add -fno-tree-fre for GCC 4.4* versions (adapted from autoconf)
#  *  add -fno-tree-fre for GCC 4.3.4 version (detected at NERSC)
    GCCVERSION="`$CC -dumpversion 2>&1`"
    gcc44=`echo $GCCVERSION | grep -c '^4\.4'`
    if test $gcc44 -gt 0; then
	CFLAGS="$CFLAGS -fno-tree-fre"
    fi


    GCCLONGVERSION="`$CC --version 2>&1`"
    gcc434=`echo $GCCLONGVERSION | grep -c '4\.3\.4'`
    if test $gcc434 -gt 0; then
	CFLAGS="$CFLAGS -fno-tree-fre"
    fi

}
# -----------------------------------------------------------------
checkCParall () {
# check that C compiler actually support parallel (OpenMP) compilation
tmpfile=to_be_removed
suffix=.c

${CAT} > ${tmpfile}${suffix} <<EOF
#include <omp.h>
#include <stdio.h>
int main() {
#pragma omp parallel
printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}
EOF

${CC} ${PRCFLAGS} -c ${tmpfile}${suffix} -o ${tmpfile}.o  >  ${DEVNULL} 2>&1
if [ ! -s ${tmpfile}.o ]; then
    echo
    echo "WARNING: the C compiler "
    echo "   ${CC} ${PRCFLAGS} "
    echo "  does not support OpenMP."
    echo "  The code will not be as fast as it could be on multiprocessor architectures."
    PRCFLAGS=""
    ${RM} ${tmpfile}.*
    #crashAndBurn
fi
${RM} ${tmpfile}.*

}
# -----------------------------------------------------------------
IdentifyCParallCompiler () {
# add OpenMP flag for C compiler (currently only gcc, icc and clang-omp)
# http://openmp.org/wp/openmp-compilers/
    nicc=`$CC -V 2>&1          | ${GREP} -i intel | ${WC} -l`
    ngcc=`$CC --version 2>&1   | ${GREP} -i 'GCC' | ${GREP} -v pgcc | ${WC} -l`
    nclang=`$CC --version 2>&1 | ${GREP}  'clang' | ${WC} -l`
    #npgc=`$CC -V 2>&1          | ${GREP} -i portland | ${WC} -l` # portland C
    npgc=`$CC -V 2>&1          | ${GREP} PGI | ${WC} -l` # PGI C
    npath=`$CC -v 2>&1         | ${GREP} -i ekopath  | ${WC} -l` # pathscale EKOPath
    PRCFLAGS=""
    if [ $nicc != 0 ] ; then
	PRCFLAGS='-qopenmp' # -openmp-report0
    elif [ $ngcc != 0 ] ; then
	PRCFLAGS='-fopenmp'
    elif [ $nclang != 0 ] ; then
	PRCFLAGS='-fopenmp'
    elif [ $npgc != 0 ] ; then
	PRCFLAGS='-mp'
    elif [ $npath != 0 ] ; then
	PRCFLAGS='-mp'
    else
	echo "$CC: Unknown C compiler"
	echo "Enter flags for C compilation with OpenMP"
	read answer
	[ "x$answer" != "x" ] && PRCFLAGS="$answer"
    fi
    checkCParall
}
# -----------------------------------------------------------------
ExtendCFLAGS () {
tmp=$CFLAGS
for i in $1 ; do
    match=0
    for j in $tmp ; do
	[ "$i" = "$j" ] && match=1
    done
    [ "$match" = "0" ] && tmp="$tmp $i"
done
CFLAGS=$tmp
}
# -----------------------------------------------------------------

IdentifyCCompiler () {
#    ngcc=`$CC --version 2>&1   | ${GREP} '(GCC)'     | ${WC} -l` # gcc
    ngcc=`$CC --version 2>&1   | ${GREP} -i 'GCC'  | grep -v pgcc   | ${WC} -l` # gcc
    nicc=`$CC -V 2>&1          | ${GREP} -i intel    | ${WC} -l` # intel C compiler
    nclang=`$CC --version 2>&1 | ${GREP} clang       | ${WC} -l` # clang
    #npgc=`$CC -V 2>&1          | ${GREP} -i portland | ${WC} -l` # portland C
    npgc=`$CC -V 2>&1          | ${GREP} PGI | ${WC} -l` # PGI C
    npath=`$CC -v 2>&1         | ${GREP} -i ekopath  | ${WC} -l` # pathscale EKOPath
	ExtendCFLAGS "-I\$(HEALPIX)/include"
    if [ $ngcc != 0 ] ; then
	echo "$CC: GCC compiler"
	ExtendCFLAGS "-O3 -std=c99"
    elif [ $nicc != 0 ] ; then
	echo "$CC: Intel C compiler"
	ExtendCFLAGS "-O3 -std=c99"
    elif [ $nclang != 0 ] ; then
	echo "$CC: clang C compiler"
	ExtendCFLAGS "-O3 -std=c99"
    elif [ $npgc != 0 ] ; then
	#echo "$CC: PGI C compiler"
	echo "$CC: PGI C compiler"
	CFLAGS=`echo $CFLAGS | ${SED} 's|-std=c99||g'`
	ExtendCFLAGS "-O3 -c99"
    elif [ $npath != 0 ] ; then
	echo "$CC: Pathscale EKOPath C compiler"
	ExtendCFLAGS "-O3"
    else
	echo "$CC: unknown C compiler"
	ExtendCFLAGS "-O"
    fi
}
# -----------------------------------------------------------------

IdentifyF90Compiler () {
# For Linux and Darwin
# Lahey and Fujitsu still have to be tested
	DO_F90_SHARED=0 # do NOT know how to create a shared library
        nima=`$FC -V 2>&1 | ${GREP} -i imagine1 | ${WC} -l`
        nnag=`$FC -V 2>&1 | ${GREP} -i nagware  | ${WC} -l`
        nifc=`$FC -V 2>&1 | ${GREP} -i intel    | ${WC} -l`
        #npgf=`$FC -V 2>&1 | ${GREP} -i portland | ${WC} -l`
        npgf=`$FC -V 2>&1 | ${GREP} PGI | ${WC} -l`
	nlah=`$FC --version 2>&1 | ${GREP} -i lahey | ${WC} -l`
	nfuj=`$FC -V 2>&1 | ${GREP} -i fujitsu | ${WC} -l`
	nxlf=`$FC --help 2>&1 | ${HEAD} -15 | ${GREP} XL | ${WC} -l`
	nabs=`$FC -V 2>&1 | ${GREP} 'Pro Fortran' | ${WC} -l`
	ng95=`$FC -dumpversion 2>&1 | ${GREP} 'g95' | ${WC} -l`
	ngfortran=`$FC --version 2>&1 | ${GREP} 'GNU Fortran' | ${WC} -l`
	npath=`$FC -v 2>&1 | ${GREP} -i ekopath | ${WC} -l`
	nflang=`$FC --version | ${GREP} clang | ${WC} -l`
        if [ $nima != 0 ] ; then
                FCNAME="Imagine F compiler"
                FFLAGS="$FFLAGS -w -dusty -mismatch_all"
		echo "$FCNAME is not supported yet"
		crashAndBurn
        elif [ $nnag != 0 ] ; then
                FCNAME="NAGWare compiler"
		PPFLAGS="-fpp"
# compiler flags for very thorough checking. use for debugging
#                FFLAGS="$FFLAGS -strict95 -g -gline -C=all -u -colour"
# standard flags
                FFLAGS="$FFLAGS -strict95"
		FI8FLAG="-double" # change default INTEGER and FLOAT to 64 bits
		MODDIR="-mdir " # output location of modules
        elif [ $nifc != 0 ] ; then
		ifc_modules
                FCNAME="Intel Fortran Compiler"
		junk=`$FC -v 2>&1 | grep -i version | sed "s|[ifort,version,Version, ]||g"`
		intelversion=`echo $junk | awk -F. '{print $1}'`
		#if [ $intelversion -le 17 ] ; then # old syntax, supported up to version 17 included
		if [ $intelversion -lt 15 ] ; then # old syntax, supported up to version 17 included
		    FFLAGS="$IFCINC -cm -w -sox -vec_report0"
		    PRFLAGS="-openmp -openmp_report0" # Open MP enabled # June 2007
		else # new syntax, supported since version 15, required in 18
		    FFLAGS="$IFCINC -cm -w -sox -qopt-report=0"
		    PRFLAGS="-qopenmp" # Open MP enabled # Sept 2017
		fi
		MOD="$IFCMOD"
		FTYPE="$IFCVERSION"
  		OFLAGS="-O3"
		FI8FLAG="-i8" # change default INTEGER to 64 bits
##		FI8FLAG="-integer-size 64" # change default INTEGER to 64 bits
		CFLAGS="$CFLAGS -DINTEL_COMPILER" # to combine C and F90
		MODDIR="-module " # output location of modules
		[ $OS = "Linux" ]  && F90_WLRPATH="-Wl,-R"
		[ $OS = "Darwin" ] && F90_WLRPATH="-Wl,-rpath,"
		DO_F90_SHARED=1
        elif [ $npgf != 0 ] ; then
                #FCNAME="Portland Group Compiler"
                FCNAME="PGI Compiler"
  		OFLAGS="-O3"
		FI8FLAG="-i8" # change default INTEGER and LOGICAL to 64 bits
		PRFLAGS="-mp" # Open MP enabled, to be tested
		MODDIR="-module " # output location of modules
		CFLAGS="$CFLAGS -DpgiFortran" # to combine C and F90
		[ $OS = "Linux" ]  && F90_WLRPATH="-Wl,-R"
		[ $OS = "Darwin" ] && F90_WLRPATH="-Wl,-rpath,"
		DO_F90_SHARED=1
        elif [ $nlah != 0 ] ; then
                FCNAME="Lahey/Fujitsu Compiler"
#  		FFLAGS="$FFLAGS --nap --nchk --npca --ntrace --tpp --trap dio"
		FFLAGS="$FFLAGS --nap --nchk --npca --ntrace --tpp --trap" # (on trial version)
        elif [ $nfuj != 0 ] ; then
                FCNAME="Fujitsu Compiler"
  		FFLAGS="$FFLAGS -Am -X9 -static"
	elif [ $nxlf != 0 ] ; then
	    FTYPE="xlf"
	    if [ "$OS" = "AIX" ] ; then
		FC="xlf90_r"
		FCNAME="IBM XL Fortran"
		FFLAGS="$FFLAGS -qsuffix=f=f90:cpp=F90"
		OFLAGS="-O"
		CC="${CC-gcc}"   # gcc unless already defined
		CFLAGS="$CFLAGS -DRS6000" # to combine C and F90
		FPP="-WF,-D"
		PRFLAGS="-qsmp=omp" # Open MP enabled
		F90_AR="ar -rsv" # archive with index table
		FF64="-q64"
		CF64="-q64"
		AR64="-X64"
	    else
		FC="xlf90"
		FCNAME="IBM XL Fortran for Mac OS"
		FFLAGS="$FFLAGS -qfree=f90 -qsuffix=f=f90:cpp=F90"
		OFLAGS="-O"
		CC="${CC-gcc}"   # gcc unless already defined
		CFLAGS="$CFLAGS -DRS6000" # to combine C and F90
		#### FPP="-WF,-D"
		PRFLAGS="-qsmp=omp" # Open MP enabled
	    fi
	    MODDIR="-qmoddir " # output location of modules
	    FI8FLAG="-qintsize=8" # change default INTEGER to 64 bits
	elif [ $nabs != 0 ] ; then
	        FCNAME="Absoft Pro Compiler"
		FFLAGS=`echo $FFLAGS | ${SED} "s/-I/-p/g"`
		FFLAGS="$FFLAGS -YEXT_NAMES=LCS -YEXT_SFX=_ -q"
		OFLAGS="-O3 -cpu:host"
		F90_LDFLAGS="$F90_LDFLAGS -lU77"
		CFLAGS="$CFLAGS -DAbsoftProFortran"  # to combine C and F90
		CC="${CC-gcc}"   # gcc unless already defined
		MODDIR="-YMOD_OUT_DIR=" # output location of modules
	elif [ $ng95 != 0 ] ; then
	        FCNAME="g95 compiler"
		FFLAGS="$FFLAGS -DGFORTRAN -DG95 -w -ffree-form -fno-second-underscore"
		OFLAGS="-O3"
		CC="${CC-gcc}"   # gcc unless already defined
		CFLAGS="$CFLAGS -DgFortran" # to combine C and F90
		FI8FLAG="-i8" # change default INTEGER to 64 bits
		[ $OS = "Linux" ]  && F90_WLRPATH="-Wl,-R"
		[ $OS = "Darwin" ] && F90_WLRPATH="-Wl,-rpath,"
		MODDIR="-fmod=" # output location of modules
		DO_F90_SHARED=1
	elif [ $ngfortran != 0 ] ; then
	        FCNAME="gfortran compiler"
		OFLAGS="-O3"
		PRFLAGS="-fopenmp" # Open MP enabled
		CC="${CC-gcc}"   # gcc unless already defined
		FI8FLAG="-fdefault-integer-8" # change default INTEGER to 64 bits
		[ $OS = "Linux" ]  && F90_WLRPATH="-Wl,-R"
		[ $OS = "Darwin" ] && F90_WLRPATH="-Wl,-rpath,"
		MODDIR="-J" # output location of modules
		DO_F90_SHARED=1
	elif [ $nflang != 0 ] ; then
	        FCNAME="flang compiler"
		OFLAGS="-O3"
		PRFLAGS="-fopenmp" # Open MP enabled
		CC="${CC-gcc}"   # gcc unless already defined
		FI8FLAG="-fdefault-integer-8" # change default INTEGER to 64 bits
		[ $OS = "Linux" ]  && F90_WLRPATH="-Wl,-R"
		[ $OS = "Darwin" ] && F90_WLRPATH="-Wl,-rpath,"
		MODDIR="-J" # output location of modules
		DO_F90_SHARED=1
	elif [ $npath != 0 ] ; then
	        FCNAME="PathScale EKOPath compiler"
		FFLAGS="$FFLAGS"
		OFLAGS="-O"
		#CC="pathcc"
		CC="${CC-pathcc}"   # pathcc unless already defined
		PRFLAGS="-mp" # Open MP enabled
		FI8FLAG="-i8" # change default INTEGER to 64 bits
		#FI8FLAG="-default64" # change default INTEGER and FLOAT to 64 bits
		MODDIR="-module " # output location of modules
        else
	    nvas=`$FC | ${GREP} -i sierra | ${WC} -l`
            if [ $nvas != 0 ] ; then
                FCNAME="Pacific/Sierra Compiler"
		echo "$FCNAME is not supported"
		crashAndBurn
	    else
                echo "$FC: Unknown compiler"
                crashAndBurn
	    fi
        fi
}

# -----------------------------------------------------------------
add64bitF90Flags () {

    if [ "x$FF64$CF64$AR64" != "x" ]; then
	echo "Do you want to make a 64 bit compilation ? [y/N]"
	read answer
	if [ "x$answer" = "xy" -o "x$answer" = "xY" ]; then
	    FFLAGS="$FFLAGS $FF64"
	    CFLAGS="$CFLAGS $CF64"
	    F90_AR="$F90_AR $AR64"
	fi
    fi
}
# -----------------------------------------------------------------
countF90Bits () {
    # count bits in F90 compilation
    tmpfile=to_be_removed
    suffix=.f90
    ${CAT} > ${tmpfile}${suffix} <<EOF
program test
end program test
EOF
    FFLAGS_=`par2curly  ${FFLAGS}`  # $() -> ${}
    FFLAGS_=`eval echo ${FFLAGS_}` # ${} -> value
    ${FC} ${FFLAGS_} ${tmpfile}${suffix} -o ${tmpfile} 1>${DEVNULL} 2>&1
    f90_ok=0
    [ -s ${tmpfile} ] && f90_ok=1
    f90_64=`${FILE} ${tmpfile} | ${GREP} 64 | ${WC} -l`
    ${RM}  ${tmpfile}*
}
# -----------------------------------------------------------------
countCBits () {
    # count bits in C compilation
    tmpfile=to_be_removed
    suffix=.c
    ${CAT} > ${tmpfile}${suffix} <<EOF
int main(){
}
EOF

    $CC $CFLAGS ${tmpfile}${suffix} -o ${tmpfile} 1>${DEVNULL} 2>&1
    c_ok=0
    [ -s ${tmpfile} ] && c_ok=1
    c_64=`${FILE} ${tmpfile} | ${GREP} 64 | ${WC} -l`
    ${RM}  ${tmpfile}*
}
# -----------------------------------------------------------------
checkF90Compilation () {
    # check that F90 compiler actually works
    # requires compilation and execution of F90 code
    tmpfile=./to_be_removed
    suffix=.f90
    ${CAT} > ${tmpfile}${suffix} <<EOF
program test
    print*,'hello'
end program test
EOF
    canrun=0
    cancompile=0
    FFLAGS_=`par2curly  ${FFLAGS}`  # $() -> ${}
    FFLAGS_=`eval echo ${FFLAGS_}` # ${} -> value
    ${FC} ${FFLAGS_} ${tmpfile}${suffix} -o ${tmpfile}  1>${DEVNULL} 2>&1
    [ -s ${tmpfile} ] && cancompile=1
    if [ -x ${tmpfile} ] ; then
	canrun=`${tmpfile} | grep hello | ${WC} -l`
    fi
    ${RM} ${tmpfile}*

    if [ $cancompile -eq 0 ]; then
	echo
	echo "  ERROR: Compilation with "
	echo "${FC} ${FFLAGS}"
	echo "currently fails."
	echo
	echo "Please check that this compiler is supported by your system"
	crashAndBurn
    fi
    if [ $canrun -eq 0 ]; then
	echo
	echo "  WARNING: Currently the codes compiled with "
	echo "${FC} ${FFLAGS}"
	echo "can not be executed."
	echo "Most likely, some compiler related dynamic libraries are not found."
	echo "(Check the LD_LIBRAY_PATH variable and read the compiler documentation.)"
	echo
	echo "That can affect negatively the result of this configuration script."
	echo
	echo
    fi
}
# -----------------------------------------------------------------
checkF90LongLong () {
    # check that F90 support 8 byte integers
    # requires compilation and execution of F90 code
    tmpfile=./to_be_removed
    suffix=.f90

    # conduct real test
    ${CAT} > ${tmpfile}${suffix} <<EOF
program test
    if (selected_int_kind(16) > selected_int_kind(9)) print*,'OK'
end program test
EOF
    longlong=0
    FFLAGS_=`par2curly  ${FFLAGS}`  # $() -> ${}
    FFLAGS_=`eval echo ${FFLAGS_}` # ${} -> value
    ${FC} ${FFLAGS_} ${tmpfile}${suffix} -o ${tmpfile}  1>${DEVNULL} 2>&1
    if [ -x ${tmpfile} ] ; then
	longlong=`${tmpfile} | grep OK | ${WC} -l`
    fi
    ${RM} ${tmpfile}*

}
# -----------------------------------------------------------------

askUserF90 () {
    echoLn "enter name of your F90 compiler [$FC]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer    
    [ "x$answer" != "x" ] && FC="$answer"
}

# -----------------------------------------------------------------

showDefaultDirs () {
    echo " compiled Healpix products will be:"
    echo "F90_BINDIR =  ${F90_BINDIR}[suffix]"
    echo "F90_INCDIR =  ${F90_INCDIR}[suffix]"
    echo "F90_LIBDIR =  ${F90_LIBDIR}[suffix]"
    echo "F90_BUILDDIR =  ${F90_BUILDDIR}[suffix]"
#    echo " and the Makefile will be copied into Makefile[suffix]"
}

updateDirs () {
    F90_BINDIR=${F90_BINDIR}${F_DIRSUFF}
    F90_INCDIR=${F90_INCDIR}${F_DIRSUFF}
    F90_LIBDIR=${F90_LIBDIR}${F_DIRSUFF}
    F90_BUILDDIR=${F90_BUILDDIR}${F_DIRSUFF}
}

showActualDirs () {
    echo " compiled Healpix products will be:"
    echo "F90_BINDIR =  ${F90_BINDIR}"
    echo "F90_INCDIR =  ${F90_INCDIR}"
    echo "F90_LIBDIR =  ${F90_LIBDIR}"
    echo "F90_BUILDDIR =  ${F90_BUILDDIR}"
}
# -----------------------------------------------------------------
askUserMisc () {
    echo "  Note: your Fortran compiler is $FCNAME"

    echoLn " "

    add64bitF90Flags

    showDefaultDirs
    echoLn "enter suffix for directories [$F_DIRSUFF]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && F_DIRSUFF="$answer"
    updateDirs
    showActualDirs

    checkDir $F90_BINDIR $F90_INCDIR $F90_LIBDIR $F90_BUILDDIR
    fullPath F90_BINDIR  F90_INCDIR  F90_LIBDIR  F90_BUILDDIR

    echoLn " "

    echoLn "enter compilation flags for $FC compiler [$FFLAGS]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && FFLAGS="$answer"

    echoLn "enter optimisation flags for $FC compiler [$OFLAGS]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && OFLAGS="$answer"

    checkF90Compilation

    checkF90LongLong
    if [ ${longlong} = 0 ] ; then
	echo "Your compiler does not seem to support 8-byte integers"
	echo "The compilation flag ${FPP}NO64BITS will be added to prevent their usage."
	FFLAGS="${FFLAGS} ${FPP}NO64BITS"
    fi

    FFLAGS="$OFLAGS $FFLAGS"
    echo "  Fortran code will be compiled with $FC $FFLAGS"


    echoLn "enter name of your C compiler [$CC]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CC="$answer"
    IdentifyCCompiler

    echoLn "enter compilation/optimisation flags for C compiler [$CFLAGS]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && CFLAGS="$answer"

    countF90Bits
    countCBits
    if [ ${f90_ok} -eq 0 -o ${c_ok} -eq 0 ] ; then
	echo
	echo "ERROR!:"
	if [ ${f90_ok} -eq 0 ] ; then
	    echo " F90 Compilation ($FC $FFLAGS) failed!"
	fi
	if [ ${c_ok} -eq 0 ] ; then
	    echo " C Compilation ($CC $CFLAGS) failed!"
	fi
	echo "Check the availability of the compiler(s) and the options used."
	crashAndBurn
    fi
    if [ $c_64 != $f90_64 ] ; then
	echo "Warning: "
	if [ $f90_64 != 0 ] ; then
	    echo "F90 compiler generates 64 bit code, "
	    echo "while C compiler generates 32 bit code"
	else
	    echo "F90 compiler generates 32 bit code, "
	    echo "while C compiler generates 64 bit code"
	fi
	echoLn "you may want to change the C compilation options [$CFLAGS]: "
	read answer
	[ "x$answer" != "x" ] && CFLAGS="$answer"
	echoLn "or the F90 compilations options ($FFLAGS): "
	read answer
	[ "x$answer" != "x" ] && FFLAGS="$answer"
	echo "you also may have to recompile cfitsio with the correct options to ensure that its C routines and Fortran wrappers are consistent with each other and with Healpix"
    fi
    echo "  C subroutines will be compiled with $CC $CFLAGS"

    echoLn "enter command for library archiving [$F90_AR]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && F90_AR="$answer"

    echoLn "enter full name of cfitsio library [lib${LIBFITS}.a]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && LIBFITS=`${BASENAME} $answer ".a" | ${SED} "s/^lib//"`

    findFITSLib $LIBDIR $FITSDIR
    echoLn "enter location of cfitsio library [$FITSDIR]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    [ "x$answer" != "x" ] && FITSDIR=$answer
    fullPath FITSDIR

    lib="${FITSDIR}/lib${LIBFITS}.a"
    if [ ! -r $lib ]; then
	echo
	echo "error: fits library $lib not found"
	echo
	crashAndBurn
    fi

    # add option on where to search runtime libraries, on compilers supporting it
    if [ "x$F90_WLRPATH" != "x" ] ; then
	F90_WLRPATH_="${F90_WLRPATH}${FITSDIR} ${F90_WLRPATH}${SHARPLDIR} ${F90_WLRPATH}${F90_LIBDIR}" # expand $FITSDIR,$SHARPLDIR,F90_LIBDIR
	F90_WLRPATH="${F90_WLRPATH}\$(FITSDIR) ${F90_WLRPATH}\$(SHARPLDIR) ${F90_WLRPATH}\$(F90_LIBDIR)"  # keep $(FITSDIR), ...
	F90_LDFLAGS="${F90_LDFLAGS} ${F90_WLRPATH}"
    fi

    checkF90Fitsio  ${lib}
    checkFitsioCurl ${lib}
    F90_LDFLAGS="$F90_LDFLAGS $CFITSIOCURL"
    checkF90FitsioLink
    checkF90FitsioVersion

}


# -----------------------------------------------------------------
askPgplot () {

    echo "  The generator of non-gaussian CMB maps (ng_sims) can optionally "
    echo "produce plots of the maps Prob. Dens. Function using PGPLOT."
    echo "Do you want to enable this option ? "
#    echo " (the ng_sims code will run just as well without it,)  "
    echoLn "(this assumes that PGPLOT is already installed on your computer) (y|n) [n]: "
    read answer
    if [ "x$answer" = "xy"  -o "x$answer" = "xY" ]; then
	PGFLAG="${FPP}PGPLOT" # set preprocessing variable
	PGLIBS=${PGLIBSDEF}
	echo " Enter the options necessary to link to PGPLOT:"
	echoLn "( $PGLIBSDEF ) :"
	read answer
	[ "x$answer" != "x" ] && PGLIBS="$answer"
    fi
    echo " "

}

# -----------------------------------------------------------------

editF90Makefile () {

    echoLn "Editing top Makefile for F90 ..."
#    [ -r Makefile ] && mv Makefile Makefile.bak

    F90_BINDIR_H=`echo ${F90_BINDIR} | sed "s|${HEALPIX}|\\$(HEALPIX)|g"`
    F90_INCDIR_H=`echo ${F90_INCDIR} | sed "s|${HEALPIX}|\\$(HEALPIX)|g"`
    F90_LIBDIR_H=`echo ${F90_LIBDIR} | sed "s|${HEALPIX}|\\$(HEALPIX)|g"`
    F90_BUILDDIR_H=`echo ${F90_BUILDDIR} | sed "s|${HEALPIX}|\\$(HEALPIX)|g"`
    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^F90_FC.*$|F90_FC	= $FC|" |\
	${SED} "s|^F90_FFLAGS.*$|F90_FFLAGS	= $FFLAGS|" |\
	${SED} "s|^F90_LDFLAGS.*$|F90_LDFLAGS	= $F90_LDFLAGS|" |\
	${SED} "s|^F90_CC.*$|F90_CC	= $CC|" |\
	${SED} "s|^F90_CFLAGS.*$|F90_CFLAGS	= $CFLAGS|" |\
	${SED} "s|^HEALPIX=.*$|HEALPIX	= $HEALPIX|" |\
	${SED} "s|^FITSDIR.*$|FITSDIR	= $FITSDIR|" |\
	${SED} "s|^LIBFITS.*$|LIBFITS	= $LIBFITS|" |\
	${SED} "s|^F90_BINDIR.*$|F90_BINDIR	= $F90_BINDIR_H|" |\
	${SED} "s|^F90_INCDIR.*$|F90_INCDIR	= $F90_INCDIR_H|" |\
	${SED} "s|^F90_LIBDIR.*$|F90_LIBDIR	= $F90_LIBDIR_H|" |\
	${SED} "s|^F90_BUILDDIR.*$|F90_BUILDDIR	= $F90_BUILDDIR_H|" |\
	${SED} "s|^F90_AR.*$|F90_AR        = $F90_AR|" |\
	${SED} "s|^F90_MODDIR[[:space:]=].*$|F90_MODDIR	= \"$MODDIR\"|" |\
	${SED} "s|^F90_MOD[[:space:]=].*$|F90_MOD	= $MOD|" |\
	${SED} "s|^F90_PPFLAGS.*$|F90_PPFLAGS	= $PPFLAGS|" |\
	${SED} "s|^F90_PGFLAG.*$|F90_PGFLAG  = $PGFLAG|" |\
	${SED} "s|^F90_PGLIBS.*$|F90_PGLIBS  = $PGLIBS|" |\
	${SED} "s|^F90_OS.*$|F90_OS	= $OS|" |\
	${SED} "s|^F90_I8FLAG.*$|F90_I8FLAG  = $FI8FLAG|" |\
	${SED} "s|^F90_LIBSUFFIX.*$|F90_LIBSUFFIX = $F90_LIBSUFFIX|" |\
	${SED} "s|^F90_FLAGNAMELIB.*$|F90_FLAGNAMELIB = $F90_FLAGNAMELIB|" |\
	${SED} "s|^ALL\(.*\) f90-void\(.*\)|ALL\1 f90-all\2|" |\
	${SED} "s|^TESTS\(.*\) f90-void\(.*\)|TESTS\1 f90-test\2|" |\
	${SED} "s|^CLEAN\(.*\) f90-void\(.*\)|CLEAN\1 f90-clean\2|" |\
	${SED} "s|^DISTCLEAN\(.*\) f90-void\(.*\)|DISTCLEAN\1 f90-distclean\2|" |\
	${SED} "s|^TIDY\(.*\) f90-void\(.*\)|TIDY\1 f90-tidy\2|" > Makefile


# 	if [ "x$F_DIRSUFF" != "x" ] ; then
# 	    if [ "x$F_DIRSUFF" != "x.in" ] ; then
# 		${CP} Makefile Makefile$F_DIRSUFF
# 	    fi
# 	fi

    echo " done."
    edited_makefile=1

}

# -----------------------------------------------------------------

generateConfF90File () {
	echo "Generating $HPX_CONF_F90"

    	echo "# F90 configuration for HEALPix `date`" > $HPX_CONF_F90
# put HEALPIX variable back into F90_BINDIR
    dollar="$"
    F90_BINDIR_SHORT=`echo $F90_BINDIR | sed "s|$HEALPIX|{HEALPIX}|g"`
    F90_BINDIR_SHORT="${dollar}${F90_BINDIR_SHORT}"

    case $SHELL in
    csh|tcsh)
	${CAT} <<EOF >>$HPX_CONF_F90
setenv HEXE    ${F90_BINDIR_SHORT}
setenv PATH    \${HEXE}:\${PATH}
EOF
	;;
    sh|ksh|bash|zsh)
	${CAT} <<EOF >>$HPX_CONF_F90
HEXE=${F90_BINDIR_SHORT}
PATH="\${HEXE}:\${PATH}"
export HEXE PATH
EOF
	;;
    *)
	echo "Shell $SHELL not supported yet."
	${RM}  $HPX_CONF_F90
	;;
    esac
}

# -----------------------------------------------------------------

offerF90Compilation () {
    echo "F90 Configuration finished."
    echo "You can run \"(GNU)make\" to build the package,"
    echo "        and \"(GNU)make test\" to test it."
    echoLn "You can also choose to build the package right now from here (Y|n): "
    read answer
    if [ "x$answer" != "xn" -a  "x$answer" != "xN" ]; then
    # find out make command
	askMake
    # make compilation
	${MAKE}   || crashAndBurn
	${MAKE} test || crashAndBurn
    #
	echo
	echo
	echo "F90 package installed !"
	echo
    fi
}

# -----------------------------------------------------------------
f90_shared () {
    echo " "
    echo " Experimental feature: "
    echoLn "A static library is produced by default. Do you rather want a shared/dynamic library ?"
    if [ `isTrue ${F_SHARED}` -eq 1 ]; then
	echoLn " (y|n) [y]: "
	read answer
	[ `isFalse ${answer}` -eq 1 ] && DO_F90_SHARED=0  || DO_F90_SHARED=1
    else
	echoLn " (y|n) [n]: "
	read answer
	[ `isTrue ${answer}` -eq 1 ] && DO_F90_SHARED=1  || DO_F90_SHARED=0
    fi
    [ $INTERACTIVE -eq 0 ] && echo $answer
#     echo "============================"
#     echo "F_SHARED =  ${F_SHARED}"
#     echo "answer=${answer}."  
#     echo "DO_F90_SHARED =  ${DO_F90_SHARED}"
#     #echo `isTrue ${answer}` `isFalse ${answer}`
#     echo "============================"
    if [ ${DO_F90_SHARED} -eq 1 ]; then
	case $OS in
	    Darwin)
		F90_LIBSUFFIX=".dylib"
		if [ `isTrue ${USE_ATRPATH}` -eq 1 ]; then
		    F90_AR="${FC} ${F90PIC} -dynamiclib -Wl,-undefined,dynamic_lookup -Wl,-rpath,\$(F90_LIBDIR) -o "
		    F90_FLAGNAMELIB="-Wl,-install_name,@rpath/"
		else
		    F90_AR="${FC} ${F90PIC} -dynamiclib -Wl,-undefined,dynamic_lookup -o "
		    F90_FLAGNAMELIB=""
		fi
		;;
	    Linux)
		F90_AR="${FC} ${F90PIC} -shared -o "
		F90_FLAGNAMELIB="-Wl,-soname,"
		F90_LIBSUFFIX=".so";;
	    *)
		F90_AR="${FC} ${F90PIC} -shared -o "
		F90_FLAGNAMELIB="-Wl,-soname,"
		F90_LIBSUFFIX=".so";;
	esac
    fi
}
# -----------------------------------------------------------------
writeF90pkgconfigFile (){
    #pkgconfigFile=${HEALPIX}/lib${F_DIRSUFF}/healpix.pc
    pkgconfigFile=${HEALPIX}/lib${F_DIRSUFF}/pkgconfig/healpix.pc
    ${MKDIR} -p `${DIRNAME} ${pkgconfigFile}`
    echo
    echo "Writing pkgconfig file: ${pkgconfigFile}"
    ${CAT}<<EOF > ${pkgconfigFile}
# HEALPix/F90 pkg-config file
# compiled with ${FC}

compiler=${FC}
prefix=${HEALPIX}
suffix=${F_DIRSUFF}
exec_prefix=\${prefix}/bin\${suffix}
libdir=\${prefix}/lib\${suffix}
sharpdir=\${prefix}/lib
includedir=\${prefix}/include\${suffix}

Name: HEALPix
Description: F90 library for HEALPix (Hierarchical Equal-Area iso-Latitude) pixelisation of the sphere
Version: ${HPXVERSION}
URL: https://healpix.sourceforge.io
Requires: cfitsio >= ${CFITSIOVREQ}
Libs: -L\${libdir} -L\${sharpdir} -lhealpix -lhpxgif -lsharp
Cflags: -I\${includedir} ${PRFLAGS} ${F90PIC}

EOF

    ${CAT} ${pkgconfigFile}
    echo "           -------------- "
    echo

}
# -----------------------------------------------------------------

f90_config () {
    test_Sharp
    if [ $sharp_configured = "0" ] ; then
	echo "Configuring the libsharp library first, since this is a dependency:"
	Sharp_config
	Sharp_install
	echo
	echo "Now configuring Healpix F90 itself:"
    fi
    HPX_CONF_F90=$1
    setF90Defaults
    askUserF90
    GuessF90Compiler
    askUserMisc
    askPgplot
    askOpenMP
    askF90PIC
    patchF90
    f90_shared
    #makeProfile
    generateConfF90File
    editF90Makefile
    [ $NOPROFILEYET = 1 ] && installProfile
    writeF90pkgconfigFile
#    offerF90Compilation
}


#=====================================
#=========== Check Configuration ===========
#=====================================
checkConfFiles () {

    echo "Currently, the configuration files created are :"
    echo "__________________________________________________________________"
    for conffile in ${HPX_CONF_DIR}/*; do
	echo "${conffile} : "
	${CAT} ${conffile}
	echo
    done
    echo "__________________________________________________________________"

}

#=====================================
#=========== Top package ===========
#=====================================
#   mainMenu:
#   installProfile: modify user's shell profile if agreed
#   makeTopConf:
#   readyTopMakefile:
#   setTopDefaults:
#   setConfDir:
#-------------
mainMenu () {

    echo
    echo "Do you want to:"
    echo "(0): exit"
    echo "(1): configure Healpix IDL package"
    echo "(2): configure Healpix C   package, and edit Makefile"
    echo "(3): configure Healpix F90 package, and edit Makefile"
    echo "(4): configure Healpix C++ package, and edit Makefile"
    echo "(5): configure Healpix Python (healpy) package, and edit Makefile"
    echo "(7): configure, compile and install libsharp package, and edit Makefile"
    echo "(8): see what configuration files have been created so far"
    echo "(9): edit your shell configuration file to have easier access to Healpix codes"
    echo "(-1): reset"
    echo "     (will *REMOVE* the Makefile and configuration files, and exit)"
    echo "(0): exit"
    echo
    echoLn "Enter your choice (configuration of packages can be done in any order) [0]: "
    read answer
    [ $INTERACTIVE -eq 0 ] && echo $answer
    case x$answer in
	x1)
	  eval idlconffile=$HPX_CONF_IDL
	  eval gdlconffile=$HPX_CONF_GDL
	  eval flconffile=$HPX_CONF_FL
          idl_config $idlconffile $gdlconffile $flconffile;;
        x2)
           C_config;;
	x3)
	   eval f90conffile=$HPX_CONF_F90
	   f90_config $f90conffile;;
	x4)
	   eval cppconffile=$HPX_CONF_CPP
	   Cpp_config $cppconffile;;
 	x5)
 	   Healpy_config;;
	x7)
	  Sharp_config
	  Sharp_install;;
	x8)
	   checkConfFiles;;
	x9)
	   installProfile;;
	x0)
	   goodBye;;
	x)
	   goodBye;;
	"x-1")
	   restartFromScratch
	   goodBye;;
	x*)
	   echo "unknown answer !"
	   echo
	   crashAndBurn;;
    esac
}
#-------------
installProfile () {
    # will modity user's configuration file to invoke Healpix configuration

    case $SHELL in
    sh|ksh|bash|zsh)
	prof="${HOME}/.profile"
	comd="[ -r ${HPX_CONF_MAIN} ] && . ${HPX_CONF_MAIN}";;
    csh)
	prof="${HOME}/.cshrc"
	comd="if ( -e ${HPX_CONF_MAIN} ) source ${HPX_CONF_MAIN}";;
    tcsh)
	prof="${HOME}/.tcshrc"
	[ ! -r $prof -a -r "${HOME}/.cshrc" ] && prof="${HOME}/.cshrc"
	comd="if ( -e ${HPX_CONF_MAIN} ) source ${HPX_CONF_MAIN}";;
    *) ;;
    esac
    [ ! -r $prof ] && touch $prof
    # do not do edition if it was previously done
    check=`${GREP} ${HPX_CONF_MAIN} $prof | ${WC} -l`
    if [ $check -eq 0 ]; then
	${CAT} <<EOF

The following line should be inserted into your home shell profile ($prof):

  $comd

 Where the file ${HPX_CONF_MAIN} contains:
EOF
${CAT} ${HPX_CONF_MAIN}

	echo ""
	if [ `isFalse ${PROFILE_EDIT}` -eq 1 ] ; then
	    echoLn "Do you want this modification to be done (y|n)? [n]: "
	    read answer
	    [ `isTrue ${answer}` -eq 1 ] && edit_prof=1 || edit_prof=0
	else
	    echoLn "Do you want this modification to be done (y|n)? [y]: "
	    read answer
	    [ `isFalse ${answer}` -eq 1 ] && edit_prof=0 || edit_prof=1
        fi
	[ $INTERACTIVE -eq 0 ] && echo $answer
	if [ ${edit_prof} -eq 1 ] ; then
	    ${CP} $prof ${prof}".save"
	    echo "" >> $prof
	    echo "# modifications by HEALPixAutoConf ${HPXVERSION}" >> $prof
	    echo $comd >> $prof
	    echo "Modification done and previous shell profile saved"
	    echo "as ${prof}.save."
	fi
    else
	echo "Your home shell profile ($prof)"
	echo "has already been edited."
	if [ ${INTERACTIVE} -eq 0 ] ; then
	    read answer  # just to keep number of questions unchanged
	fi
    fi
    NOPROFILEYET=0
}
#-------------
makeTopConf(){

    ${MKDIR} -p ${HPX_CONF_DIR}

    case $SHELL in
    sh|ksh|bash|zsh)
	${CAT}<<EOF >| ${HPX_CONF_MAIN}
# configuration for Healpix $HPXVERSION
HEALPIX=${HEALPIX} ; export HEALPIX
HPX_CONF_DIR=${HPX_CONF_DIR}
if [ -r ${HPX_CONF_IDL} ] ; then . ${HPX_CONF_IDL} ; fi
if [ -r ${HPX_CONF_GDL} ] ; then . ${HPX_CONF_GDL} ; fi
if [ -r ${HPX_CONF_FL} ]  ; then . ${HPX_CONF_FL}  ; fi
if [ -r ${HPX_CONF_F90} ] ; then . ${HPX_CONF_F90} ; fi
if [ -r ${HPX_CONF_CPP} ] ; then . ${HPX_CONF_CPP} ; fi
if [ -r ${HPX_CONF_C} ] ;   then . ${HPX_CONF_C} ;   fi
EOF
    echo ' ' ;;
    csh|tcsh)
	${CAT}<<EOF >| ${HPX_CONF_MAIN}
# configuration for Healpix $HPXVERSION
setenv HEALPIX $HEALPIX
setenv HPX_CONF_DIR ${HPX_CONF_DIR}
if ( -e ${HPX_CONF_IDL} ) source ${HPX_CONF_IDL}
if ( -e ${HPX_CONF_GDL} ) source ${HPX_CONF_GDL}
if ( -e ${HPX_CONF_FL}  ) source ${HPX_CONF_FL}
if ( -e ${HPX_CONF_F90} ) source ${HPX_CONF_F90}
if ( -e ${HPX_CONF_CPP} ) source ${HPX_CONF_CPP}
if ( -e ${HPX_CONF_C} )   source ${HPX_CONF_C}
EOF
    echo ' '  ;;
    *) ;;
    esac

}
#-------------
readyTopMakefile () {

    # backup name
    sdate=`date +%s`
    [ "x${sdate}" = "x" ] && sdate="1"
    mkbk=Makefile_bk${sdate}
    [ -s ${mkbk} ] && mkbk="${mkbk}a"

    if [ -s Makefile ] ; then
	${CP} -f Makefile ${mkbk}
    else
	if [ ! -r Makefile.in ] ; then
	    echo "top makefile template (Makefile.in) was not found. Can not proceed."
	    crashAndBurn
	fi
	${CP} -f Makefile.in Makefile
    fi

}
#-------------
restartFromScratch () {

    echo "Removing Main Makefile"
    ${RM} Makefile
    echo "Removing configuration files in " ${HPX_CONF_DIR}
    for hfile in ${HPX_CONF_MAIN} ${HPX_CONF_IDL} ${HPX_CONF_GDL} ${HPX_CONF_FL} ${HPX_CONF_F90} ${HPX_CONF_CPP} ${HPX_CONF_C} ; do
	eval thisfile=${hfile}
        ${RM} ${thisfile}
    done
    echo "Removing configuration directory: " ${HPX_CONF_DIR}
    ${RMDIR} ${HPX_CONF_DIR}

}
#-------------
setTopDefaults() {

    AWK="awk"
    BASENAME="basename"
    CAT="cat"
    CC="${CC-gcc}"   # gcc unless already defined
    CXX="${CXX-g++}" # g++ unless already defined
    CP="cp"
    DEVNULL="/dev/null"
    DIRNAME="dirname"
    FILE="file"
    GREP="grep"
    HEAD="head" # introduced 2008-11-21
    LS="ls"
    MAKE="${MAKE-make}" # make unless already defined
    MKDIR="mkdir"
    NM="nm"
    PRINTF="printf"
    #PWD="pwd"
    PYTHON="${PYTHON-python}" # python unless already defined
    RM="/bin/rm -f"
    RMDIR="rmdir"
    SED="sed"
    WC="wc"
    OS=`uname -s`
    WHEREIS="whereis"

    HEALPIX=`pwd`
    #echo "HEALPIX ${HEALPIX}"
    SHARPPREFIX=$HEALPIX
    SHARPLDIR=${SHARPPREFIX}/lib
    CXXPREFIX=$HEALPIX
    CXXDIR=${HEALPIX}/src/cxx

    USE_ATRPATH="${USE_ATRPATH-0}" # if set (on MacOS), use @rpath for location of libsharp.dylib, libhealpix_cxx.dylib, libhealpix.dylib
    NOPROFILEYET=1
    SHELL=`${BASENAME} ${SHELL-/bin/sh}`

    MAKESET=0

    F_SHARED="${F_SHARED-0}" # 0 unless already defined
    F_PARAL="${F_PARAL-1}"   # 1 unless already defined
    LIBFITS="cfitsio"
    # look for libcfitsio.* in various directories starting with /usr/local/lib or user defined value
    findFITSLib     ${FITSDIR-/usr/local/lib}
    # look for fitsio.h in various directories starting with /usr/local/lib or user defined value
    findFITSInclude ${FITSINC-/usr/local/include}
    # define FITSDIR and FITSINC with /usr/local/lib and /usr/local/include as initial guess
    # findFITSLib /usr/local/lib
    # findFITSInclude /usr/local/include `cd ${FITSDIR}/../include ; pwd`
    # FITSDIR="${FITSDIR-/usr/local/lib}"     # /usr/local/lib     unless already defined
    # FITSINC="${FITSINC-/usr/local/include}" # /usr/local/include unless already defined
    FITSPREFIX="/usr/local"

    C_SHARED=${C_SHARED-0} # 0 unless already defined
    C_FITS="${C_FITS-1}" # 1 unless already defined
    #C_WITHOUT_CFITSIO="${C_WITHOUT_CFITSIO-0}" # 0 unless already defined

    PROFILE_EDIT="${PROFILE_EDIT-n}" # n unless already defined

    edited_makefile=0


    HPX_VERSION=`echo $HPXVERSION | ${SED} "s|\.|_|g"`
    HPX_CONF_DIR_HOME=${HOME}/.healpix/${HPX_VERSION}_${OS}
    HPX_CONF_DIR_INPLACE=${HEALPIX}/confdir/${HPX_VERSION}_${OS}

}
#-------------
setConfDir () {

    case $SHELL in
    sh|ksh|bash|zsh)
        suffix=sh;;
    csh|tcsh)
	suffix=csh;;
    *) ;;
    esac

    HPX_CONF_MAIN=$HPX_CONF_DIR/config
    HPX_CONF_IDL=\${HPX_CONF_DIR}/idl.${suffix}
    HPX_CONF_GDL=\${HPX_CONF_DIR}/gdl.${suffix}
    HPX_CONF_FL=\${HPX_CONF_DIR}/fl.${suffix}
    HPX_CONF_F90=\${HPX_CONF_DIR}/f90.${suffix}
    HPX_CONF_CPP=\${HPX_CONF_DIR}/cpp.${suffix}
    HPX_CONF_C=\${HPX_CONF_DIR}/c.${suffix}


}

#----------------------------------------
