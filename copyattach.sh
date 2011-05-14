#!/bin/bash
# Attach the copyright and version information at the top of files.
# If $1= -i use individual arguments. If $1 = anything else. Use that.
if [ ".$1" == . ] ; then 
    echo Getting the overall sceptic version from CVS.
    VERSION=`cvs log sceptic.F | grep head: | sed -e "s/head://"`
fi
for file in *.f *.F ; do
    echo -n "Versioning file $file   "
    cp $file $file.bak
    sed -e "/c___c/,/c___c/ d" $file > temp.tmp
    cat >$file <<EOF
c___________________________________________________________________________
c
c     SCEPTIC3D
c
c     This code is copyright (c)
c              Ian H Hutchinson    hutch@psfc.mit.edu.
c              Leonardo Patacchini patacchi@mit.edu
c
c     It may be used freely with the stipulation that any scientific or
c     scholarly publication concerning work that uses the code must give
c     an acknowledgement referring to the relevant papers
c
c     I.H. Hutchinson, Plasma Physics and Controlled Fusion, vol 44, p
c     1953 (2002), vol 45, p 1477 (2003).
c
c     L. Patacchini and I.H. Hutchinson, Plasma Physics and Controlled
c     Fusion, vol 49, p1193 (2007), vol 49, p 1719 (2007).
c
c     I.H. Hutchinson and L. Patacchini, Physics of Plasmas, vol 14,
c     p013505 (2007)
c
c     The code may not be redistributed except in its original package.
c
c     No warranty, explicit or implied, is given. If you choose to build
c     or run the code, you do so at your own risk.
c___________________________________________________________________________

c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c
c     This code is copyright (c)
c              Ian H Hutchinson    hutch@psfc.mit.edu.
c              Leonardo Patacchini patacchi@mit.edu
c
c     It may be used freely with the stipulation that any scientific or
c     scholarly publication concerning work that uses the code must give
c     an acknowledgement referring to the relevant papers
c
c     I.H. Hutchinson, Plasma Physics and Controlled Fusion, vol 44, p
c     1953 (2002), vol 45, p 1477 (2003).
c
c     L. Patacchini and I.H. Hutchinson, Plasma Physics and Controlled
c     Fusion, vol 49, p1193 (2007), vol 49, p 1719 (2007).
c
c     I.H. Hutchinson and L. Patacchini, Physics of Plasmas, vol 14,
c     p013505 (2007)
c
c     The code may not be redistributed except in its original package.
c
c     No warranty, explicit or implied, is given. If you choose to build
c     or run the code, you do so at your own risk.
c
EOF
if [ ".$1" != . ] ; then 
    if [ "$1" == "-i" ] ; then
	VERSION=`cvs log $file | grep head: | sed -e "s/head://"`
    else
	VERSION="$1" 
    fi
fi
echo Version: $VERSION
echo "c     Version: $VERSION   `date`" >> $file
echo c >> $file
echo c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___ >> $file
    cat temp.tmp >> $file
    rm $file.bak
done
rm -f temp.tmp
