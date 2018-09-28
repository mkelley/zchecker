ZDIR=/oort/msk
ZDATA=/oort1/ZTF
ZWEB=${ZDATA}/web
zchecker clean-eph $1
zchecker eph $*
zchecker search $*
zchecker cutouts
zproject
zstack
"(cd $ZDATA; sqlite3 <${ZDIR}/zbrowser/scripts/dump-foundobs.sql && mv /tmp/foundobs.db $ZWEB/)"
python3 ${ZDIR}/zbrowser/scripts/stack2web.py ${ZWEB}/stacks --desg=$1 -f
