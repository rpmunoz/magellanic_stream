
Visit http://archive.eso.org/eso/eso_archive_main.html
Fill "Program ID" field, write 10000 in "Return max" field and click Search button 

Program : 098.B-0510

sed -i '.ORIG' 's/data_with_raw_calibs/raw/g' downloadRequest302668script.sh 

chmod u+x downloadRequest302668script.sh

./downloadRequest302668script.sh

./downloadRequest302668script.sh -d "--no-check-certificate"